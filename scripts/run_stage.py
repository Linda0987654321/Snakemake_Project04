"""Run one pair stage and record failures without blocking other pairs."""
import argparse
import json
import re
import shutil
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

STAGE_DIRS = {"clean_structures": "cleaned", "lightdock_setup": "setup",
              "lightdock_simulate": "simulation", "lightdock_rank": "ranked"}


def require_file(path, protein=False):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty file: {path}")
    if protein and not any(line.startswith("ATOM  ")
                           for line in path.read_text().splitlines()):
        raise ValueError(f"No protein ATOM records: {path}")
    return path


def run_command(args, cwd, log):
    print("Command:", " ".join(map(str, args)), file=log, flush=True)
    subprocess.run(list(map(str, args)), cwd=cwd, stdout=log, stderr=log, check=True)


def relocate(value, old, new):
    # Update only absolute paths inside the copied setup directory.
    if isinstance(value, dict):
        return {k: relocate(v, old, new) for k, v in value.items()}
    if isinstance(value, list):
        return [relocate(v, old, new) for v in value]
    if isinstance(value, str):
        try:
            p = Path(value)
            if p.is_absolute():
                return str(new / p.relative_to(old))
        except ValueError:
            pass
    return value


def perform(stage, pair, root, out, previous, settings, threads,
            receptor, ligand, log, status, execute):
    helper = root / "scripts"
    if stage == "clean_structures":
        # Validate both partners before preparing either one.
        sources = [require_file(p, protein=True) for p in (receptor, ligand)]
        for source, partner in zip(sources, ("rec", "lig")):
            dest = out / f"{pair}_{partner}_clean.pdb"
            execute([sys.executable, helper / "fix_and_clean_pdb.py", source, dest], root, log)
            require_file(dest, protein=True)
    elif stage == "lightdock_setup":
        structures = out / "structures"
        structures.mkdir()
        names = [f"{pair}_{p}_clean.pdb" for p in ("rec", "lig")]
        for name in names:
            shutil.copy2(require_file(previous / name), structures / name)
        execute(["lgd_setup.py", *(str(Path("structures") / n) for n in names),
                 "-s", settings["swarms"], "-g", settings["glowworms"],
                 "--noxt", "--noh", "--now", "-anm"], out, log)
        setup = json.loads(require_file(out / "setup.json").read_text())
        if not (out / "init").is_dir():
            raise ValueError("LightDock setup produced no init directory.")
        for name in names:
            require_file(structures / ("lightdock_" + name), protein=True)
        if setup.get("use_anm"):
            for p in ("rec", "lig"):
                if setup.get(f"anm_{p}", 0) > 0:
                    require_file(out / f"lightdock_{p}.nm.npy")
    elif stage == "lightdock_simulate":
        # Keep upstream setup outputs immutable in a self-contained run copy.
        shutil.copytree(previous, out, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns("status.json", ".snakemake_timestamp"))
        setup_path = out / "setup.json"
        setup = relocate(json.loads(setup_path.read_text()), previous, out)
        setup_path.write_text(json.dumps(setup, indent=2) + "\n")
        execute(["lgd_run.py", "setup.json", settings["steps"],
                 "-c", threads, "-s", "dfire2"], out, log)
        # Verify final files even if the external command exited normally.
        for swarm in range(settings["swarms"]):
            gso = require_file(out / f"swarm_{swarm}/gso_{settings['steps']}.out")
            if not any(line.startswith("(") for line in gso.read_text().splitlines()):
                raise ValueError(f"No final docking poses in {gso}")
    elif stage == "lightdock_rank":
        ranking = previous / "rank_by_scoring.list"
        ranking.unlink(missing_ok=True)
        execute(["lgd_rank.py", settings["swarms"], settings["steps"]], previous, log)
        require_file(ranking)
        shutil.copy2(ranking, out / f"{pair}_rankings.gso")
        top = out / f"{pair}_top_1.pdb"
        execute([sys.executable, helper / "export_top_pose.py", ranking,
                 previous / f"structures/{pair}_rec_clean.pdb",
                 previous / f"structures/{pair}_lig_clean.pdb",
                 previous / "setup.json", top, settings["steps"]], previous, log)
        require_file(top, protein=True)
        metadata = json.loads(require_file(top.with_suffix(".json")).read_text())
        rows = [line for line in ranking.read_text().splitlines() if line.strip()]
        status["metrics"] = {"top_swarm_id": metadata["swarm"],
                             "top_glowworm_id": metadata["glowworm"],
                             "top_dfire2_score": float(rows[1].split()[-1])}
    elif stage == "visualize_top_pose":
        for suffix in (".pdb", ".json"):
            source = require_file(previous / f"{pair}_top_1{suffix}")
            shutil.copy2(source, out / source.name)
        status["artifacts"]["docked_structure"] = f"{pair}_top_1.pdb"
        status["artifacts"]["selection_metadata"] = f"{pair}_top_1.json"
        simulation = root / "work" / pair / "simulation"
        png = out / f"{pair}_top_pose.png"
        try:
            execute([sys.executable, helper / "plot_top_pose.py",
                     out / f"{pair}_top_1.pdb", png,
                     simulation / f"structures/lightdock_{pair}_rec_clean.pdb",
                     simulation / f"structures/lightdock_{pair}_lig_clean.pdb"], root, log)
            require_file(png)
        except Exception:
            png.unlink(missing_ok=True)
            raise
        status["artifacts"]["figure"] = png.name


def run_stage(stage, pair, root, output, log_path, previous=None, settings=None,
              threads=1, receptor=None, ligand=None, execute=run_command):
    root, out, log_path = [Path(p).resolve() for p in (root, output, log_path)]
    previous = Path(previous).resolve() if previous else None
    if stage not in (*STAGE_DIRS, "visualize_top_pose"):
        raise ValueError(f"Unknown stage: {stage}")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", pair):
        raise ValueError(f"Unsafe pair ID: {pair}")
    expected = (root / "results" / pair if stage == "visualize_top_pose"
                else root / "work" / pair / STAGE_DIRS[stage]).resolve()
    if out != expected or out == root:
        raise ValueError(f"Refusing to clear an unowned output directory: {out}")
    # Clear only this stage's exact owned directory, never an input folder.
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    status = dict(pair=pair, stage=stage, status="success", failed_stage="",
                  error="", log_path="", completed_stages=[], metrics={}, artifacts={})
    # Missing/corrupt upstream reports are infrastructure errors and remain fatal.
    if previous:
        upstream = json.loads((previous / "status.json").read_text())
        status["completed_stages"] = list(upstream["completed_stages"])
        status["metrics"] = dict(upstream["metrics"])
        if upstream["status"] != "success":
            status.update(status="skipped", failed_stage=upstream["failed_stage"],
                          error=upstream["error"], log_path=upstream["log_path"])
    with log_path.open("w", encoding="utf-8") as log:
        if status["status"] == "skipped":
            print(f"SKIPPED {pair}: upstream {status['failed_stage']} failed. "
                  f"See {status['log_path']}", file=log)
        else:
            try:
                perform(stage, pair, root, out, previous, settings or {}, threads,
                        receptor, ligand, log, status, execute)
                status["completed_stages"].append(stage)
                print(f"SUCCESS {pair}: {stage}", file=log)
            except Exception as error:
                status.update(status="failed", failed_stage=stage,
                              error=f"{type(error).__name__}: {error}",
                              log_path=str(log_path.relative_to(root)))
                traceback.print_exc(file=log)
                print(f"FAILED {pair}: {stage}; other pairs can continue.", file=log)
    status["updated_at"] = datetime.now(timezone.utc).isoformat()
    # Reporting errors and Ctrl+C are not swallowed.
    temp = out / "status.json.tmp"
    temp.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    temp.replace(out / "status.json")
    print(f"{pair}: {stage}: {status['status']}" +
          (f" (see {status['log_path']})" if status["failed_stage"] else ""), flush=True)
    return status


if __name__ == "__main__":
    if "snakemake" in globals():
        sm = snakemake
        run_stage(sm.params.stage, sm.wildcards.pair, sm.params.root,
                  sm.output.stage_dir, sm.log[0],
                  previous=getattr(sm.input, "previous", None),
                  settings=json.loads(getattr(sm.params, "settings", "{}")),
                  threads=sm.threads,
                  receptor=getattr(sm.params, "receptor", None),
                  ligand=getattr(sm.params, "ligand", None))
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("stage", choices=STAGE_DIRS)
        for name in ("pair", "root", "previous", "output", "log", "settings"):
            parser.add_argument("--" + name, required=True)
        parser.add_argument("--threads", type=int, default=1)
        args = parser.parse_args()
        run_stage(args.stage, args.pair, args.root, args.output, args.log,
                  previous=args.previous, settings=json.loads(args.settings),
                  threads=args.threads)
