"""Export the exact winning glowworm with LightDock's ANM-aware exporter."""
import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path


def winning_pose(ranking_path, steps):
    with ranking_path.open() as handle:
        next(handle)
        row = next((line for line in handle if line.strip()), None)
    if row is None:
        raise ValueError("The ranking contains no poses.")
    prefix, _, ranked_pose = row.partition("(")
    swarm, glowworm = map(int, prefix.split())
    source = ranking_path.parent / f"swarm_{swarm}" / f"gso_{steps}.out"
    poses = [line for line in source.read_text().splitlines()
             if line.startswith("(")]
    if glowworm < 0 or glowworm >= len(poses):
        raise ValueError(f"Glowworm {glowworm} is absent from {source}")
    pose = poses[glowworm]
    # Confirm that the ranking points to this pose, allowing ranking rounding.
    ranked_coordinates = [float(x) for x in ranked_pose.split(")")[0].split(",")]
    coordinates = [float(x) for x in pose.split(")")[0][1:].split(",")]
    if len(coordinates) != len(ranked_coordinates):
        raise ValueError("Ranking and swarm coordinate lengths differ.")
    if any(abs(a - b) > 0.00051 for a, b in zip(coordinates, ranked_coordinates)):
        raise ValueError("Ranking and selected swarm pose disagree.")
    if abs(float(pose.split(")")[1].split()[-1]) - float(row.split()[-1])) > 0.00051:
        raise ValueError("Ranking and selected pose scores disagree.")
    return swarm, glowworm, pose


def export_top_pose(ranking, receptor, ligand, setup_file, output, steps):
    ranking, receptor, ligand, setup_file, output = [
        Path(path).resolve() for path in (ranking, receptor, ligand, setup_file, output)
    ]
    setup = json.loads(setup_file.read_text())
    swarm, glowworm, pose = winning_pose(ranking, steps)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".export-", dir=output.parent) as scratch:
        stage = Path(scratch)
        # The exporter discovers normal modes next to its receptor input.
        # Stage structures, masks, and modes together to satisfy that lookup.
        for original in (receptor, ligand):
            processed = original.with_name("lightdock_" + original.name)
            shutil.copy2(processed, stage / processed.name)
            shutil.copy2(original, stage / original.name)
            mask = processed.with_name(processed.stem + "_mask.npy")
            if mask.exists():
                shutil.copy2(mask, stage / mask.name)

        if setup.get("use_anm"):
            for partner in ("rec", "lig"):
                if setup.get(f"anm_{partner}", 0) > 0:
                    modes = setup_file.parent / f"lightdock_{partner}.nm.npy"
                    if not modes.is_file():
                        raise FileNotFoundError(f"Missing ANM modes: {modes}")
                    shutil.copy2(modes, stage / modes.name)

        selected = stage / "selected.out"
        selected.write_text(pose + "\n")
        subprocess.run([
            "lgd_generate_conformations.py", "--setup", str(setup_file),
            str(stage / receptor.name), str(stage / ligand.name),
            str(selected), "1",
        ], check=True)
        generated = stage / "lightdock_0.pdb"
        if not generated.is_file() or not any(
            line.startswith("ATOM  ") for line in generated.read_text().splitlines()
        ):
            raise ValueError("The pose exporter did not produce a protein PDB.")
        generated.replace(output)
    provenance = {
        "swarm": swarm, "glowworm": glowworm,
        "ranking": str(ranking), "steps": steps,
        "anm_applied": bool(setup.get("use_anm")),
    }
    output.with_suffix(".json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("ranking", "receptor", "ligand", "setup", "output"):
        parser.add_argument(name)
    parser.add_argument("steps", type=int)
    args = parser.parse_args()
    export_top_pose(args.ranking, args.receptor, args.ligand,
                    args.setup, args.output, args.steps)
