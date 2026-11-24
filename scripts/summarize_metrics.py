"""Combine every pair report, including failed and skipped pairs, in one CSV."""
import contextlib
import csv
import json
from pathlib import Path

FIELDS = ["pair", "status", "failed_stage", "error", "log_path",
          "top_swarm_id", "top_glowworm_id", "top_dfire2_score",
          "simulate_wall_time_s", "simulate_max_rss_mb",
          "docked_structure", "figure"]


def summarize(pairs, reports, benchmarks, output):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    failures = 0
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for pair, folder, benchmark in zip(pairs, reports, benchmarks, strict=True):
            folder = Path(folder)
            report = json.loads((folder / "status.json").read_text())
            if report["pair"] != pair:
                raise ValueError(f"Mismatched pair report: {folder}")
            failed = report["status"] != "success"
            failures += int(failed)
            row = {"pair": pair, "status": "failed" if failed else "success",
                   "failed_stage": report["failed_stage"], "error": report["error"],
                   "log_path": report["log_path"], **report["metrics"]}
            if "lightdock_simulate" in report["completed_stages"]:
                with Path(benchmark).open() as source:
                    stats = next(csv.DictReader(source, delimiter="\t"))
                row.update(simulate_wall_time_s=round(float(stats["s"]), 2),
                           simulate_max_rss_mb=round(float(stats["max_rss"]), 2))
            for kind in ("docked_structure", "figure"):
                if kind in report["artifacts"]:
                    artifact = folder / report["artifacts"][kind]
                    if not artifact.is_file():
                        raise ValueError(f"Reported artifact is absent: {artifact}")
                    row[kind] = str(artifact)
            writer.writerow(row)
    message = f"Batch report: {len(pairs) - failures} successful, {failures} failed pairs."
    print(message)
    if failures:
        print("WARNING: Snakemake finished the batch, but some pairs failed. "
              "Check status, failed_stage and log_path in the CSV.")
    return message


if __name__ == "__main__":
    Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
    with open(snakemake.log[0], "w", encoding="utf-8") as log:
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            message = summarize(snakemake.params.pairs, snakemake.input.reports,
                                snakemake.input.benchmarks, snakemake.output.summary)
    print(message + " See " + str(snakemake.output.summary), flush=True)
