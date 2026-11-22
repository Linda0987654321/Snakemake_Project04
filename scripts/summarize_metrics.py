import csv
import os
import sys

def parse_ranking(gso_path):
    with open(gso_path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    # Line 0 is header, Line 1 is rank #1
    if len(lines) < 2:
        return "NA", "NA", "NA"
    parts = lines[1].split()
    swarm_id = parts[0]
    glowworm_id = parts[1]
    scoring = parts[-1]
    return swarm_id, glowworm_id, scoring

def parse_benchmark(bench_path):
    with open(bench_path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    if len(lines) < 2:
        return "NA", "NA"
    headers = lines[0].split("\t")
    values = lines[1].split("\t")
    row = dict(zip(headers, values))
    wall_time_s = round(float(row.get("s", 0.0)), 2)
    max_rss_mb = round(float(row.get("max_rss", 0.0)), 2)
    return wall_time_s, max_rss_mb

def main():
    output_csv = snakemake.output.summary
    pairs = snakemake.params.pairs
    rankings_files = snakemake.input.rankings
    benchmark_files = snakemake.input.benchmarks

    os.makedirs(os.path.dirname(output_csv), exist_ok=True)

    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "pair",
            "top_swarm_id",
            "top_glowworm_id",
            "top_dfire2_score",
            "simulate_wall_time_s",
            "simulate_max_rss_mb"
        ])

        for pair, rank_file, bench_file in zip(pairs, rankings_files, benchmark_files):
            swarm_id, glowworm_id, score = parse_ranking(rank_file)
            wall_time, max_rss = parse_benchmark(bench_file)
            writer.writerow([pair, swarm_id, glowworm_id, score, wall_time, max_rss])

if __name__ == "__main__":
    import contextlib
    import traceback
    os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
    with open(snakemake.log[0], "w") as log_file:
        with contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file):
            try:
                main()
                print(f"Wrote {snakemake.output.summary}: {len(snakemake.params.pairs)} pairs")
            except Exception:
                traceback.print_exc()
                raise