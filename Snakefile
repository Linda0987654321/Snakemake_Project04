from pathlib import Path
import json
import re
from snakemake.utils import min_version
from snakemake.exceptions import WorkflowError

min_version("9.27.1")
configfile: "config/config.yaml"

# Step 1: Discover IDs using either partner. Completeness is checked per pair in
# Step 2, so an incomplete pair is reported without preventing other pairs running.
ROOT = str(Path.cwd().resolve())
INPUT_DIR = str(Path(config["input_dir"]))
input_path = Path(INPUT_DIR).resolve()
for generated in (Path("work").resolve(), Path("results").resolve()):
    if input_path == generated or generated in input_path.parents:
        raise WorkflowError("input_dir must be outside generated work/ and results/.")
REC_IDS = set(glob_wildcards(INPUT_DIR + "/{pair}_u1.pdb").pair)
LIG_IDS = set(glob_wildcards(INPUT_DIR + "/{pair}_u2.pdb").pair)
PAIRS = sorted(REC_IDS | LIG_IDS)
if not PAIRS:
    raise WorkflowError("No protein pairs found in " + INPUT_DIR +
                        ". Add <pair>_u1.pdb and <pair>_u2.pdb, or select data/example.")
for pair in PAIRS:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", pair):
        raise WorkflowError("Invalid pair ID: " + pair)

SETTINGS = config["docking_params"]
for name in ("swarms", "glowworms", "steps"):
    if type(SETTINGS[name]) is not int or SETTINGS[name] < 1:
        raise WorkflowError(name + " must be a positive integer.")
SETTINGS_JSON = json.dumps({name: SETTINGS[name] for name in
                            ("swarms", "glowworms", "steps")})

def existing_partners(wildcards):
    # Track every existing source, but let Step 2 report a missing partner.
    return [str(p) for n in (1, 2)
            if (p := Path(INPUT_DIR) / f"{wildcards.pair}_u{n}.pdb").is_file()]

wildcard_constraints:
    pair = r"[A-Za-z0-9][A-Za-z0-9_-]*"

# Every pair's final folder always has status.json. Successful pairs also have
# real structures/figures. A failed pair never receives dummy scientific outputs.
# Directory outputs allow optional artifacts without hiding a failure in the CSV.
rule all:
    input:
        expand("results/{pair}", pair=PAIRS),
        "results/summary_metrics.csv"

# Step 2: Validate both files, repair structures with PDBFixer and record status.
# The stage runner catches pair-level errors; its log keeps the full traceback.
rule clean_structures:
    input:
        sources = existing_partners,
        helper = "scripts/fix_and_clean_pdb.py"
    output:
        stage_dir = directory("work/{pair}/cleaned")
    params:
        stage = "clean_structures",
        root = ROOT,
        receptor = lambda wc: str(Path(INPUT_DIR) / f"{wc.pair}_u1.pdb"),
        ligand = lambda wc: str(Path(INPUT_DIR) / f"{wc.pair}_u2.pdb")
    log:
        "logs/{pair}/clean_structures.log"
    conda:
        "envs/workflow_tools.yaml"
    script:
        "scripts/run_stage.py"

# Step 3: lgd_setup.py (from LightDock) creates starting poses and ANM modes.
# The runner skips this pair if preparation failed, recording the original error.
rule lightdock_setup:
    input:
        previous = "work/{pair}/cleaned",
        runner = "scripts/run_stage.py"
    output:
        stage_dir = directory("work/{pair}/setup")
    params:
        root = ROOT,
        settings = SETTINGS_JSON
    log:
        "logs/{pair}/lightdock_setup.log"
    conda:
        "envs/workflow_tools.yaml"
    shell:
        "python {input.runner:q} lightdock_setup --pair {wildcards.pair:q} "
        "--root {params.root:q} --previous {input.previous:q} "
        "--output {output.stage_dir:q} --log {log:q} --settings {params.settings:q}"

# Step 4: lgd_run.py performs the actual dfire2 docking optimization.
# Its final swarm files are verified inside the stage directory before success.
rule lightdock_simulate:
    input:
        previous = "work/{pair}/setup",
        runner = "scripts/run_stage.py"
    output:
        stage_dir = directory("work/{pair}/simulation")
    params:
        root = ROOT,
        settings = SETTINGS_JSON
    threads: 4
    log:
        "logs/{pair}/lightdock_simulate.log"
    benchmark:
        "benchmarks/{pair}/simulate.txt"
    conda:
        "envs/workflow_tools.yaml"
    shell:
        "python {input.runner:q} lightdock_simulate --pair {wildcards.pair:q} "
        "--root {params.root:q} --previous {input.previous:q} "
        "--output {output.stage_dir:q} --log {log:q} --settings {params.settings:q} "
        "--threads {threads}"

# Step 5: lgd_rank.py ranks candidates. Our exporter reconstructs the exact
# winning glowworm, retaining its ANM deformation and selection metadata.
rule lightdock_rank:
    input:
        previous = "work/{pair}/simulation",
        runner = "scripts/run_stage.py",
        exporter = "scripts/export_top_pose.py"
    output:
        stage_dir = directory("work/{pair}/ranked")
    params:
        root = ROOT,
        settings = SETTINGS_JSON
    log:
        "logs/{pair}/lightdock_rank.log"
    benchmark:
        "benchmarks/{pair}/rank.txt"
    conda:
        "envs/workflow_tools.yaml"
    shell:
        "python {input.runner:q} lightdock_rank --pair {wildcards.pair:q} "
        "--root {params.root:q} --previous {input.previous:q} "
        "--output {output.stage_dir:q} --log {log:q} --settings {params.settings:q}"

# Step 6: Copy the docked model to results and render it using PyMOL.
# If rendering fails, retain the valid model and report the missing figure.
rule visualize_top_pose:
    input:
        previous = "work/{pair}/ranked",
        simulation = "work/{pair}/simulation",
        renderer = "scripts/plot_top_pose.py"
    output:
        stage_dir = directory("results/{pair}")
    params:
        stage = "visualize_top_pose",
        root = ROOT
    log:
        "logs/{pair}/visualize.log"
    conda:
        "envs/workflow_tools.yaml"
    script:
        "scripts/run_stage.py"

# Step 7: Wait for all pair reports and write one CSV, including failed pairs.
# Metrics are blank when their stage did not complete successfully.
rule summarize_results:
    input:
        reports = expand("results/{pair}", pair=PAIRS),
        benchmarks = expand("benchmarks/{pair}/simulate.txt", pair=PAIRS)
    output:
        summary = "results/summary_metrics.csv"
    params:
        pairs = PAIRS
    log:
        "logs/summary_metrics.log"
    script:
        "scripts/summarize_metrics.py"
