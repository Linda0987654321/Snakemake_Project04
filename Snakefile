from pathlib import Path
import re
from snakemake.utils import min_version
from snakemake.exceptions import WorkflowError

min_version("9.27.1")
configfile: "config/config.yaml"

# Step 1: Discover local receptor/ligand pairs before building the job graph.
# u1 is the receptor, u2 the ligand. The prefix is the reusable {pair} wildcard.
INPUT_DIR = str(config["input_dir"]).rstrip("/")
REC_IDS = set(glob_wildcards(INPUT_DIR + "/{pair}_u1.pdb").pair)
LIG_IDS = set(glob_wildcards(INPUT_DIR + "/{pair}_u2.pdb").pair)
PAIRS = sorted(REC_IDS | LIG_IDS)
if not PAIRS:
    raise WorkflowError("No protein pairs found in " + INPUT_DIR +
                        ". Add <pair>_u1.pdb and <pair>_u2.pdb.")
for pair in PAIRS:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", pair):
        raise WorkflowError("Invalid pair ID: " + pair +
                            ". Use letters, digits, underscores or hyphens.")
    for partner in (1, 2):
        path = Path(INPUT_DIR) / f"{pair}_u{partner}.pdb"
        if not path.is_file() or path.stat().st_size == 0:
            raise WorkflowError("Missing or empty partner file: " + str(path))

# Read settings once. Validate early to avoid partial runs with invalid settings.
SWARMS = config["docking_params"]["swarms"]
GLOWWORMS = config["docking_params"]["glowworms"]
STEPS = config["docking_params"]["steps"]
for name, value in [("swarms", SWARMS), ("glowworms", GLOWWORMS), ("steps", STEPS)]:
    if type(value) is not int or value < 1:
        raise WorkflowError(name + " must be a positive integer.")

wildcard_constraints:
    pair = r"[A-Za-z0-9][A-Za-z0-9_-]*"

# Default target: request a structure, metadata and figure for every detected pair.
# expand() fills in the pair IDs. The summary is produced by a separate rule.
rule all:
    input:
        expand("results/{pair}/{pair}_top_1.pdb", pair=PAIRS),
        expand("results/{pair}/{pair}_top_1.json", pair=PAIRS),
        expand("results/{pair}/{pair}_top_pose.png", pair=PAIRS),
        "results/summary_metrics.csv"

# Step 2: Read the local input PDBs and write repaired structures under work/.
# script: makes named input/output/log paths available to the Python script.
rule clean_structures:
    input:
        rec = INPUT_DIR + "/{pair}_u1.pdb",
        lig = INPUT_DIR + "/{pair}_u2.pdb"
    output:
        rec = "work/{pair}/cleaned/{pair}_rec_clean.pdb",
        lig = "work/{pair}/cleaned/{pair}_lig_clean.pdb"
    log:
        "logs/{pair}/clean_structures.log"
    conda:
        "envs/workflow_tools.yaml"
    script:
        "scripts/fix_and_clean_pdb.py"

# Step 3: Create starting poses, processed PDBs, masks and ANM normal modes.
# Water/hydrogen/OXT are excluded for scoring after PDBFixer preparation.
# The complete tool environment is reused by preparation, docking and plotting.
rule lightdock_setup:
    input:
        rec = "work/{pair}/cleaned/{pair}_rec_clean.pdb",
        lig = "work/{pair}/cleaned/{pair}_lig_clean.pdb"
    output:
        setup = "work/{pair}/setup.json",
        init_dir = directory("work/{pair}/init"),
        rec = "work/{pair}/cleaned/lightdock_{pair}_rec_clean.pdb",
        lig = "work/{pair}/cleaned/lightdock_{pair}_lig_clean.pdb"
    params:
        swarms = SWARMS,
        glowworms = GLOWWORMS
    log:
        "logs/{pair}/lightdock_setup.log"
    conda:
        "envs/workflow_tools.yaml"
    shell:
        """
        mkdir -p work/{wildcards.pair} logs/{wildcards.pair}
        REC=$(realpath {input.rec:q})
        LIG=$(realpath {input.lig:q})
        LOG=$(realpath {log:q})
        RUN=$(realpath work/{wildcards.pair})
        # Clear only this pair's derived setup/simulation files on a setup rerun.
        rm -f "$RUN"/cleaned/lightdock_* "$RUN"/cleaned/*.points
        rm -rf "$RUN"/swarm_* "$RUN"/init
        rm -f "$RUN"/setup.json "$RUN"/*.npy
        cd "$RUN"
        lgd_setup.py "$REC" "$LIG" -s {params.swarms} -g {params.glowworms} \
            --noxt --noh --now -anm > "$LOG" 2>&1
        """

# Step 4: Run the actual optimization and declare every final swarm file.
# Double braces preserve {pair} while expand() inserts swarm IDs and steps.
rule lightdock_simulate:
    input:
        setup = "work/{pair}/setup.json",
        init_dir = "work/{pair}/init"
    output:
        gso_files = expand("work/{{pair}}/swarm_{swarm}/gso_{steps}.out",
                           swarm=range(SWARMS), steps=[STEPS])
    params:
        steps = STEPS
    threads: 4
    log:
        "logs/{pair}/lightdock_simulate.log"
    benchmark:
        "benchmarks/{pair}/simulate.txt"
    conda:
        "envs/workflow_tools.yaml"
    shell:
        """
        LOG=$(realpath {log:q})
        cd work/{wildcards.pair}
        lgd_run.py setup.json {params.steps} -c {threads} -s dfire2 > "$LOG" 2>&1
        """

# Step 5: Rank all candidate poses and reconstruct the exact winning model.
# The helper retains the winning glowworm and stages ANM modes/masks correctly.
rule lightdock_rank:
    input:
        gso_files = expand("work/{{pair}}/swarm_{swarm}/gso_{steps}.out",
                           swarm=range(SWARMS), steps=[STEPS]),
        rec = "work/{pair}/cleaned/{pair}_rec_clean.pdb",
        lig = "work/{pair}/cleaned/{pair}_lig_clean.pdb",
        setup = "work/{pair}/setup.json",
        exporter = "scripts/export_top_pose.py"
    output:
        rankings = "work/{pair}/{pair}_rankings.gso",
        top_pdb = "results/{pair}/{pair}_top_1.pdb",
        provenance = "results/{pair}/{pair}_top_1.json"
    params:
        swarms = SWARMS,
        steps = STEPS
    log:
        "logs/{pair}/lightdock_rank.log"
    benchmark:
        "benchmarks/{pair}/rank.txt"
    conda:
        "envs/workflow_tools.yaml"
    shell:
        """
        mkdir -p results/{wildcards.pair}
        REC=$(realpath {input.rec:q})
        LIG=$(realpath {input.lig:q})
        LOG=$(realpath {log:q})
        EXPORTER=$(realpath {input.exporter:q})
        RANKINGS=$(realpath -m {output.rankings:q})
        TOP=$(realpath -m {output.top_pdb:q})
        cd work/{wildcards.pair}
        lgd_rank.py {params.swarms} {params.steps} > "$LOG" 2>&1
        cp rank_by_scoring.list "$RANKINGS"
        python "$EXPORTER" rank_by_scoring.list "$REC" "$LIG" \
            setup.json "$TOP" {params.steps} >> "$LOG" 2>&1
        """

# Step 6: Render the exported structure. Known processed atom counts identify
# both complete partners, including multiple chains, and preserve cartoon bonds.
rule visualize_top_pose:
    input:
        top_pdb = "results/{pair}/{pair}_top_1.pdb",
        rec = "work/{pair}/cleaned/lightdock_{pair}_rec_clean.pdb",
        lig = "work/{pair}/cleaned/lightdock_{pair}_lig_clean.pdb",
        renderer = "scripts/plot_top_pose.py"
    output:
        image = "results/{pair}/{pair}_top_pose.png"
    log:
        "logs/{pair}/visualize.log"
    conda:
        "envs/workflow_tools.yaml"
    shell:
        """
        python {input.renderer:q} {input.top_pdb:q} {output.image:q} \
            {input.rec:q} {input.lig:q} > {log:q} 2>&1
        """

# Step 7: Combine one ranking and simulation benchmark from each detected pair.
# One shared output table contains one row per pair, not one table per pair.
rule summarize_results:
    input:
        rankings = expand("work/{pair}/{pair}_rankings.gso", pair=PAIRS),
        benchmarks = expand("benchmarks/{pair}/simulate.txt", pair=PAIRS)
    output:
        summary = "results/summary_metrics.csv"
    params:
        pairs = PAIRS
    log:
        "logs/summary_metrics.log"
    script:
        "scripts/summarize_metrics.py"
