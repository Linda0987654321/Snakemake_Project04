from snakemake.utils import min_version
min_version("9.27.1")

configfile: "config/config.yaml"

# Read the pair IDs and shared docking settings from the config file.
# The same rules are used for each pair; {pair} identifies its own files.
PAIRS = list(config["pairs"].keys())
SWARMS = config["docking_params"]["swarms"]
GLOWWORMS = config["docking_params"]["glowworms"]
STEPS = config["docking_params"]["steps"]

# Final targets. expand() requests one figure per configured pair.
# The summary CSV requires the rankings and benchmarks from every pair.
# Snakemake starts here and finds the upstream jobs needed to make them.
rule all:
    input:
        expand("results/{pair}/visualization/{pair}_top_pose.png", pair=PAIRS),
        "results/summary_metrics.csv"

# Step 1: Download both original protein structures for a pair.
# The lambda functions look up receptor and ligand URLs using the pair
# wildcard. shell: runs curl and checks that both downloaded files have content.
rule download_pdbs:
    output:
        rec = "data/{pair}_rec.pdb",
        lig = "data/{pair}_lig.pdb"
    params:
        rec_url = lambda wildcards: config["pairs"][wildcards.pair]["receptor_url"],
        lig_url = lambda wildcards: config["pairs"][wildcards.pair]["ligand_url"]
    log:
        "logs/{pair}/download.log"
    shell:
        """
        mkdir -p data logs/{wildcards.pair}
        set -euo pipefail

        curl -fSL "{params.rec_url}" -o {output.rec} 2> {log}
        curl -fSL "{params.lig_url}" -o {output.lig} 2>> {log}
        
        test -s {output.rec}
        test -s {output.lig}
        """

# Step 2: Clean and preprocess both proteins with PDBFixer.
# script: passes the named input/output paths to the Python program
# through snakemake.input and snakemake.output. It removes heterogens,
# models detected missing residues, repairs atoms and adds hydrogens.
rule clean_structures:
    input:
        rec = "data/{pair}_rec.pdb",
        lig = "data/{pair}_lig.pdb"
    output:
        rec = "results/{pair}/cleaned/{pair}_rec_clean.pdb",
        lig = "results/{pair}/cleaned/{pair}_lig_clean.pdb"
    log:
        "logs/{pair}/clean_structures.log"
    conda:
        "envs/lightdock_env.yaml"
    script:
        "scripts/fix_and_clean_pdb.py"

# Step 3: Set up LightDock using the cleaned receptor and ligand.
# Setup creates the starting poses, processed PDBs, masks and ANM modes.
# ANM enables backbone flexibility. Water, hydrogen and OXT atoms are
# excluded from the scoring inputs. Absolute paths work after cd below.
rule lightdock_setup:
    input:
        rec = "results/{pair}/cleaned/{pair}_rec_clean.pdb",
        lig = "results/{pair}/cleaned/{pair}_lig_clean.pdb"
    output:
        setup = "results/{pair}/setup.json",
        init_dir = directory("results/{pair}/init")
    params:
        swarms = SWARMS,
        glowworms = GLOWWORMS
    log:
        "logs/{pair}/lightdock_setup.log"
    conda:
        "envs/lightdock_env.yaml"
    shell:
        """
        REC=$(realpath {input.rec})
        LIG=$(realpath {input.lig})
        LOG=$(realpath {log})
        mkdir -p results/{wildcards.pair} logs/{wildcards.pair}
        
	# Remove any stale LightDock setup files from previous runs/retries
        rm -f results/{wildcards.pair}/cleaned/lightdock_* results/{wildcards.pair}/cleaned/*.points
        rm -rf results/{wildcards.pair}/swarm_* results/{wildcards.pair}/init results/{wildcards.pair}/setup.json results/{wildcards.pair}/*.npy
	
	cd results/{wildcards.pair}
        lgd_setup.py "$REC" "$LIG" -s {params.swarms} -g {params.glowworms} --noxt --noh --now -anm > "$LOG" 2>&1
        """

# Step 4: Optimize the docking poses using dfire2 and ANM flexibility.
# expand() fills in the swarm IDs and final step number. {{pair}} keeps
# the pair wildcard for Snakemake to resolve separately for each job.
# threads requests up to four cores; benchmark records runtime and memory.
rule lightdock_simulate:
    input:
        setup = "results/{pair}/setup.json",
        init_dir = "results/{pair}/init"
    output:
        gso_files = expand(
            "results/{{pair}}/swarm_{swarm}/gso_{steps}.out",
            swarm=range(SWARMS),
            steps=[STEPS]
        )
    params:
        steps = STEPS
    threads: 4
    log:
        "logs/{pair}/lightdock_simulate.log"
    benchmark:
        "benchmarks/{pair}/simulate.txt"
    conda:
        "envs/lightdock_env.yaml"
    shell:
        """
        LOG=$(realpath {log})
        cd results/{wildcards.pair}
        lgd_run.py setup.json {params.steps} -c {threads} -s dfire2 > "$LOG" 2>&1
        """

# Step 5: Rank the poses and export the exact highest-ranked model.
# The helper selects both the winning swarm and glowworm. It also stages
# the ANM modes beside the processed structures so LightDock applies
# the flexible deformation during reconstruction. JSON records the IDs.
rule lightdock_rank:
    input:
        gso_files = expand(
            "results/{{pair}}/swarm_{swarm}/gso_{steps}.out",
            swarm=range(SWARMS),
            steps=[STEPS]
        ),
        rec = "results/{pair}/cleaned/{pair}_rec_clean.pdb",
        lig = "results/{pair}/cleaned/{pair}_lig_clean.pdb",
        setup = "results/{pair}/setup.json",
        exporter = "scripts/export_top_pose.py"
    output:
        rankings = "results/{pair}/{pair}_rankings.gso",
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
        "envs/lightdock_env.yaml"
    shell:
        """
        REC=$(realpath {input.rec})
        LIG=$(realpath {input.lig})
        LOG=$(realpath {log})
        cd results/{wildcards.pair}

        lgd_rank.py {params.swarms} {params.steps} > "$LOG" 2>&1
        RANK_FILE="rank_by_scoring.list"
        cp "$RANK_FILE" {wildcards.pair}_rankings.gso

        # The exporter validates the winning IDs against the swarm output.
        python ../../scripts/export_top_pose.py "$RANK_FILE" "$REC" "$LIG" \
            setup.json {wildcards.pair}_top_1.pdb {params.steps} >> "$LOG" 2>&1
        """

# Step 6: Render the exported complex in PyMOL.
# The processed partner atom counts identify the complete receptor and
# ligand, including multichain partners. Blue shows the receptor with
# a translucent surface; orange shows the ligand.
rule visualize_top_pose:
    input:
        top_pdb = "results/{pair}/{pair}_top_1.pdb",
        rec = "results/{pair}/cleaned/{pair}_rec_clean.pdb",
        lig = "results/{pair}/cleaned/{pair}_lig_clean.pdb",
        renderer = "scripts/plot_top_pose.py"
    output:
        image = "results/{pair}/visualization/{pair}_top_pose.png"
    log:
        "logs/{pair}/visualize.log"
    conda:
        "envs/lightdock_env.yaml"
    shell:
        """
        mkdir -p results/{wildcards.pair}/visualization
        python {input.renderer:q} {input.top_pdb:q} {output.image:q} \
            "$(dirname {input.rec:q})/lightdock_$(basename {input.rec:q})" \
            "$(dirname {input.lig:q})/lightdock_$(basename {input.lig:q})" > {log:q} 2>&1
        """

# Step 7: Summarize the top docking score and simulation benchmark.
# expand() supplies one ranking and one benchmark for every pair.
# The script writes one row per pair to results/summary_metrics.csv.
rule summarize_results:
    input:
        rankings = expand("results/{pair}/{pair}_rankings.gso", pair=PAIRS),
        benchmarks = expand("benchmarks/{pair}/simulate.txt", pair=PAIRS)
    output:
        summary = "results/summary_metrics.csv"
    params:
        pairs = PAIRS
    script:
        "scripts/summarize_metrics.py"
