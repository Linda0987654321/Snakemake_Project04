# Reproducible Multi-Pair Protein-Protein Docking Pipeline (Snakemake + LightDock)

An end-to-end, Snakemake-based workflow designed for reproducible execution. Built with **Snakemake** and **LightDock** to download, repair, dock, rank, benchmark, and visualize multiple protein-protein complexes in parallel.

## Repository Structure

```text
protein_docking/
├── Snakefile                       # Multi-pair workflow rules & full swarm DAG dependencies
├── config/
│   └── config.yaml                 # Protein complex URLs & docking parameters
├── envs/
│   └── lightdock_env.yaml          # Conda/pip environment (LightDock, PDBFixer, PyMOL, Biopython)
├── scripts/
│   ├── fix_and_clean_pdb.py        # PDBFixer preprocessing (waters/heterogens, missing residues, H at pH 7.0)
│   ├── plot_top_pose.py            # Headless PyMOL ray-traced 3D visualization
│   └── summarize_metrics.py        # Cross-pair aggregation of dfire2 scores & benchmarks
├── benchmarks/                     # Snakemake wall-clock & memory benchmarks per complex
├── results/
│   ├── summary_metrics.csv         # Aggregated cross-pair scores & runtime metrics
│   ├── EnzymeI_HPr/                # Rankings, top pose PDB, and rendered PNG
│   └── Trypsin_BPTI/               # Rankings, top pose PDB, and rendered PNG
├── dag_multipair.png               # Full multi-pair execution DAG
├── rulegraph.png                   # High-level rule dependency graph
├── .gitignore                      # Excludes heavy caches & intermediate trajectory files
└── README.md
```

## Workflow Architecture (`dag_multipair.png`)

![Multi-Pair Workflow DAG](dag_multipair.png)

1. **`download_pdbs`**: Fetches unbound receptor and ligand structures directly from RCSB PDB (`1EZA`/`1HDN` and `1TPO`/`4PTI`).
2. **`clean_structures`**: Uses **PDBFixer** (`OpenMM`) to strip crystallographic waters and heterogens, rebuild missing residues and missing heavy atoms, and add hydrogens at physiological pH (`7.0`).
3. **`lightdock_setup`**: Initializes simulation swarms (`10` swarms, `20` glowworms/swarm) with Anisotropic Network Model (`-anm`) backbone flexibility (`--noxt --noh --now -anm`).
4. **`lightdock_simulate`**: Runs Glowworm Swarm Optimization (GSO) for `20` steps across `4` CPU cores using the `dfire2` scoring function, explicitly declaring all swarm outputs (`swarm_0..9/gso_20.out`) in the DAG.
5. **`lightdock_rank`**: Ranks all glowworms across all `10` swarms globally (`lgd_rank.py`) and reconstructs the #1 ranked 3D complex (`{pair}_top_1.pdb`).
6. **`visualize_top_pose`**: Renders a high-resolution ($1600 \times 1200$, 300 DPI) ray-traced PyMOL figure (`{pair}_top_pose.png`) showing the receptor (marine cartoon + translucent surface) and bound ligand (orange cartoon).
7. **`summarize_results`**: Aggregates top `dfire2` scores, winning swarm/glowworm IDs, and wall-clock/memory benchmarks into `results/summary_metrics.csv`.

## Configuration (`config/config.yaml`)

```yaml
pairs:
  EnzymeI_HPr:
    receptor_url: "https://files.rcsb.org/download/1EZA.pdb"
    ligand_url: "https://files.rcsb.org/download/1HDN.pdb"
  Trypsin_BPTI:
    receptor_url: "https://files.rcsb.org/download/1TPO.pdb"
    ligand_url: "https://files.rcsb.org/download/4PTI.pdb"

docking_params:
  swarms: 10
  glowworms: 20
  steps: 20
```

## How to Run

Ensure Conda/Mamba and Snakemake are installed, then execute:

```bash
# Dry-run to verify DAG and target files
snakemake -n --use-conda

# Execute full workflow across 4 cores
snakemake --cores 4 --use-conda
```

## Benchmark & Docking Results (`results/summary_metrics.csv`)

| Complex | Receptor (PDB) | Ligand (PDB) | Top Swarm | Top Glowworm | Top `dfire2` Score | Simulation Wall Time (s) | Max RSS Memory (MB) |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **`EnzymeI_HPr`** | *E. coli* Enzyme I N-terminal domain (`1EZA`) | *E. coli* HPr (`1HDN`) | `2` | `0` | `-406.546` | `12.64` | `634.22` |
| **`Trypsin_BPTI`** | Bovine $\beta$-trypsin (`1TPO`) | Bovine pancreatic trypsin inhibitor (`4PTI`) | `6` | `16` | `-364.047` | `8.84` | `614.22` |

## Rendered Top Docking Poses

| Enzyme I (N-terminal domain) – HPr (`1EZA` + `1HDN`) | $\beta$-Trypsin – BPTI (`1TPO` + `4PTI`) |
| :---: | :---: |
| ![EnzymeI_HPr Top Pose](results/EnzymeI_HPr/visualization/EnzymeI_HPr_top_pose.png) | ![Trypsin_BPTI Top Pose](results/Trypsin_BPTI/visualization/Trypsin_BPTI_top_pose.png) |