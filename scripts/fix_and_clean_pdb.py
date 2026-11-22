import os
import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile


def fix_and_clean_pdb(input_pdb: str, output_pdb: str, ph: float = 7.0) -> None:
    print(f"Preparing {input_pdb}", flush=True)
    fixer = PDBFixer(filename=input_pdb)

    # 1. Remove water molecules and non-standard heteroatoms
    fixer.removeHeterogens(keepWater=False)

    # 2. Identify and rebuild missing residues & missing heavy atoms
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # 3. Add hydrogens at physiological pH
    fixer.addMissingHydrogens(ph)

    os.makedirs(os.path.dirname(output_pdb), exist_ok=True)
    with open(output_pdb, "w") as out_f:
        PDBFile.writeFile(fixer.topology, fixer.positions, out_f, keepIds=True)
    print(f"Wrote {output_pdb}", flush=True)


if __name__ == "__main__":
    # Works both via Snakemake `script:` directive and standalone CLI
    if "snakemake" in globals():
        import contextlib
        import traceback
        os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
        with open(snakemake.log[0], "w") as log_file:
            with contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file):
                try:
                    fix_and_clean_pdb(snakemake.input.rec, snakemake.output.rec)
                    fix_and_clean_pdb(snakemake.input.lig, snakemake.output.lig)
                    print("Preparation complete.", flush=True)
                except Exception:
                    traceback.print_exc()
                    raise
    else:
        fix_and_clean_pdb(sys.argv[1], sys.argv[2])