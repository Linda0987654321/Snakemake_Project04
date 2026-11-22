import os
import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile


def fix_and_clean_pdb(input_pdb: str, output_pdb: str, ph: float = 7.0) -> None:
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


if __name__ == "__main__":
    # Works both via Snakemake `script:` directive and standalone CLI
    if "snakemake" in globals():
        sys.stderr = open(snakemake.log[0], "w")
        sys.stdout = sys.stderr
        fix_and_clean_pdb(snakemake.input.rec, snakemake.output.rec)
        fix_and_clean_pdb(snakemake.input.lig, snakemake.output.lig)
    else:
        fix_and_clean_pdb(sys.argv[1], sys.argv[2])