"""Render partners using their known processed atom counts, including all chains."""
import argparse
from pathlib import Path


def split_partners(top_pdb, receptor_pdb, ligand_pdb):
    def atoms(path):
        return [line for line in Path(path).read_text().splitlines(keepends=True)
                if line.startswith(("ATOM  ", "HETATM"))]

    receptor_atoms = atoms(receptor_pdb)
    ligand_atoms = atoms(ligand_pdb)
    model_atoms = atoms(top_pdb)
    if not receptor_atoms or not ligand_atoms:
        raise ValueError("A processed docking partner contains no atoms.")
    if len(model_atoms) != len(receptor_atoms) + len(ligand_atoms):
        raise ValueError("Exported model atom count does not match its partners.")

    # LightDock writes every receptor atom first, then every ligand atom.
    # Preserve internal chain boundaries even when chain IDs are repeated.
    def pdb_block(lines):
        block = []
        previous_chain = None
        for line in lines:
            chain = line[21:22]
            if previous_chain is not None and chain != previous_chain:
                block.append("TER\n")
            block.append(line)
            previous_chain = chain
        block.extend(["TER\n", "END\n"])
        return "".join(block)

    count = len(receptor_atoms)
    return pdb_block(model_atoms[:count]), pdb_block(model_atoms[count:])


def render_pose(top_pdb, output_png, receptor_pdb, ligand_pdb):
    import pymol
    from pymol import cmd

    receptor, ligand = split_partners(top_pdb, receptor_pdb, ligand_pdb)
    pymol.finish_launching(["pymol", "-cq"])
    cmd.reinitialize()
    cmd.read_pdbstr(receptor, "receptor")
    cmd.read_pdbstr(ligand, "ligand")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.hide("everything")
    cmd.show("cartoon", "receptor")
    cmd.show("surface", "receptor")
    cmd.set("transparency", 0.72, "receptor")
    cmd.color("marine", "receptor")
    cmd.show("cartoon", "ligand")
    cmd.color("tv_orange", "ligand")
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("ray_shadows", 0)
    cmd.orient()
    cmd.zoom("all", buffer=5)
    cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)
    cmd.quit()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("top_pdb", "output_png", "receptor_pdb", "ligand_pdb"):
        parser.add_argument(name)
    args = parser.parse_args()
    render_pose(args.top_pdb, args.output_png, args.receptor_pdb, args.ligand_pdb)
