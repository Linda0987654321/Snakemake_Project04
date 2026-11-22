import sys
import tempfile
import os
import pymol
from pymol import cmd

def render_pose(top_pdb, output_png):
    # Split top_1.pdb at the first TER record into separate receptor and ligand PDB blocks
    with open(top_pdb, "r") as f:
        lines = f.readlines()

    rec_lines, lig_lines = [], []
    hit_ter = False
    for line in lines:
        if line.startswith("TER") and not hit_ter:
            hit_ter = True
            rec_lines.append(line)
            continue
        if line.startswith(("ATOM", "HETATM")):
            if not hit_ter:
                rec_lines.append(line)
            else:
                lig_lines.append(line)

    # Fallback if no TER record exists: split when residue sequence number drops
    if not lig_lines:
        rec_lines, lig_lines = [], []
        prev_res = -99999
        in_lig = False
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                res_num = int(line[22:26])
                if not in_lig and res_num < prev_res:
                    in_lig = True
                prev_res = res_num
                (lig_lines if in_lig else rec_lines).append(line)

    with tempfile.TemporaryDirectory() as tmpdir:
        rec_tmp = os.path.join(tmpdir, "rec.pdb")
        lig_tmp = os.path.join(tmpdir, "lig.pdb")
        with open(rec_tmp, "w") as f:
            f.writelines(rec_lines)
        with open(lig_tmp, "w") as f:
            f.writelines(lig_lines)

        pymol.finish_launching(['pymol', '-cq'])
        cmd.reinitialize()

        cmd.load(rec_tmp, "receptor")
        cmd.load(lig_tmp, "ligand")
        cmd.unbond("all", "all")

        # Visual styling with opaque white background
        cmd.bg_color("white")
        cmd.set("ray_opaque_background", 1)
        cmd.hide("everything")

        # Receptor: marine cartoon + translucent surface
        cmd.show("cartoon", "receptor")
        cmd.show("surface", "receptor")
        cmd.set("transparency", 0.72, "receptor")
        cmd.color("marine", "receptor")

        # Ligand: high-contrast orange cartoon
        cmd.show("cartoon", "ligand")
        cmd.color("tv_orange", "ligand")

        # Smooth ribbon rendering & framing
        cmd.set("cartoon_smooth_loops", 1)
        cmd.set("ray_shadows", 0)
        cmd.orient()
        cmd.zoom("all", buffer=5)
        cmd.png(output_png, width=1600, height=1200, dpi=300, ray=1)
        cmd.quit()

if __name__ == "__main__":
    render_pose(sys.argv[1], sys.argv[2])