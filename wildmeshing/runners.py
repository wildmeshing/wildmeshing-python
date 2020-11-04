import argparse
import json
import wildmeshing as wm
import os.path
import numpy as np


def triangulate():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, default="",
                        help="Input segments in .obj format.")

    parser.add_argument("-o", "--output", type=str,
                        default="",  help="Output path.")
    parser.add_argument("--feature-input", type=str, default="",
                        help="Input feature json file.")
    parser.add_argument("--stop-energy", type=float, default=10,
                        help="Specify max AMIPS energy for stopping mesh optimization.")
    parser.add_argument("-e", "--epsr", type=float, default=1e-3,
                        help="relative envelope epsilon_r. Absolute epsilonn = epsilon_r * diagonal_of_bbox")
    parser.add_argument("--feature-envelope",  type=float, default=0.001,
                        help="Relative feature envelope mu_r. Absolute mu = mu_r * diagonal_of_bbox")
    parser.add_argument("-l", "--lr", type=float, default=0.05,
                        help="ideal_edge_length = diag_of_bbox * L. (double, optional, default: 0.05)")
    parser.add_argument("--mute-log", type=bool,
                        default=False, help="Mute prints.")
    parser.add_argument("--cut-outside", type=bool, default=True,
                        help="Remove \"outside part\".")
    parser.add_argument("--skip-eps", type=bool,
                        default=True, help="Skip saving eps.")

    parser.add_argument("--cut-holes", type=str, default="",
                        help="Input a .xyz file for specifying points inside holes you want to remove.")

    args = parser.parse_args()

    extension = os.path.splitext(args.input)[1]

    if extension == ".svg":
        with np.load(args.cut_holes) as holes:
            V_out, F_out, nodes, F_nodes = wm.triangulate_svg(args.input,
                                                              stop_quality=args.stop_energy,
                                                              epsilon=args.epsr,
                                                              feature_epsilon=args.feature_envelope,
                                                              target_edge_len=args.lr,
                                                              hole_pts=holes,
                                                              cut_outside=args.cut_outside,
                                                              skip_eps=args.skip_eps,
                                                              mute_log=args.mute_log)
            np.savetxt(args.output + "_V.txt", V_out)
            np.savetxt(args.output + "_F.txt", F_out)
            np.savetxt(args.output + "_Vn.txt", nodes)
            np.savetxt(args.output + "_Fn.txt", F_nodes)
    else:
        wm.triangulate(args.input,
                       feature_input=args.feature_input,
                       stop_quality=args.stop_energy,
                       epsilon=args.epsr,
                       feature_epsilon=args.feature_envelope,
                       target_edge_len=args.lr,
                       hole_file=args.cut_holes,
                       cut_outside=args.cut_outside,
                       output=args.output,
                       skip_eps=args.skip_eps,
                       mute_log=args.mute_log)


def tetrahedralize():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input surface mesh INPUT in .off/.obj/.stl/.ply format. (string, required)")
    parser.add_argument("-o", "--output", type=str,
                        help="Output tetmesh OUTPUT in .msh format. (string, optional, default: input_file+postfix+'.msh')")

    parser.add_argument("-l", "--lr", type=float, default=0.05,
                        help="ideal_edge_length = diag_of_bbox * L. (double, optional, default: 0.05)")
    parser.add_argument("-e", "--epsr", type=float, default=1e-3,
                        help="epsilon = diag_of_bbox * EPS. (double, optional, default: 1e-3)")

    parser.add_argument("--stop-energy", type=float, default=10,
                        help="Stop optimization when max energy is lower than this.")

    parser.add_argument("--level", type=int, default=2,
                        help="Log level (0 = most verbose, 6 = off).")

    parser.add_argument("--skip-simplify", type=bool,
                        default=False, help="skip preprocessing.")
    parser.add_argument("--no-binary", type=bool,
                        default=False, help="export meshes as ascii")

    parser.add_argument("--smooth-open-boundary", type=bool,
                        default=False, help="Smooth the open boundary.")
    parser.add_argument("--manifold-surface", type=bool,
                        default=False, help="Force the output to be manifold.")
    parser.add_argument("--coarsen", type=bool, default=True,
                        help="Coarsen the output as much as possible.")
    parser.add_argument("--csg", type=str, default="",
                        help="json file containg a csg tree")

    parser.add_argument("--disable-filtering", type=bool, default=False,
                        help="Disable filtering out outside elements.")
    parser.add_argument("--use-floodfill", type=bool, default=False,
                        help="Use flood-fill to extract interior volume.")
    parser.add_argument("--use-input-for-wn", type=bool,
                        default=False, help="Use input surface for winding number.")

    parser.add_argument("--bg-mesh", type=str, default="",
                        help="Background mesh for sizing field (.msh file).")
    parser.add_argument("--epsr-tags", type=str, default="",
                        help="List of envelope size for each input faces.")

    args = parser.parse_args()

    tetra = wm.Tetrahedralizer(stop_quality=args.stop_energy,
                               epsilon=args.epsr,
                               edge_length_r=args.lr,
                               skip_simplify=args.skip_simplify,
                               coarsen=args.coarsen)
    tetra.set_log_level(args.level)

    if(len(args.bg_mesh) > 0):
        tetra.load_sizing_field(args.bg_mesh)

    if len(args.csg) > 0:
        with open(args.csg, "r") as f:
            data = json.load(f)
        tetra.load_csg_tree(data)
    else:
        tetra.load_mesh(args.input)
    tetra.tetrahedralize()

    tetra.save(args.output,
               all_mesh=args.disable_filtering,
               smooth_open_boundary=args.smooth_open_boundary,
               floodfill=args.use_floodfill,
               use_input_for_wn=args.use_input_for_wn,
               manifold_surface=args.manifold_surface,
               binary=not args.no_binary)
