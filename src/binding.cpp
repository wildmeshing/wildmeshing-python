#include "triwild/optimization.h"
#include "triwild/feature.h"
#include "triwild/triangulation.h"

#include "triwild/meshio.hpp"

#include <igl/Timer.h>
#include <igl/writeSTL.h>


#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>

namespace py = pybind11;




namespace {
    void init_globals()
    {
        static bool initialized = false;

        if(!initialized)
        {
#ifndef WIN32
            setenv("GEO_NO_SIGNAL_HANDLER", "1", 1);
#endif

            GEO::initialize();

            // Import standard command line arguments, and custom ones
            GEO::CmdLine::import_arg_group("standard");
            GEO::CmdLine::import_arg_group("pre");
            GEO::CmdLine::import_arg_group("algo");

            initialized = true;
        }
    }

}

PYBIND11_MODULE(wildmeshing, m) {
    // m.doc() = "..."


    m.def("triangulate", [](const std::string &input, const std::string &feature_input, const std::string &output,
        double stop_quality, int max_its, int stage,
        double epsilon, double feature_epsilon,
        double target_edge_len, double edge_length_r,
        double flat_feature_angle,
        bool cut_outside,
        bool skip_eps,
        const std::string &hole_file,
        bool mute_log) {
            init_globals();

            using namespace triwild;
            using namespace std;

            args.input = input;
            // args.postfix = postfix;
            args.feature_input = feature_input;
            args.stop_quality = stop_quality;
            args.max_its = max_its;
            args.stage = stage;
            args.epsilon = epsilon;
            args.feature_epsilon = feature_epsilon;
            args.target_edge_len = target_edge_len;
            args.edge_length_r = edge_length_r;
            args.flat_feature_angle = flat_feature_angle;

            args.mute_log = mute_log;

            // app.add_option("--log-file", args.log_file, "Output a log file.");

            // bool diffusion_curves = false;
            // app.add_flag("--diffusion-curves", diffusion_curves, "Export for diffusion curves");

            igl::Timer main_timer;
            double load_mesh_time = 0;
            double preprocess_time = 0;
            double bsp_time = 0;
            double optimization_time = 0;
            double curving_time = 0;
            double cut_and_hole_time = 0;


            if (args.mute_log) {
                std::streambuf *orig_buf = cout.rdbuf();
                cout.rdbuf(NULL);
            }

            // if (args.output == "")
                // args.output = args.input + "_" + args.postfix;

            igl::Timer igl_timer;
            double t;
            ///////////////////
            cout << "Loading and preprocessing ..." << endl;
            igl_timer.start();

            Eigen::MatrixXd V;
            std::vector<std::array<int, 2>> edges;
            triangulation::load_input(args.input, V, edges);
            int input_v_cnt = V.rows();
            int input_e_cnt = edges.size();
            load_mesh_time = igl_timer.getElapsedTime();
            igl_timer.start();

            if (feature::init(args.feature_input))
                args.is_preserving_feature = true;

            if (args.stop_quality < 0)
                args.stop_quality = args.is_preserving_feature ? 20 : 10;
            if (args.epsilon < 0)
                // args.epsilon = args.is_preserving_feature ? 5e-3 : 1e-3;
                args.epsilon = args.is_preserving_feature ? 2e-3 : 1e-3;

            double line_width = 0.3 * args.diagonal_len;
            double f_line_width = 0.5 * args.diagonal_len;
            double s_f_line_width = 0.4 * args.diagonal_len;
            // double draw_points = true;
            double draw_points = false;

            double point_size = 0.0005 * args.diagonal_len;
            double f_point_size = 0.001 * args.diagonal_len;
            double s_f_point_size = 0.003 * args.diagonal_len;

            std::string line_col = "0 0 0";
            std::string f_line_col = "0.9 0.3 0.2";
            std::string s_f_line_col = "0.1 0.7 0.6";

            std::string point_col = "0 0 0";
            std::string f_point_col = "0.9 0.3 0.2";
            std::string s_f_point_col = "0.1 0.7 0.6";

            GEO::Mesh b_mesh;
            triangulation::preprocessing(V, edges, b_mesh);
            GEO::MeshFacetsAABB b_tree(b_mesh);
            t = igl_timer.getElapsedTime();
            preprocess_time = t;
            cout << "Loaded and preprocessed." << endl;
            cout << "time = " << t << "s" << endl << endl;

            ///////////////////
            cout << "BSP subdivision..." << endl;
            igl_timer.start();
            MeshData mesh;
            std::vector<std::vector<int>> tag_boundary_es;
            triangulation::BSP_subdivision(V, edges, mesh, tag_boundary_es);
            t = igl_timer.getElapsedTime();
            bsp_time = t;
            cout << "BSP subdivision done." << endl;
            cout << "time = " << t << "s" << endl << endl;

            ///////////////////
            cout << "Mesh optimization..." << endl;
            igl_timer.start();
            optimization::init(V, edges, tag_boundary_es, mesh, b_tree);

            optimization::refine(mesh, b_tree, std::array<int, 4>({1, 1, 1, 1}));
            t = igl_timer.getElapsedTime();
            optimization_time = t;
            cout << "Mesh optimization done." << endl;
            cout << "time = " << t << "s" << endl;

            // if (!skip_eps) {
            //     export_eps(mesh,
            //         line_width, line_col, point_size, point_col,
            //         f_line_width, f_line_col, f_point_size, f_point_col,
            //         s_f_line_width, s_f_line_col, s_f_point_size, s_f_point_col,
            //         draw_points, args.output + "_lin.eps");
            // }

            ///////////////////
            if (!args.is_preserving_feature) {
                optimization::output_mesh(mesh);

                if (args.log_file != "") {
                    std::ofstream file(args.log_file);
                    if (!file.fail()) {
                        file << "load_mesh_time " << load_mesh_time << "\n";
                        file << "preprocess_time " << preprocess_time << "\n";
                        file << "bsp_time " << bsp_time << "\n";
                        file << "optimization_time " << optimization_time << "\n";
                        file << "curving_time " << curving_time << "\n";
                        file << "cut_and_hole_time " << cut_and_hole_time << "\n";
                        file << "input_v_cnt " << input_v_cnt  << "\n";
                        file << "input_e_cnt " << input_e_cnt  << "\n";
                        optimization::output_stats(mesh, file);
                    }
                    file.close();
                }

                return;
            }

            cout << "Curving..." << endl;
            igl_timer.start();
            feature::curving(mesh, b_tree);
            t = igl_timer.getElapsedTime();
            curving_time = t;
            cout << "Curving done." << endl;
            cout << "time = " << t << "s" << endl;

            igl_timer.start();
            if (cut_outside)
                optimization::erase_outside(mesh);
            if (hole_file != "")
                optimization::erase_holes(mesh, hole_file);
            cut_and_hole_time = igl_timer.getElapsedTime();

            optimization::output_mesh(mesh);

            // if (diffusion_curves)
                // write_msh_DiffusionCurve(mesh, args.output + "DC");

            if (!skip_eps) {
               export_eps(mesh,
                           line_width, line_col, point_size, point_col,
                           f_line_width, f_line_col, f_point_size, f_point_col,
                           s_f_line_width, s_f_line_col, s_f_point_size, s_f_point_col,
                           draw_points, args.output + ".eps");

                // const double N = 10;
                // for(int i = 0; i <= N; ++i)
                // {
                // const double t = i/N;
                //     const double t = i/N;
                //     export_eps(mesh,
                //      line_width, line_col, point_size, point_col,
                //      f_line_width, f_line_col, f_point_size, f_point_col,
                //      s_f_line_width, s_f_line_col, s_f_point_size, s_f_point_col,
                //      draw_points, args.output + "_" + std::to_string(i) + ".eps", t);
                // }
            }

            if (args.log_file != "") {
                std::ofstream file(args.log_file);
                if (!file.fail()) {
                    file << "load_mesh_time " << load_mesh_time << "\n";
                    file << "preprocess_time " << preprocess_time << "\n";
                    file << "bsp_time " << bsp_time << "\n";
                    file << "optimization_time " << optimization_time << "\n";
                    file << "curving_time " << curving_time << "\n";
                    file << "cut_and_hole_time " << cut_and_hole_time << "\n";
                    file << "input_v_cnt " << input_v_cnt  << "\n";
                    file << "input_e_cnt " << input_e_cnt  << "\n";
                    optimization::output_stats(mesh, file);
                    feature::output_stats(mesh, file);
                }
                file.close();
            }
    },
    "Robust Triangulation",
    py::arg("input"), // "Input segments in .obj format"
    // py::arg("postfix") = "" //"Add postfix into outputs' file name"
    py::arg("feature_input") = "", // "Input feature json file"
    py::arg("output") = "", //"Output path"
    py::arg("stop_quality") = -1, // "Specify max AMIPS energy for stopping mesh optimization"
    py::arg("max_its") = 80, // "Max number of mesh optimization iterations"
    py::arg("stage") = 1, // "Specify envelope stage"
    py::arg("epsilon") = -1, // "relative envelope epsilon_r. Absolute epsilonn = epsilon_r * diagonal_of_bbox"
    py::arg("feature_epsilon") = 1e-3, // "Relative feature envelope mu_r. Absolute mu = mu_r * diagonal_of_bbox"
    py::arg("target_edge_len") = -1, // "Absolute target edge length l"
    py::arg("edge_length_r") = 1./20., // "Relative target edge length l_r. Absolute l = l_r * diagonal_of_bbox"
    py::arg("flat_feature_angle") = 10., // "Desired minimal angle"
    py::arg("cut_outside") = false, // "Remove \"outside part\""
    py::arg("skip_eps") = false, // "Skip saving eps"
    py::arg("hole_file") = "", // "Input a .xyz file for specifying points inside holes you want to remove"
    py::arg("mute_log") = false // "Mute prints");
    );
}