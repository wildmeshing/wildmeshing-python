#include "triwild/optimization.h"
#include "triwild/feature.h"
#include "triwild/triangulation.h"
#include "triwild/do_triwild.h"

#include "triwild/meshio.hpp"

#include <igl/Timer.h>
#include <igl/writeSTL.h>




#ifdef USE_TBB
#include <tbb/task_scheduler_init.h>
#include <thread>
#endif

#include <floattetwild/Mesh.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/FloatTetCutting.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/AABBWrapper.h>
#include <floattetwild/Statistics.h>

#include <floattetwild/Logger.hpp>
#include <Eigen/Dense>



#include <geogram/mesh/mesh.h>
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

m.def("triangulate_data", [](const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const py::object &feature_info,
        double stop_quality, int max_its, int stage,
        double epsilon, double feature_epsilon,
        double target_edge_len, double edge_length_r,
        double flat_feature_angle,
        bool cut_outside,
        bool skip_eps,
        const std::string &hole_file,
        bool mute_log) {
            json jfeature_info;
            if(!feature_info.is(py::none())){
                const std::string json_string = py::str(feature_info);
			    jfeature_info = json::parse(json_string);
            }

            Eigen::MatrixXd V_out;
            Eigen::MatrixXi F_out;
            Eigen::MatrixXd nodes;
            std::vector<std::vector<int>> F_nodes;

            triwild::do_triwild(V, E, jfeature_info,
                    V_out, F_out, nodes, F_nodes,
                     stop_quality, max_its, stage,
                     epsilon,  feature_epsilon,
                     target_edge_len,  edge_length_r,
                     flat_feature_angle,
                    cut_outside,
                    hole_file,
                    mute_log);
            return py::make_tuple(V_out, F_out, nodes, F_nodes);
        },
    "Robust Triangulation",
    py::arg("V"), // "Input vertices"
    py::arg("E"), // "Input edges"
    py::arg("feature_info") = py::none(), // "Json string containing the features"
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
            args.output = output;
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

#ifndef WIN32
            if (args.mute_log) {
                std::streambuf *orig_buf = cout.rdbuf();
                cout.rdbuf(NULL);
            }
#endif

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
    )


    .def("tetrahedralize", [](const std::string &input, const std::string &output,
        double stop_quality, int max_its, int stage, int stop_p,
        double epsilon,
        double edge_length_r,
        bool mute_log,
        bool skip_simplify
        ) {
            using namespace floatTetWild;
            using namespace Eigen;
            init_globals();



            Mesh mesh;
            Parameters &params = mesh.params;

            params.input_path = input;
            params.output_path = output;

            params.stop_energy = stop_quality;
            params.max_its = max_its;
            params.stage = stage;
            params.stop_p = stop_p;

            params.eps_rel = epsilon;
            params.ideal_edge_length_rel = edge_length_r;

            params.is_quiet = mute_log;


            // command_line.add_option("--tag", params.tag_path, "");
            // const int UNION = 0;
            // const int INTERSECTION = 1;
            // const int DIFFERENCE = 2;
            int boolean_op = -1;
            // command_line.add_option("--op", boolean_op, "");
            // command_line.add_option("--postfix", params.postfix, "");
            params.log_level = 2;
            // command_line.add_option("--log", params.log_path, "Log info to given file.");
            // command_line.add_option("--level", params.log_level, "Log level (0 = most verbose, 6 = off).");

            unsigned int max_threads = std::numeric_limits<unsigned int>::max();
            // #ifdef USE_TBB
            //         command_line.add_option("--max-threads", max_threads, "maximum number of threads used");
            // #endif

            #ifdef USE_TBB
            const size_t MB = 1024 * 1024;
            const size_t stack_size = 64 * MB;
            unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
            num_threads = std::min(max_threads, num_threads);
            params.num_threads = num_threads;
            std::cout << "TBB threads " << num_threads << std::endl;
            tbb::task_scheduler_init scheduler(num_threads, stack_size);
            #endif

            Logger::init(!params.is_quiet, params.log_path);
            params.log_level = std::max(0, std::min(6, params.log_level));
            spdlog::set_level(static_cast<spdlog::level::level_enum>(params.log_level));
            spdlog::flush_every(std::chrono::seconds(3));

            // GEO::Logger *geo_logger = GEO::Logger::instance();
            // geo_logger->unregister_all_clients();
            // geo_logger->register_client(new GeoLoggerForward(logger().clone("geogram")));
            // geo_logger->set_pretty(false);


            if (params.output_path.empty())
                params.output_path = params.input_path;
            if (params.log_path.empty())
                params.log_path = params.output_path;


            std::string output_mesh_name = params.output_path;
            if (params.output_path.size() > 3
                && params.output_path.substr(params.output_path.size() - 3, params.output_path.size()) == "msh")
                output_mesh_name = params.output_path;
            else if (params.output_path.size() > 4
                     && params.output_path.substr(params.output_path.size() - 4, params.output_path.size()) == "mesh")
                output_mesh_name = params.output_path;
            else
                output_mesh_name = params.output_path + "_" + params.postfix + ".msh";

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            std::vector<Vector3> input_vertices;
            std::vector<floatTetWild::Vector3i> input_faces;
            std::vector<int> input_tags;

            if(!params.tag_path.empty()) {
                input_tags.reserve(input_faces.size());
                std::string line;
                std::ifstream fin(params.tag_path);
                if (fin.is_open()) {
                    while (getline(fin, line)) {
                        input_tags.push_back(std::stoi(line));
                    }
                    fin.close();
                }
            }


            igl::Timer timer;

            GEO::Mesh sf_mesh;
            if (!MeshIO::load_mesh(params.input_path, input_vertices, input_faces, sf_mesh, input_tags)) {
                logger().error("Unable to load mesh at {}", params.input_path);
                MeshIO::write_mesh(output_mesh_name, mesh, false);
                return;
            } else if(input_vertices.empty() || input_faces.empty()){
                MeshIO::write_mesh(output_mesh_name, mesh, false);
                return;
            }

            if (input_tags.size() != input_faces.size()) {
                input_tags.resize(input_faces.size());
                std::fill(input_tags.begin(), input_tags.end(), 0);
            }
            AABBWrapper tree(sf_mesh);

            if (!params.init(tree.get_sf_diag())) {
                return;
            }

            Statistics::stats().states.push_back(StateInfo(StateInfo::init_id,
                                                       0, input_vertices.size(), input_faces.size(), -1, -1));

            timer.start();
            simplify(input_vertices, input_faces, input_tags, tree, params, skip_simplify);
            tree.init_b_mesh_and_tree(input_vertices, input_faces);
            logger().info("preprocessing {}s", timer.getElapsedTimeInSec());
            logger().info("");
            Statistics::stats().states.push_back(StateInfo(StateInfo::preprocessing_id,
                                                           timer.getElapsedTimeInSec(), input_vertices.size(),
                                                           input_faces.size(),
                                                           -1, -1));
            if(params.log_level<=1)
                output_component(input_vertices, input_faces, input_tags);

            timer.start();
            std::vector<bool> is_face_inserted(input_faces.size(), false);
            FloatTetDelaunay::tetrahedralize(input_vertices, input_faces, tree, mesh, is_face_inserted);
            logger().info("#v = {}", mesh.get_v_num());
            logger().info("#t = {}", mesh.get_t_num());
            logger().info("tetrahedralizing {}s", timer.getElapsedTimeInSec());
            logger().info("");
            Statistics::stats().states.push_back(StateInfo(StateInfo::tetrahedralization_id,
                                                           timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
                                                           -1, -1));

            timer.start();
            cutting(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree);
            logger().info("cutting {}s", timer.getElapsedTimeInSec());
            logger().info("");
            Statistics::stats().states.push_back(StateInfo(StateInfo::cutting_id,
                                                           timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
                                                           mesh.get_max_energy(), mesh.get_avg_energy(),
                                                           std::count(is_face_inserted.begin(), is_face_inserted.end(), false)));

            timer.start();
            optimization(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, {{1, 1, 1, 1}});
            logger().info("mesh optimization {}s", timer.getElapsedTimeInSec());
            logger().info("");
            Statistics::stats().states.push_back(StateInfo(StateInfo::optimization_id,
                                                           timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
                                                           mesh.get_max_energy(), mesh.get_avg_energy()));

            timer.start();
            if(boolean_op<0)
                filter_outside(mesh);
            else
                boolean_operation(mesh, boolean_op);
            Statistics::stats().states.push_back(StateInfo(StateInfo::wn_id,
                                                           timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
                                                           mesh.get_max_energy(), mesh.get_avg_energy()));
            logger().info("after winding number");
            logger().info("#v = {}", mesh.get_v_num());
            logger().info("#t = {}", mesh.get_t_num());
            logger().info("winding number {}s", timer.getElapsedTimeInSec());
            logger().info("");


               // if (params.output_path.size() > 3
               //     && params.output_path.substr(params.output_path.size() - 3, params.output_path.size()) == "msh")
               //     MeshIO::write_mesh(params.output_path, mesh, false);
               // else if (params.output_path.size() > 4
               //          && params.output_path.substr(params.output_path.size() - 4, params.output_path.size()) == "mesh")
               //     MeshIO::write_mesh(params.output_path, mesh, false);
               // else
               //     MeshIO::write_mesh(params.output_path + "_" + params.postfix + ".msh", mesh, false);

            MeshIO::write_mesh(output_mesh_name, mesh, false);
            MeshIO::write_surface_mesh(params.output_path + "_" + params.postfix + "_sf.obj", mesh, false);

            std::ofstream fout(params.log_path + "_" + params.postfix + ".csv");
            if (fout.good())
                fout << Statistics::stats();
            fout.close();
    },
    "Robust Tetrahedralization, this is an alpha developement version of TetWild. For a stable release refer to the C++ version https://github.com/Yixin-Hu/TetWild",
    py::arg("input"), // "Input surface mesh INPUT in .off/.obj/.stl/.ply format. (string, required)"
    // py::arg("postfix") = "" //"Add postfix into outputs' file name"
    py::arg("output") = "", // "Output tetmesh OUTPUT in .msh format. (string, optional, default: input_file+postfix+'.msh')"
    py::arg("stop_quality") = 10, // "Specify max AMIPS energy for stopping mesh optimization"
    py::arg("max_its") = 80, // "Max number of mesh optimization iterations"
    py::arg("stage") = 2, // "Specify envelope stage"
    py::arg("stop_p") = -1, //
    py::arg("epsilon") = 1e-3, // "relative envelope epsilon_r. Absolute epsilonn = epsilon_r * diagonal_of_bbox"
    py::arg("edge_length_r") = 1./20., // "Relative target edge length l_r. Absolute l = l_r * diagonal_of_bbox"
    py::arg("mute_log") = false, // "Mute prints");
    py::arg("skip_simplify") = false //
    );

}