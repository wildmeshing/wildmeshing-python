#include "tetrahedralize.hpp"

#include <floattetwild/Mesh.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/FloatTetCutting.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/AABBWrapper.h>
#include <floattetwild/Statistics.h>
#include <floattetwild/TriangleInsertion.h>

#include <floattetwild/Logger.hpp>
#include <Eigen/Dense>

#include <igl/Timer.h>
#include <igl/writeSTL.h>

#include <Eigen/Dense>

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <thread>
#endif

namespace wildmeshing_binding
{

void tetrahedralize(py::module &m)
{
    m.def("tetrahedralize", [](const std::string &input, const std::string &output, double stop_quality, int max_its, int stage, int stop_p, double epsilon, double edge_length_r, bool mute_log, bool skip_simplify) {
        using namespace floatTetWild;
        using namespace Eigen;

        wildmeshing_binding::init_globals();

        int boolean_op = -1;

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

        unsigned int max_threads = std::numeric_limits<unsigned int>::max();

#ifdef FLOAT_TETWILD_USE_TBB
        const size_t MB = 1024 * 1024;
        const size_t stack_size = 64 * MB;
        unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
        num_threads = std::min(max_threads, num_threads);
        params.num_threads = num_threads;
        std::cout << "TBB threads " << num_threads << std::endl;
        tbb::task_scheduler_init scheduler(num_threads, stack_size);
#endif

        //    if(params.is_quiet){
        //        std::streambuf *orig_buf = cout.rdbuf();
        //        cout.rdbuf(NULL);
        //    }

        static bool initialized = false;
        if (!initialized)
        {
            Logger::init(!params.is_quiet, params.log_path);
            initialized = true;
        }

        params.log_level = std::max(0, std::min(6, params.log_level));
        spdlog::set_level(static_cast<spdlog::level::level_enum>(params.log_level));
        spdlog::flush_every(std::chrono::seconds(3));

        if (params.output_path.empty())
            params.output_path = params.input_path;
        if (params.log_path.empty())
            params.log_path = params.output_path;

        std::string output_mesh_name = params.output_path;
        if (params.output_path.size() > 3 && params.output_path.substr(params.output_path.size() - 3, params.output_path.size()) == "msh")
            output_mesh_name = params.output_path;
        else if (params.output_path.size() > 4 && params.output_path.substr(params.output_path.size() - 4, params.output_path.size()) == "mesh")
            output_mesh_name = params.output_path;
        else
            output_mesh_name = params.output_path + "_" + params.postfix + ".msh";

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        std::vector<Vector3> input_vertices;
        std::vector<Vector3i> input_faces;
        std::vector<int> input_tags;

        if (!params.tag_path.empty())
        {
            input_tags.reserve(input_faces.size());
            std::string line;
            std::ifstream fin(params.tag_path);
            if (fin.is_open())
            {
                while (getline(fin, line))
                {
                    input_tags.push_back(std::stoi(line));
                }
                fin.close();
            }
        }

        igl::Timer timer;

        GEO::Mesh sf_mesh;
        if (!MeshIO::load_mesh(params.input_path, input_vertices, input_faces, sf_mesh, input_tags))
        {
            logger().error("Unable to load mesh at {}", params.input_path);
            MeshIO::write_mesh(output_mesh_name, mesh, false);
            return false;
        }
        else if (input_vertices.empty() || input_faces.empty())
        {
            MeshIO::write_mesh(output_mesh_name, mesh, false);
            return false;
        }

        if (input_tags.size() != input_faces.size())
        {
            input_tags.resize(input_faces.size());
            std::fill(input_tags.begin(), input_tags.end(), 0);
        }
        AABBWrapper tree(sf_mesh);

        if (!params.init(tree.get_sf_diag()))
        {
            return false;
        }

        stats().record(StateInfo::init_id, 0, input_vertices.size(), input_faces.size(), -1, -1);

        timer.start();
        simplify(input_vertices, input_faces, input_tags, tree, params, skip_simplify);
        tree.init_b_mesh_and_tree(input_vertices, input_faces);
        logger().info("preprocessing {}s", timer.getElapsedTimeInSec());
        logger().info("");
        stats().record(StateInfo::preprocessing_id, timer.getElapsedTimeInSec(), input_vertices.size(),
                       input_faces.size(), -1, -1);
        if (params.log_level <= 1)
            output_component(input_vertices, input_faces, input_tags);

        timer.start();
        std::vector<bool> is_face_inserted(input_faces.size(), false);
        FloatTetDelaunay::tetrahedralize(input_vertices, input_faces, tree, mesh, is_face_inserted);
        logger().info("#v = {}", mesh.get_v_num());
        logger().info("#t = {}", mesh.get_t_num());
        logger().info("tetrahedralizing {}s", timer.getElapsedTimeInSec());
        logger().info("");
        stats().record(StateInfo::tetrahedralization_id, timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
                       -1, -1);

        timer.start();
        insert_triangles(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, false);
        logger().info("cutting {}s", timer.getElapsedTimeInSec());
        logger().info("");
        stats().record(StateInfo::cutting_id, timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
                       mesh.get_max_energy(), mesh.get_avg_energy(),
                       std::count(is_face_inserted.begin(), is_face_inserted.end(), false));

        timer.start();
        optimization(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, {{1, 1, 1, 1}});
        logger().info("mesh optimization {}s", timer.getElapsedTimeInSec());
        logger().info("");
        stats().record(StateInfo::optimization_id, timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
                       mesh.get_max_energy(), mesh.get_avg_energy());

        timer.start();
        correct_tracked_surface_orientation(mesh, tree);
        logger().info("correct_tracked_surface_orientation done");
        if (boolean_op < 0)
            filter_outside(mesh);
        else
            boolean_operation(mesh, boolean_op);
        stats().record(StateInfo::wn_id, timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
                       mesh.get_max_energy(), mesh.get_avg_energy());
        logger().info("after winding number");
        logger().info("#v = {}", mesh.get_v_num());
        logger().info("#t = {}", mesh.get_t_num());
        logger().info("winding number {}s", timer.getElapsedTimeInSec());
        logger().info("");

        //    if (params.output_path.size() > 3
        //        && params.output_path.substr(params.output_path.size() - 3, params.output_path.size()) == "msh")
        //        MeshIO::write_mesh(params.output_path, mesh, false);
        //    else if (params.output_path.size() > 4
        //             && params.output_path.substr(params.output_path.size() - 4, params.output_path.size()) == "mesh")
        //        MeshIO::write_mesh(params.output_path, mesh, false);
        //    else
        //        MeshIO::write_mesh(params.output_path + "_" + params.postfix + ".msh", mesh, false);

        // std::ofstream fout(params.log_path + "_" + params.postfix + ".csv");
        // if (fout.good())
        //     fout << stats();
        // fout.close();

        if (!params.envelope_log.empty())
        {
            std::ofstream fout(params.envelope_log);
            fout << envelope_log_csv;
        }

        return true;
    },
          "Robust Tetrahedralization, this is an alpha developement version of TetWild. For a stable release refer to the C++ version https://github.com/Yixin-Hu/TetWild",
          py::arg("input"), // "Input surface mesh INPUT in .off/.obj/.stl/.ply format. (string, required)"
          // py::arg("postfix") = "" //"Add postfix into outputs' file name"
          py::arg("output") = "",              // "Output tetmesh OUTPUT in .msh format. (string, optional, default: input_file+postfix+'.msh')"
          py::arg("stop_quality") = 10,        // "Specify max AMIPS energy for stopping mesh optimization"
          py::arg("max_its") = 80,             // "Max number of mesh optimization iterations"
          py::arg("stage") = 2,                // "Specify envelope stage"
          py::arg("stop_p") = -1,              //
          py::arg("epsilon") = 1e-3,           // "relative envelope epsilon_r. Absolute epsilonn = epsilon_r * diagonal_of_bbox"
          py::arg("edge_length_r") = 1. / 20., // "Relative target edge length l_r. Absolute l = l_r * diagonal_of_bbox"
          py::arg("mute_log") = false,         // "Mute prints");
          py::arg("skip_simplify") = false     //
    );
}
}