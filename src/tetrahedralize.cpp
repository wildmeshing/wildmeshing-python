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

#include <memory>

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <thread>
#endif

using namespace floatTetWild;
using namespace Eigen;

namespace wildmeshing_binding
{
    namespace {
        class Tetrahedralizer
        {
            public:
                Mesh mesh;
                int boolean_op = -1;

                std::vector<Vector3> input_vertices;
                std::vector<Vector3i> input_faces;
                std::vector<int> input_tags;
                GEO::Mesh sf_mesh;
                std::unique_ptr<AABBWrapper> tree;
                bool skip_simplify;

                Tetrahedralizer(
                    double stop_quality, int max_its, int stage, int stop_p,
                double epsilon, double edge_length_r,
                bool skip_simplify, bool smooth_open_boundary, bool manifold_surface, bool correct_surface_orientation):
                skip_simplify(skip_simplify)
                {
                    Parameters &params = mesh.params;

                    // params.input_path = input;
                    // params.output_path = output;

                    params.stop_energy = stop_quality;
                    params.max_its = max_its;
                    params.stage = stage;
                    params.stop_p = stop_p;

                    params.eps_rel = epsilon;
                    params.ideal_edge_length_rel = edge_length_r;

                    params.correct_surface_orientation = correct_surface_orientation;
                    params.smooth_open_boundary = smooth_open_boundary;
                    params.manifold_surface = manifold_surface;

                    params.log_level = 0;
                }

                void set_num_threads(int num_threads)
                {
                    Parameters &params = mesh.params;
                    params.num_threads = num_threads;
                }

                void set_log_level(int level)
                {
                    Parameters &params = mesh.params;
                    params.log_level = std::max(0, std::min(6, level));
                    spdlog::set_level(static_cast<spdlog::level::level_enum>(params.log_level));
                    spdlog::flush_every(std::chrono::seconds(3));
                }

                bool load_mesh(const std::string &path, const std::string &tag_path = "")
                {
                    Parameters &params = mesh.params;
                    params.input_path = path;
                    params.tag_path = tag_path;
                    igl::Timer timer;

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

                    if (!MeshIO::load_mesh(params.input_path, input_vertices, input_faces, sf_mesh, input_tags))
                    {
                        logger().error("Unable to load mesh at {}", params.input_path);
                        return false;
                    }
                    else if (input_vertices.empty() || input_faces.empty())
                    {
                        return false;
                    }

                    if (input_tags.size() != input_faces.size())
                    {
                        input_tags.resize(input_faces.size());
                        std::fill(input_tags.begin(), input_tags.end(), 0);
                    }
                    tree = std::make_unique<AABBWrapper>(sf_mesh);

                    if (!params.init(tree->get_sf_diag()))
                    {
                        return false;
                    }

                    stats().record(StateInfo::init_id, 0, input_vertices.size(), input_faces.size(), -1, -1);
                    return true;
                }

                void tetrahedralize()
                {
                    Parameters &params = mesh.params;
                    igl::Timer timer;

                    timer.start();
                    simplify(input_vertices, input_faces, input_tags, *tree, params, skip_simplify);
                    tree->init_b_mesh_and_tree(input_vertices, input_faces);
                    logger().info("preprocessing {}s", timer.getElapsedTimeInSec());
                    logger().info("");
                    stats().record(StateInfo::preprocessing_id, timer.getElapsedTimeInSec(), input_vertices.size(),
                                   input_faces.size(), -1, -1);
                    if (params.log_level <= 1)
                        output_component(input_vertices, input_faces, input_tags);

                    timer.start();
                    std::vector<bool> is_face_inserted(input_faces.size(), false);
                    FloatTetDelaunay::tetrahedralize(input_vertices, input_faces, *tree, mesh, is_face_inserted);
                    logger().info("#v = {}", mesh.get_v_num());
                    logger().info("#t = {}", mesh.get_t_num());
                    logger().info("tetrahedralizing {}s", timer.getElapsedTimeInSec());
                    logger().info("");
                    stats().record(StateInfo::tetrahedralization_id, timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
                                   -1, -1);

                    timer.start();
                    insert_triangles(input_vertices, input_faces, input_tags, mesh, is_face_inserted, *tree, false);
                    logger().info("cutting {}s", timer.getElapsedTimeInSec());
                    logger().info("");
                    stats().record(StateInfo::cutting_id, timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
                                   mesh.get_max_energy(), mesh.get_avg_energy(),
                                   std::count(is_face_inserted.begin(), is_face_inserted.end(), false));

                    timer.start();
                    optimization(input_vertices, input_faces, input_tags, is_face_inserted, mesh, *tree, {{1, 1, 1, 1}});
                    logger().info("mesh optimization {}s", timer.getElapsedTimeInSec());
                    logger().info("");
                    stats().record(StateInfo::optimization_id, timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
                                   mesh.get_max_energy(), mesh.get_avg_energy());

                    timer.start();
                    correct_tracked_surface_orientation(mesh, *tree);
                    logger().info("correct_tracked_surface_orientation done");
                }

                void save(const std::string &path)
                {
                    igl::Timer timer;

                    Parameters &params = mesh.params;
                    params.output_path = path;

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

                    if (boolean_op < 0)
                    {
                        if (params.smooth_open_boundary)
                        {
                            floatTetWild::smooth_open_boundary(mesh, *tree);
                            for (auto &t : mesh.tets)
                            {
                                if (t.is_outside)
                                    t.is_removed = true;
                            }
                        }
                        else
                            filter_outside(mesh);
                    }
                    else
                        boolean_operation(mesh, boolean_op);
                    if (params.manifold_surface)
                    {
                        //        MeshIO::write_mesh(params.output_path + "_" + params.postfix + "_non_manifold.msh", mesh, false);
                        floatTetWild::manifold_surface(mesh);
                    }
                    stats().record(StateInfo::wn_id, timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
                                   mesh.get_max_energy(), mesh.get_avg_energy());
                    logger().info("after winding number");
                    logger().info("#v = {}", mesh.get_v_num());
                    logger().info("#t = {}", mesh.get_t_num());
                    logger().info("winding number {}s", timer.getElapsedTimeInSec());
                    logger().info("");

                    if (params.output_path.size() > 3 && params.output_path.substr(params.output_path.size() - 3, params.output_path.size()) == "msh")
                        MeshIO::write_mesh(params.output_path, mesh, false);
                    else if (params.output_path.size() > 4 && params.output_path.substr(params.output_path.size() - 4, params.output_path.size()) == "mesh")
                        MeshIO::write_mesh(params.output_path, mesh, false);
                    else
                        MeshIO::write_mesh(params.output_path + "_" + params.postfix + ".msh", mesh, false);
                }

                std::string get_log() {
                    std::stringstream ss;
                    ss << envelope_log_csv;
                    return ss.str();
                }
        };

    }

void tetrahedralize(py::module &m)
{
    m.def("tetrahedralize", [](const std::string &input, const std::string &output, double stop_quality, int max_its, int stage, int stop_p, double epsilon, double edge_length_r, bool mute_log, bool skip_simplify, bool smooth_open_boundary, bool manifold_surface, bool correct_surface_orientation) {
        wildmeshing_binding::init_globals();

        unsigned int max_threads = std::numeric_limits<unsigned int>::max();
        unsigned int num_threads = 1;
#ifdef FLOAT_TETWILD_USE_TBB
            const size_t MB = 1024 * 1024;
        const size_t stack_size = 64 * MB;
        num_threads = std::max(1u, std::thread::hardware_concurrency());
        num_threads = std::min(max_threads, num_threads);
        // params.num_threads = num_threads;
        std::cout << "TBB threads " << num_threads << std::endl;
        tbb::task_scheduler_init scheduler(num_threads, stack_size);
#endif

        static bool initialized = false;
        if (!initialized)
        {
            Logger::init(!mute_log);
            initialized = true;
        }

        Tetrahedralizer tetra(stop_quality, max_its, stage, stop_p, epsilon, edge_length_r, skip_simplify, smooth_open_boundary, manifold_surface, correct_surface_orientation);
        tetra.set_num_threads(num_threads);
        if(!tetra.load_mesh(input))
            return false;

        tetra.tetrahedralize();
        tetra.save(output);

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
          py::arg("skip_simplify") = false,    //
          py::arg("smooth_open_boundary") = false,
          py::arg("manifold_surface") = false,
          py::arg("correct_surface_orientation") = false
    );
}
}