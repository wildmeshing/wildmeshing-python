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
#include <floattetwild/CSGTreeParser.hpp>

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
namespace
{
class Tetrahedralizer
{
public:
    Mesh mesh;

    std::vector<Vector3> input_vertices;
    std::vector<Vector3i> input_faces;
    std::vector<int> input_tags;
    GEO::Mesh sf_mesh;
    std::unique_ptr<AABBWrapper> tree;
    bool skip_simplify;
    json tree_with_ids;
    bool has_json_csg = false;

    Tetrahedralizer(
        double stop_quality, int max_its, int stage, int stop_p,
        double epsilon, double edge_length_r,
        bool skip_simplify) : skip_simplify(skip_simplify)
    {
        wildmeshing_binding::init_globals();

        Parameters &params = mesh.params;

        // params.input_path = input;
        // params.output_path = output;

        params.stop_energy = stop_quality;
        params.max_its = max_its;
        params.stage = stage;
        params.stop_p = stop_p;

        params.eps_rel = epsilon;
        params.ideal_edge_length_rel = edge_length_r;

        params.log_level = 6;

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
        set_num_threads(num_threads);
    }

private:
    void set_num_threads(int num_threads)
    {
        Parameters &params = mesh.params;
        params.num_threads = num_threads;
    }

public:
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
            else
            {
                throw std::invalid_argument("Invalid mesh format tag at " + params.tag_path);
            }

        }

        if (!MeshIO::load_mesh(params.input_path, input_vertices, input_faces, sf_mesh, input_tags))
        {
            throw std::invalid_argument("Invalid mesh path at " + params.input_path);
            return false;
        }
        else if (input_vertices.empty() || input_faces.empty())
        {
            throw std::invalid_argument("Invalid mesh path at " + params.input_path);
            return false;
        }

        return load_mesh_aux();
    }

    bool boolean_operation(const std::string &json_string)
    {
        json csg_tree;
        std::ifstream file(json_string);

        if (file.is_open())
        {
            file >> csg_tree;
            file.close();
        }
        else
            csg_tree = json::parse(json_string);

        std::vector<std::string> meshes;

        CSGTreeParser::get_meshes(csg_tree, meshes, tree_with_ids);
        has_json_csg = true;

        bool ok = CSGTreeParser::load_and_merge(meshes, input_vertices, input_faces, sf_mesh, input_tags);
        if(!ok)
        {
            throw std::invalid_argument("Invalid mesh path in the json");
        }

        load_mesh_aux();
        return ok;
    }

    void set_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
    {
        if (F.cols() != 3)
            throw std::invalid_argument("Mesh format not supported, F should have 3 cols");
        if (V.cols() != 3)
            throw std::invalid_argument("Mesh format not supported, V should have 3 cols");

        input_vertices.resize(V.rows());
        for (int i = 0; i < V.rows(); ++i)
        {
            input_vertices[i][0] = V(i, 0);
            input_vertices[i][1] = V(i, 1);
            input_vertices[i][2] = V(i, 2);
        }

        input_faces.resize(F.rows());
        for (int i = 0; i < F.rows(); ++i)
        {
            input_faces[i][0] = F(i, 0);
            input_faces[i][1] = F(i, 1);
            input_faces[i][2] = F(i, 2);
        }

        MeshIO::load_mesh(input_vertices, input_faces, sf_mesh, input_tags);

        load_mesh_aux();
    }

    void set_meshes(const std::vector<Eigen::MatrixXd> &V, const std::vector<Eigen::MatrixXi> &F)
    {
        std::vector<std::vector<Vector3>> vs(V.size());
        std::vector<std::vector<Vector3i>> fs(F.size());

        if(V.size() != F.size())
            throw std::invalid_argument("V and F must have the same size");

        for (int j = 0; j < V.size(); ++j)
        {
            if (V[j].cols() != 3)
                throw std::invalid_argument("Mesh format not supported, V should have 3 cols");

            vs[j].resize(V[j].rows());
            for (int i = 0; i < V[j].rows(); ++i)
            {
                vs[j][i][0] = V[j](i, 0);
                vs[j][i][1] = V[j](i, 1);
                vs[j][i][2] = V[j](i, 2);
            }
        }

        for (int j = 0; j < F.size(); ++j){
            if (F[j].cols() != 3)
                throw std::invalid_argument("Mesh format not supported, F should have 3 cols");

            fs[j].resize(F[j].rows());
            for (int i = 0; i < F[j].rows(); ++i)
            {
                fs[j][i][0] = F[j](i, 0);
                fs[j][i][1] = F[j](i, 1);
                fs[j][i][2] = F[j](i, 2);
            }
        }

        CSGTreeParser::merge(vs, fs, input_vertices, input_faces, sf_mesh, input_tags);

        load_mesh_aux();
    }

private:
    bool load_mesh_aux()
    {
        Parameters &params = mesh.params;

        if (input_tags.size() != input_faces.size())
        {
            input_tags.resize(input_faces.size());
            std::fill(input_tags.begin(), input_tags.end(), 0);
        }
        tree = std::make_unique<AABBWrapper>(sf_mesh);

        if (!params.init(tree->get_sf_diag()))
        {
            throw std::invalid_argument("Unable to initialize the tree, probably a problem with the mesh");
            return false;
        }

        stats().record(StateInfo::init_id, 0, input_vertices.size(), input_faces.size(), -1, -1);
        return true;
    }

public:
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

    void save(const std::string &path, bool smooth_open_boundary, bool manifold_surface, bool correct_surface_orientation, int boolean_op = -1)
    {
        igl::Timer timer;

        Mesh mesh_copy = mesh;

        Parameters &params = mesh_copy.params;
        params.output_path = path;
        params.correct_surface_orientation = correct_surface_orientation;
        params.smooth_open_boundary = smooth_open_boundary;
        params.manifold_surface = manifold_surface;

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

        if (has_json_csg)
            floatTetWild::boolean_operation(mesh, tree_with_ids);
        else if (boolean_op >= 0)
            floatTetWild::boolean_operation(mesh, boolean_op);
        else
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
        if (params.manifold_surface)
        {
            floatTetWild::manifold_surface(mesh_copy);
        }

        stats().record(StateInfo::wn_id, timer.getElapsedTimeInSec(), mesh_copy.get_v_num(), mesh_copy.get_t_num(),
                       mesh_copy.get_max_energy(), mesh_copy.get_avg_energy());
        logger().info("after winding number");
        logger().info("#v = {}", mesh_copy.get_v_num());
        logger().info("#t = {}", mesh_copy.get_t_num());
        logger().info("winding number {}s", timer.getElapsedTimeInSec());
        logger().info("");

        if (params.output_path.size() > 3 && params.output_path.substr(params.output_path.size() - 3, params.output_path.size()) == "msh")
            MeshIO::write_mesh(params.output_path, mesh_copy, false);
        else if (params.output_path.size() > 4 && params.output_path.substr(params.output_path.size() - 4, params.output_path.size()) == "mesh")
            MeshIO::write_mesh(params.output_path, mesh_copy, false);
        else
            MeshIO::write_mesh(params.output_path + "_" + params.postfix + ".msh", mesh_copy, false);
    }

    void get_tet_mesh(bool smooth_open_boundary, bool manifold_surface, bool correct_surface_orientation, Eigen::MatrixXd &V, Eigen::MatrixXi &T, int boolean_op = -1)
    {
        igl::Timer timer;

        Mesh mesh_copy = mesh;

        const auto skip_tet    = [&mesh_copy](const int i) { return mesh_copy.tets[i].is_removed; };
        const auto skip_vertex = [&mesh_copy](const int i) { return mesh_copy.tet_vertices[i].is_removed; };
        std::vector<int> t_ids(mesh_copy.tets.size());
        std::iota(std::begin(t_ids), std::end(t_ids), 0);

        Parameters &params = mesh_copy.params;
        params.correct_surface_orientation = correct_surface_orientation;
        params.smooth_open_boundary = smooth_open_boundary;
        params.manifold_surface = manifold_surface;

        if (has_json_csg)
            floatTetWild::boolean_operation(mesh_copy, tree_with_ids);
        else if (boolean_op >= 0)
            floatTetWild::boolean_operation(mesh_copy, boolean_op);
        else
        {
            if (params.smooth_open_boundary)
            {
                floatTetWild::smooth_open_boundary(mesh_copy, *tree);
                for (auto &t : mesh_copy.tets)
                {
                    if (t.is_outside)
                        t.is_removed = true;
                }
            }
            else
                filter_outside(mesh_copy);
        }

        if (params.manifold_surface)
        {
            floatTetWild::manifold_surface(mesh_copy);
        }
        stats().record(StateInfo::wn_id, timer.getElapsedTimeInSec(), mesh_copy.get_v_num(), mesh_copy.get_t_num(),
                       mesh_copy.get_max_energy(), mesh_copy.get_avg_energy());
        logger().info("after winding number");
        logger().info("#v = {}", mesh_copy.get_v_num());
        logger().info("#t = {}", mesh_copy.get_t_num());
        logger().info("winding number {}s", timer.getElapsedTimeInSec());
        logger().info("");

        int cnt_v = 0;
        std::map<int, int> old_2_new;
        for (int i = 0; i < mesh_copy.tet_vertices.size(); i++)
        {
            if (!skip_vertex(i))
            {
                old_2_new[i] = cnt_v;
                cnt_v++;
            }
        }
        int cnt_t = 0;
        for (const int i : t_ids)
        {
            if (!skip_tet(i))
                cnt_t++;
        }

        V.resize(cnt_v, 3);
        int index = 0;
        for (size_t i = 0; i < mesh_copy.tet_vertices.size(); i++)
        {
            if (skip_vertex(i))
                continue;
            V.row(index++) << mesh_copy.tet_vertices[i][0], mesh_copy.tet_vertices[i][1], mesh_copy.tet_vertices[i][2];
        }

        T.resize(cnt_t, 4);
        index = 0;

        for (const int i : t_ids)
        {
            if (skip_tet(i))
                continue;
            for (int j = 0; j < 4; j++)
            {
                T(index, j) = old_2_new[mesh_copy.tets[i][j]];
            }
            index++;
        }
    }

    std::string get_stats() const
    {
        std::stringstream ss;
        ss << stats();
        return ss.str();
    }
};

} // namespace

void tetrahedralize(py::module &m)
{
    auto &tetra = py::class_<Tetrahedralizer>(m, "Tetrahedralizer")
                      .def(py::init<
                               double, int, int, int,
                               double, double,
                               bool>(),
                           py::arg("stop_quality") = 10,        // "Specify max AMIPS energy for stopping mesh optimization"
                           py::arg("max_its") = 80,             // "Max number of mesh optimization iterations"
                           py::arg("stage") = 2,                // "Specify envelope stage"
                           py::arg("stop_p") = -1,              //
                           py::arg("epsilon") = 1e-3,           // "relative envelope epsilon_r. Absolute epsilonn = epsilon_r * diagonal_of_bbox"
                           py::arg("edge_length_r") = 1. / 20., // "Relative target edge length l_r. Absolute l = l_r * diagonal_of_bbox"
                           py::arg("skip_simplify") = false     //
                           )

                      .def("set_log_level", [](Tetrahedralizer &t, int level) { t.set_log_level(level); }, "sets log level, valid value between 0 (all logs) and 6 (no logs)", py::arg("level"))

                      .def("load_mesh", [](Tetrahedralizer &t, const std::string &path) { t.load_mesh(path); }, "loads a mesh", py::arg("path"))
                      .def("set_mesh", [](Tetrahedralizer &t, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) { t.set_mesh(V, F); }, "sets a mesh", py::arg("V"), py::arg("F"))
                      .def("set_meshes", [](Tetrahedralizer &t, const std::vector<Eigen::MatrixXd> &V, const std::vector<Eigen::MatrixXi> &F) { t.set_meshes(V, F); }, "sets several meshes, for boolean", py::arg("V"), py::arg("F"))
                      .def("load_csg_tree", [](Tetrahedralizer &t, const py::object &csg_tree) {  const std::string tmp = py::str(csg_tree); t.boolean_operation(tmp); }, "loads a csg tree, either from file or json", py::arg("csg_tree"))

                      .def("tetrahedralize", [](Tetrahedralizer &t) { t.tetrahedralize(); }, "tetrahedralized the mesh")

                      .def("save", [](Tetrahedralizer &t, const std::string &path, bool smooth_open_boundary, bool manifold_surface, bool correct_surface_orientation) {
                          t.save(path, smooth_open_boundary, manifold_surface, correct_surface_orientation);
                      },
                           "saves the output", py::arg("path"), py::arg("smooth_open_boundary") = false, py::arg("manifold_surface") = false, py::arg("correct_surface_orientation") = false)
                      .def("get_tet_mesh", [](Tetrahedralizer &t, bool smooth_open_boundary, bool manifold_surface, bool correct_surface_orientation) {
                          Eigen::MatrixXd V;
                          Eigen::MatrixXi T;
                          t.get_tet_mesh(smooth_open_boundary, manifold_surface, correct_surface_orientation, V, T);

                          return py::make_tuple(V, T);
                      },
                           "saves the output", py::arg("smooth_open_boundary") = false, py::arg("manifold_surface") = false, py::arg("correct_surface_orientation") = false)
                      .def("get_tet_mesh_from_csg", [](Tetrahedralizer &t, const py::object &csg_tree, bool smooth_open_boundary, bool manifold_surface, bool correct_surface_orientation) {
                          Eigen::MatrixXd V;
                          Eigen::MatrixXi T;
                          const std::string tmp = py::str(csg_tree);
                          t.tree_with_ids = json::parse(tmp);
                          t.has_json_csg = true;

                          t.get_tet_mesh(smooth_open_boundary, manifold_surface, correct_surface_orientation, V, T);

                          t.has_json_csg = false;

                          return py::make_tuple(V, T);
                      },
                           "saves the output", py::arg("csg_tree"), py::arg("smooth_open_boundary") = false, py::arg("manifold_surface") = false, py::arg("correct_surface_orientation") = false)
                      .def("get_stats", [](const Tetrahedralizer &t) { return t.get_stats(); }, "returns the stats");

    tetra.doc() = "Wildmeshing tetrahedralizer";

    m.def("tetrahedralize", [](const std::string &input, const std::string &output, double stop_quality, int max_its, int stage, int stop_p, double epsilon, double edge_length_r, bool mute_log, bool skip_simplify, bool smooth_open_boundary, bool manifold_surface, bool correct_surface_orientation) {
        wildmeshing_binding::init_globals();

        static bool initialized = false;
        if (!initialized)
        {
            Logger::init(!mute_log);
            initialized = true;
        }

        Tetrahedralizer tetra(stop_quality, max_its, stage, stop_p, epsilon, edge_length_r, skip_simplify);
        if (!tetra.load_mesh(input))
            return false;

        tetra.tetrahedralize();
        tetra.save(output, smooth_open_boundary, manifold_surface, correct_surface_orientation);

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
          py::arg("smooth_open_boundary") = false, py::arg("manifold_surface") = false, py::arg("correct_surface_orientation") = false);

    m.def("boolean_operation", [](const py::object &json, const std::string &output, double stop_quality, int max_its, int stage, int stop_p, double epsilon, double edge_length_r, bool mute_log, bool skip_simplify, bool smooth_open_boundary, bool manifold_surface, bool correct_surface_orientation) {
        wildmeshing_binding::init_globals();

        static bool initialized = false;
        if (!initialized)
        {
            Logger::init(!mute_log);
            initialized = true;
        }

        Tetrahedralizer tetra(stop_quality, max_its, stage, stop_p, epsilon, edge_length_r, skip_simplify);

        const std::string tmp = py::str(json);

        if (!tetra.boolean_operation(tmp))
            return false;

        tetra.tetrahedralize();
        tetra.save(output, smooth_open_boundary, manifold_surface, correct_surface_orientation);

        return true;
    },
          "Robust boolean operation, this is an alpha developement version",
          py::arg("json"), // "Input surface mesh INPUT in .off/.obj/.stl/.ply format. (string, required)"
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
          py::arg("smooth_open_boundary") = false, py::arg("manifold_surface") = false, py::arg("correct_surface_orientation") = false);
}
} // namespace wildmeshing_binding