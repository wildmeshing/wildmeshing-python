#include "tetrahedralize.hpp"

#include <floattetwild/MeshIO.hpp>
#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/Statistics.h>
#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/CSGTreeParser.hpp>
#include <floattetwild/MshLoader.h>

#include <floattetwild/Logger.hpp>

#include <igl/Timer.h>
#include <igl/writeSTL.h>
#include <igl/remove_unreferenced.h>

#include <memory>

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <thread>
#endif

using namespace floatTetWild;
using namespace Eigen;

namespace wildmeshing_binding
{
    Tetrahedralizer::Tetrahedralizer(
        double stop_quality, int max_its, int stage, int stop_p,
        double epsilon, double edge_length_r,
        bool skip_simplify, bool coarsen) : skip_simplify(skip_simplify)
    {
        wildmeshing_binding::init_globals();

        Parameters &params = mesh.params;

        // params.input_path = input;
        // params.output_path = output;

        params.stop_energy = stop_quality;
        params.max_its = max_its;
        params.stage = stage;
        params.stop_p = stop_p;

        params.coarsen = coarsen;

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

    void Tetrahedralizer::set_num_threads(int num_threads)
    {
        Parameters &params = mesh.params;
        params.num_threads = num_threads;
    }

    void Tetrahedralizer::set_log_level(int level)
    {
        Parameters &params = mesh.params;
        params.log_level = std::max(0, std::min(6, level));
        spdlog::set_level(static_cast<spdlog::level::level_enum>(params.log_level));
        spdlog::flush_every(std::chrono::seconds(3));
    }

    bool Tetrahedralizer::load_mesh(const std::string &path, const std::string &tag_path, const std::vector<double> &epsr_tags)
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

        params.input_epsr_tags = epsr_tags;

#ifdef NEW_ENVELOPE
        if (!MeshIO::load_mesh(params.input_path, input_vertices, input_faces, sf_mesh, input_tags, params.input_epsr_tags))
#else
        if (!MeshIO::load_mesh(params.input_path, input_vertices, input_faces, sf_mesh, input_tags))
#endif
        {
            throw std::invalid_argument("Invalid mesh path at " + params.input_path);
            return false;
        }
        else if (input_vertices.empty() || input_faces.empty())
        {
            throw std::invalid_argument("Invalid mesh path at " + params.input_path);
            return false;
        }

        if (!params.input_epsr_tags.empty())
        {
            if (params.input_epsr_tags.size() != input_vertices.size())
            {
                throw std::invalid_argument("epsr_tags need to be same size as vertices, " + std::to_string(params.input_epsr_tags.size()) + " vs " + std::to_string(input_vertices.size()));
                return false;
            }
        }

        return load_mesh_aux();
    }

    bool Tetrahedralizer::boolean_operation(const std::string &json_string)
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
        if (!ok)
        {
            throw std::invalid_argument("Invalid mesh path in the json");
        }

        load_mesh_aux();
        return ok;
    }

    void Tetrahedralizer::set_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const std::vector<double> &epsr_tags)
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

        Parameters &params = mesh.params;
        params.input_epsr_tags = epsr_tags;

        if (!params.input_epsr_tags.empty())
        {
            if (params.input_epsr_tags.size() != input_vertices.size())
            {
                throw std::invalid_argument("epsr_tags need to be same size as vertices, " + std::to_string(params.input_epsr_tags.size()) + " vs " + std::to_string(input_vertices.size()));
            }
        }

#ifdef NEW_ENVELOPE
        MeshIO::load_mesh(input_vertices, input_faces, sf_mesh, input_tags, params.input_epsr_tags);
#else
        MeshIO::load_mesh(input_vertices, input_faces, sf_mesh, input_tags);
#endif

        load_mesh_aux();
    }

    void Tetrahedralizer::set_meshes(const std::vector<Eigen::MatrixXd> &V, const std::vector<Eigen::MatrixXi> &F)
    {
        std::vector<std::vector<Vector3>> vs(V.size());
        std::vector<std::vector<floatTetWild::Vector3i>> fs(F.size());

        if (V.size() != F.size())
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

        for (int j = 0; j < F.size(); ++j)
        {
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

    void Tetrahedralizer::set_sizing_field(const std::string &path)
    {
        PyMesh::MshLoader mshLoader(path);
        Eigen::VectorXd V_in = mshLoader.get_nodes();
        Eigen::VectorXi T_in = mshLoader.get_elements();
        Eigen::VectorXd values = mshLoader.get_node_field("values");

        set_sizing_field(V_in, T_in, values);
    }

    void Tetrahedralizer::set_sizing_field(const Eigen::MatrixXd &V, Eigen::MatrixXi &T, const Eigen::VectorXd &values)
    {
        assert(V.cols() == 3);
        assert(T.cols() == 4);
        assert(values.size() == V.rows());

        if (V.cols() != 3)
            throw std::invalid_argument("V should have 3 cols");
        if (T.cols() != 4)
            throw std::invalid_argument("T should have 3 cols");
        if (values.size() != V.rows())
            throw std::invalid_argument("values should have the same length as V.rows()");

        Eigen::VectorXd V_in(V.size());
        Eigen::VectorXi T_in(T.size());

        int index = 0;
        for (int i = 0; i < V.rows(); ++i)
        {
            for (int j = 0; j < 3; ++j)
                V_in[index++] = V(i, j);
        }

        index = 0;
        for (int i = 0; i < T.rows(); ++i)
        {
            for (int j = 0; j < 4; ++j)
                T_in[index++] = T(i, j);
        }

        set_sizing_field(V_in, T_in, values);
    }

    void Tetrahedralizer::set_sizing_field(const Eigen::VectorXd &V_in, const Eigen::VectorXi &T_in, const Eigen::VectorXd &values)
    {
        Parameters &params = mesh.params;
        params.apply_sizing_field = true;
        params.get_sizing_field_value = [V_in, T_in, values](const Vector3 &p) {
            GEO::Mesh bg_mesh;
            bg_mesh.vertices.clear();
            bg_mesh.vertices.create_vertices((int)V_in.rows() / 3);
            for (int i = 0; i < V_in.rows() / 3; i++)
            {
                GEO::vec3 &p = bg_mesh.vertices.point(i);
                for (int j = 0; j < 3; j++)
                    p[j] = V_in(i * 3 + j);
            }
            bg_mesh.cells.clear();
            bg_mesh.cells.create_tets((int)T_in.rows() / 4);
            for (int i = 0; i < T_in.rows() / 4; i++)
            {
                for (int j = 0; j < 4; j++)
                    bg_mesh.cells.set_vertex(i, j, T_in(i * 4 + j));
            }

            GEO::MeshCellsAABB bg_aabb(bg_mesh, false);
            GEO::vec3 geo_p(p[0], p[1], p[2]);
            int bg_t_id = bg_aabb.containing_tet(geo_p);
            if (bg_t_id == GEO::MeshCellsAABB::NO_TET)
                return -1.;

            //compute barycenter
            std::array<Vector3, 4> vs;
            for (int j = 0; j < 4; j++)
            {
                vs[j] = Vector3(V_in(T_in(bg_t_id * 4 + j) * 3), V_in(T_in(bg_t_id * 4 + j) * 3 + 1),
                                V_in(T_in(bg_t_id * 4 + j) * 3 + 2));
            }
            double value = 0;
            for (int j = 0; j < 4; j++)
            {
                Vector3 n = ((vs[(j + 1) % 4] - vs[j]).cross(vs[(j + 2) % 4] - vs[j])).normalized();
                double d = (vs[(j + 3) % 4] - vs[j]).dot(n);
                if (d == 0)
                    continue;
                double weight = abs((p - vs[j]).dot(n) / d);
                value += weight * values(T_in(bg_t_id * 4 + (j + 3) % 4));
            }
            return value; // / mesh.params.ideal_edge_length;
        };
    }

    void Tetrahedralizer::set_sizing_field(std::function<double(const Vector3 &p)> &field)
    {
        Parameters &params = mesh.params;
        params.apply_sizing_field = true;
        params.get_sizing_field_value = field;
    }

    bool Tetrahedralizer::load_mesh_aux()
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

#ifdef NEW_ENVELOPE
        if (!params.input_epsr_tags.empty())
            tree->init_sf_tree(input_vertices, input_faces, params.input_epsr_tags, params.bbox_diag_length);
        else
            tree->init_sf_tree(input_vertices, input_faces, params.eps);
#endif

        stats().record(StateInfo::init_id, 0, input_vertices.size(), input_faces.size(), -1, -1);
        return true;
    }

    void Tetrahedralizer::tetrahedralize()
    {
        Parameters &params = mesh.params;
        igl::Timer timer;

        timer.start();
        simplify(input_vertices, input_faces, input_tags, *tree, params, skip_simplify);
        tree->init_b_mesh_and_tree(input_vertices, input_faces, mesh);
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

    void Tetrahedralizer::save(const std::string &path, bool smooth_open_boundary, bool floodfill, bool manifold_surface, bool use_input_for_wn, bool correct_surface_orientation, bool all_mesh, bool binary, int boolean_op)
    {
        igl::Timer timer;

        Mesh mesh_copy = mesh;
        Eigen::MatrixXd Vout;
        Eigen::MatrixXi Fout;

        Parameters &params = mesh_copy.params;
        params.output_path = path;
        params.correct_surface_orientation = correct_surface_orientation;
        params.smooth_open_boundary = smooth_open_boundary;
        params.manifold_surface = manifold_surface;
        params.use_input_for_wn = use_input_for_wn;

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

        if (!all_mesh)
        {
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
                {
                    if (floodfill)
                        filter_outside_floodfill(mesh_copy);
                    else
                        filter_outside(mesh_copy);
                }
            }
            if (params.manifold_surface)
            {
                floatTetWild::manifold_surface(mesh_copy, Vout, Fout);
            }
        }

        stats().record(StateInfo::wn_id, timer.getElapsedTimeInSec(), mesh_copy.get_v_num(), mesh_copy.get_t_num(),
                       mesh_copy.get_max_energy(), mesh_copy.get_avg_energy());
        logger().info("after winding number");
        logger().info("#v = {}", mesh_copy.get_v_num());
        logger().info("#t = {}", mesh_copy.get_t_num());
        logger().info("winding number {}s", timer.getElapsedTimeInSec());
        logger().info("");

        if (params.output_path.size() > 3 && params.output_path.substr(params.output_path.size() - 3, params.output_path.size()) == "msh")
            MeshIO::write_mesh(params.output_path, mesh_copy, false, std::vector<double>(), binary);
        else if (params.output_path.size() > 4 && params.output_path.substr(params.output_path.size() - 4, params.output_path.size()) == "mesh")
            MeshIO::write_mesh(params.output_path, mesh_copy, false, std::vector<double>(), binary);
        else
            MeshIO::write_mesh(params.output_path + "_" + params.postfix + ".msh", mesh_copy, false, std::vector<double>(), binary);
    }

    void Tetrahedralizer::get_tet_mesh(bool smooth_open_boundary, bool floodfill, bool manifold_surface, bool use_input_for_wn, bool correct_surface_orientation, bool all_mesh, Eigen::MatrixXd &Vs, Eigen::MatrixXi &Ts, Eigen::MatrixXd &flags, int boolean_op)
    {
        igl::Timer timer;

        Mesh mesh_copy = mesh;
        Eigen::MatrixXd V;
        Eigen::MatrixXi T;

        Eigen::MatrixXd Vout;
        Eigen::MatrixXi Fout;

        const auto skip_tet = [&mesh_copy](const int i) { return mesh_copy.tets[i].is_removed; };
        const auto skip_vertex = [&mesh_copy](const int i) { return mesh_copy.tet_vertices[i].is_removed; };
        std::vector<int> t_ids(mesh_copy.tets.size());
        std::iota(std::begin(t_ids), std::end(t_ids), 0);

        Parameters &params = mesh_copy.params;
        params.correct_surface_orientation = correct_surface_orientation;
        params.smooth_open_boundary = smooth_open_boundary;
        params.manifold_surface = manifold_surface;
        params.use_input_for_wn = use_input_for_wn;

        if (!all_mesh)
        {
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
                {
                    if (floodfill)
                        filter_outside_floodfill(mesh_copy);
                    else
                        filter_outside(mesh_copy);
                }
            }

            if (params.manifold_surface)
            {
                floatTetWild::manifold_surface(mesh_copy, Vout, Fout);
            }

            stats().record(StateInfo::wn_id, timer.getElapsedTimeInSec(), mesh_copy.get_v_num(), mesh_copy.get_t_num(),
                           mesh_copy.get_max_energy(), mesh_copy.get_avg_energy());
        }
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
        flags.resize(cnt_t, 1);
        index = 0;

        const std::array<int, 4> new_indices = {{0, 1, 3, 2}};

        for (const int i : t_ids)
        {
            if (skip_tet(i))
                continue;
            for (int j = 0; j < 4; j++)
            {
                T(index, j) = old_2_new[mesh_copy.tets[i][new_indices[j]]];
            }
            flags(index) = mesh_copy.tets[i].scalar;
            index++;
        }

        Eigen::MatrixXi I;
        igl::remove_unreferenced(V, T, Vs, Ts, I);
    }

    std::string Tetrahedralizer::get_stats() const
    {
        std::stringstream ss;
        ss << stats();
        return ss.str();
    }

#ifndef WILDMESHING_SKIP_BINDINGS
    void tetrahedralize(py::module &m)
    {
        auto &tetra = py::class_<Tetrahedralizer>(m, "Tetrahedralizer")
                          .def(py::init<
                                   double, int, int, int,
                                   double, double,
                                   bool, bool>(),
                               py::arg("stop_quality") = 10,        // "Specify max AMIPS energy for stopping mesh optimization"
                               py::arg("max_its") = 80,             // "Max number of mesh optimization iterations"
                               py::arg("stage") = 2,                // "Specify envelope stage"
                               py::arg("stop_p") = -1,              //
                               py::arg("epsilon") = 1e-3,           // "relative envelope epsilon_r. Absolute epsilonn = epsilon_r * diagonal_of_bbox"
                               py::arg("edge_length_r") = 1. / 20., // "Relative target edge length l_r. Absolute l = l_r * diagonal_of_bbox"
                               py::arg("skip_simplify") = false,    //
                               py::arg("coarsen") = true)

                          .def(
                              "set_log_level", [](Tetrahedralizer &t, int level) { t.set_log_level(level); }, "sets log level, valid value between 0 (all logs) and 6 (no logs)", py::arg("level"))

                          .def(
                              "load_mesh", [](Tetrahedralizer &t, const std::string &path, const std::vector<double> &epsr_tags) { t.load_mesh(path, "", epsr_tags); }, "loads a mesh", py::arg("path"), py::arg("epsr_tags") = std::vector<double>())
                          .def(
                              "load_sizing_field", [](Tetrahedralizer &t, const std::string &path) { t.set_sizing_field(path); }, "load sizing field", py::arg("path"))
                          .def(
                              "set_sizing_field", [](Tetrahedralizer &t, const Eigen::MatrixXd &V, Eigen::MatrixXi &T, const Eigen::VectorXd &values) { t.set_sizing_field(V, T, values); }, "set sizing field", py::arg("V"), py::arg("T"), py::arg("values"))
                          .def(
                              "set_sizing_field_from_func", [](Tetrahedralizer &t, std::function<double(const Vector3 &p)> &field) { t.set_sizing_field(field); }, "set sizing field", py::arg("field"))
                          .def(
                              "set_mesh", [](Tetrahedralizer &t, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const std::vector<double> &epsr_tags) { t.set_mesh(V, F, epsr_tags); }, "sets a mesh", py::arg("V"), py::arg("F"), py::arg("epsr_tags") = std::vector<double>())
                          .def(
                              "set_meshes", [](Tetrahedralizer &t, const std::vector<Eigen::MatrixXd> &V, const std::vector<Eigen::MatrixXi> &F) { t.set_meshes(V, F); }, "sets several meshes, for boolean", py::arg("V"), py::arg("F"))
                          .def(
                              "load_csg_tree", [](Tetrahedralizer &t, const py::object &csg_tree) {  const std::string tmp = py::str(csg_tree); t.boolean_operation(tmp); }, "loads a csg tree, either from file or json", py::arg("csg_tree"))

                          .def(
                              "tetrahedralize", [](Tetrahedralizer &t) { t.tetrahedralize(); }, "tetrahedralized the mesh")

                          .def(
                              "save", [](Tetrahedralizer &t, const std::string &path, bool smooth_open_boundary, bool floodfill, bool use_input_for_wn, bool manifold_surface, bool correct_surface_orientation, bool all_mesh, bool binary) {
                                  t.save(path, smooth_open_boundary, floodfill, use_input_for_wn, manifold_surface, correct_surface_orientation, all_mesh, binary);
                              },
                              "saves the output", py::arg("path"), py::arg("smooth_open_boundary") = false, py::arg("floodfill") = false, py::arg("use_input_for_wn") = false, py::arg("manifold_surface") = false, py::arg("correct_surface_orientation") = false, py::arg("all_mesh") = false, py::arg("binary") = true)
                          .def(
                              "get_tet_mesh", [](Tetrahedralizer &t, bool smooth_open_boundary, bool floodfill, bool use_input_for_wn, bool manifold_surface, bool correct_surface_orientation, bool all_mesh) {
                                  Eigen::MatrixXd V;
                                  Eigen::MatrixXi T;
                                  Eigen::MatrixXd tags;
                                  t.get_tet_mesh(smooth_open_boundary, floodfill, manifold_surface, use_input_for_wn, correct_surface_orientation, all_mesh, V, T, tags);

                                  return py::make_tuple(V, T, tags);
                              },
                              "gets the output", py::arg("smooth_open_boundary") = false, py::arg("floodfill") = false, py::arg("use_input_for_wn") = false, py::arg("manifold_surface") = false, py::arg("correct_surface_orientation") = false, py::arg("all_mesh") = false)
                          .def(
                              "get_tet_mesh_from_csg", [](Tetrahedralizer &t, const py::object &csg_tree, bool manifold_surface, bool use_input_for_wn, bool correct_surface_orientation) {
                                  Eigen::MatrixXd V;
                                  Eigen::MatrixXi T;
                                  Eigen::MatrixXd tags;
                                  const std::string tmp = py::str(csg_tree);
                                  t.tree_with_ids = json::parse(tmp);
                                  t.has_json_csg = true;

                                  t.get_tet_mesh(false, false, manifold_surface, use_input_for_wn, correct_surface_orientation, false, V, T, tags);

                                  t.has_json_csg = false;

                                  return py::make_tuple(V, T, tags);
                              },
                              "gets the output from a csg tree", py::arg("csg_tree"), py::arg("manifold_surface") = false, py::arg("use_input_for_wn") = false, py::arg("correct_surface_orientation") = false)
                          .def(
                              "get_stats", [](const Tetrahedralizer &t) { return t.get_stats(); }, "returns the stats");

        tetra.doc() = "Wildmeshing tetrahedralizer";

        m.def(
            "tetrahedralize", [](const std::string &input, const std::string &output, double stop_quality, int max_its, int stage, int stop_p, double epsilon, double edge_length_r, bool mute_log, bool skip_simplify, bool coarsen, bool smooth_open_boundary, bool floodfill, bool use_input_for_wn, bool manifold_surface, bool correct_surface_orientation, bool all_mesh, bool binary) {
                wildmeshing_binding::init_globals();

                static bool initialized = false;
                if (!initialized)
                {
                    Logger::init(!mute_log);
                    initialized = true;
                }

                Tetrahedralizer tetra(stop_quality, max_its, stage, stop_p, epsilon, edge_length_r, skip_simplify, coarsen);
                if (!tetra.load_mesh(input, "", std::vector<double>()))
                    return false;

                tetra.tetrahedralize();
                tetra.save(output, smooth_open_boundary, floodfill, use_input_for_wn, manifold_surface, correct_surface_orientation, all_mesh, binary);

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
            py::arg("coarsen") = true, py::arg("smooth_open_boundary") = false, py::arg("floodfill") = false, py::arg("manifold_surface") = false, py::arg("use_input_for_wn") = false, py::arg("correct_surface_orientation") = false, py::arg("all_mesh") = false, py::arg("binary") = true);

        m.def(
            "boolean_operation", [](const py::object &json, const std::string &output, double stop_quality, int max_its, int stage, int stop_p, double epsilon, double edge_length_r, bool mute_log, bool skip_simplify, bool coarsen, bool manifold_surface, bool use_input_for_wn, bool correct_surface_orientation, bool all_mesh, bool binary) {
                wildmeshing_binding::init_globals();

                static bool initialized = false;
                if (!initialized)
                {
                    Logger::init(!mute_log);
                    initialized = true;
                }

                Tetrahedralizer tetra(stop_quality, max_its, stage, stop_p, epsilon, edge_length_r, skip_simplify, coarsen);

                const std::string tmp = py::str(json);

                if (!tetra.boolean_operation(tmp))
                    return false;

                tetra.tetrahedralize();
                tetra.save(output, false, false, manifold_surface, use_input_for_wn, correct_surface_orientation, all_mesh, binary);

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
            py::arg("coarsen") = true, py::arg("manifold_surface") = false, py::arg("use_input_for_wn") = false, py::arg("correct_surface_orientation") = false, py::arg("all_mesh") = false, py::arg("binary") = true);
    }
#endif
} // namespace wildmeshing_binding