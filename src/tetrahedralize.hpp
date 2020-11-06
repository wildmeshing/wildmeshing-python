#pragma once

#include "Utils.hpp"

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>
#include <floattetwild/Types.hpp>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>

#include <Eigen/Dense>

namespace wildmeshing_binding
{
    class Tetrahedralizer
    {
    public:
        floatTetWild::Mesh mesh;

        std::vector<floatTetWild::Vector3> input_vertices;
        std::vector<floatTetWild::Vector3i> input_faces;
        std::vector<int> input_tags;
        GEO::Mesh sf_mesh;
        std::unique_ptr<floatTetWild::AABBWrapper> tree;
        bool skip_simplify;
        floatTetWild::json tree_with_ids;
        bool has_json_csg = false;

        Tetrahedralizer(
            double stop_quality, int max_its, int stage, int stop_p,
            double epsilon, double edge_length_r,
            bool skip_simplify, bool coarsen);

    private:
        void set_num_threads(int num_threads);

    public:
        void set_log_level(int level);

        bool boolean_operation(const std::string &json_string);
        void set_meshes(const std::vector<Eigen::MatrixXd> &V, const std::vector<Eigen::MatrixXi> &F);

        bool load_mesh(const std::string &path, const std::string &tag_path, const std::vector<double> &epsr_tags);
        void set_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const std::vector<double> &epsr_tags);

        void set_sizing_field(const std::string &path);
        void set_sizing_field(const Eigen::MatrixXd &V, Eigen::MatrixXi &T, const Eigen::VectorXd &values);
        void set_sizing_field(std::function<double(const floatTetWild::Vector3 &p)> &field);

    private:
        void set_sizing_field(const Eigen::VectorXd &V_in, const Eigen::VectorXi &T_in, const Eigen::VectorXd &values);
        bool load_mesh_aux();

    public:
        void tetrahedralize();

        void save(const std::string &path, bool smooth_open_boundary, bool floodfill, bool manifold_surface, bool use_input_for_wn, bool correct_surface_orientation, bool all_mesh, bool binary, int boolean_op = -1);
        void get_tet_mesh(bool smooth_open_boundary, bool floodfill, bool manifold_surface, bool use_input_for_wn, bool correct_surface_orientation, bool all_mesh, Eigen::MatrixXd &Vs, Eigen::MatrixXi &Ts, Eigen::MatrixXd &flags, int boolean_op = -1);
        std::string get_stats() const;
    };
#ifndef WILDMESHING_SKIP_BINDINGS
    void tetrahedralize(py::module &m);
#endif
} // namespace wildmeshing_binding
