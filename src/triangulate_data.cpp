#include "triangulate_data.hpp"

#include "triwild/optimization.h"
#include "triwild/feature.h"
#include "triwild/triangulation.h"
#include "triwild/do_triwild.h"


#include "triwild/meshio.hpp"

#include <Eigen/Dense>


namespace wildmeshing_binding
{

void triangulate_data(py::module &m)
{
    m.def("triangulate_data", [](const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const py::object &feature_info, double stop_quality, int max_its, int stage, double epsilon, double feature_epsilon, double target_edge_len, double edge_length_r, double flat_feature_angle, bool cut_outside, bool skip_eps, const Eigen::MatrixXd &hole_pts, bool mute_log) {
        init_globals();

        json jfeature_info;
        if (!feature_info.is(py::none()))
        {
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
                            epsilon, feature_epsilon,
                            target_edge_len, edge_length_r,
                            flat_feature_angle,
                            cut_outside,
                            hole_pts,
                            mute_log);

        triwild::feature::features.clear();
        triwild::feature::secondary_features.clear();
        return py::make_tuple(V_out, F_out, nodes, F_nodes);
    },
          "Robust Triangulation",
          py::arg("V"),                            // "Input vertices"
          py::arg("E"),                            // "Input edges"
          py::arg("feature_info") = py::none(),    // "Json string containing the features"
          py::arg("stop_quality") = -1,            // "Specify max AMIPS energy for stopping mesh optimization"
          py::arg("max_its") = 80,                 // "Max number of mesh optimization iterations"
          py::arg("stage") = 1,                    // "Specify envelope stage"
          py::arg("epsilon") = -1,                 // "relative envelope epsilon_r. Absolute epsilonn = epsilon_r * diagonal_of_bbox"
          py::arg("feature_epsilon") = 1e-3,       // "Relative feature envelope mu_r. Absolute mu = mu_r * diagonal_of_bbox"
          py::arg("target_edge_len") = -1,         // "Absolute target edge length l"
          py::arg("edge_length_r") = 1. / 20.,     // "Relative target edge length l_r. Absolute l = l_r * diagonal_of_bbox"
          py::arg("flat_feature_angle") = 10.,     // "Desired minimal angle"
          py::arg("cut_outside") = false,          // "Remove \"outside part\""
          py::arg("skip_eps") = false,             // "Skip saving eps"
          py::arg("hole_pts") = Eigen::MatrixXd(), // "Input a #n x 2 matrix of points inside holes you want to remove"
          py::arg("mute_log") = false              // "Mute prints");
    );
}
}