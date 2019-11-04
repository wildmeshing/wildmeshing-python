#include "Utils.hpp"

#include "triangulate_data.hpp"
#include "triangulate.hpp"
#include "tetrahedralize.hpp"

PYBIND11_MODULE(wildmeshing, m)
{
    // m.doc() = "..."
    wildmeshing_binding::triangulate_data(m);
    wildmeshing_binding::triangulate(m);
    wildmeshing_binding::tetrahedralize(m);
}