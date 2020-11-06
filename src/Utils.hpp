#pragma once

#ifndef WILDMESHING_SKIP_BINDINGS
#include <pybind11/pybind11.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>

namespace py = pybind11;
#endif
namespace wildmeshing_binding
{
    void init_globals();
}
