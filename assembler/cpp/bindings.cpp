#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sstream>
#include "GTest.hpp"

namespace py = pybind11;

PYBIND11_MODULE(gtest, m) {
    m.doc() = "Python bindings for Shasta GTest C++ code";

    // Bind the nested Hypothesis class
    py::class_<shasta::GTest::Hypothesis>(m, "Hypothesis")
        .def(py::init<>(), "Default constructor.")
        .def(py::init<const std::vector<std::vector<bool>>&, double>(),
             "Constructor with connectivity matrix and G-value.")
        .def_readwrite("connectivityMatrix", &shasta::GTest::Hypothesis::connectivityMatrix,
             "The connectivity matrix for this hypothesis.")
        .def_readwrite("G", &shasta::GTest::Hypothesis::G,
             "The G-test value for this hypothesis.")
        .def("isForwardInjective", &shasta::GTest::Hypothesis::isForwardInjective,
             "Return true if there is a single exit for each entrance.")
        .def("isBackwardInjective", &shasta::GTest::Hypothesis::isBackwardInjective,
             "Return true if there is a single entrance for each exit.")
        .def("__lt__", &shasta::GTest::Hypothesis::operator<, py::is_operator(),
             "Less-than comparison operator, for sorting by G-value.");

    // Bind the GTest class
    py::class_<shasta::GTest>(m, "GTest")
        // Overloaded constructors
        .def(py::init<const std::vector<std::vector<uint64_t>>&, double>(),
             "Constructor taking a tangle matrix of integers.",
             py::arg("tangleMatrix"), py::arg("epsilon"))
        .def(py::init<const std::vector<std::vector<double>>&, double>(),
             "Constructor taking a tangle matrix of doubles.",
             py::arg("tangleMatrix"), py::arg("epsilon"))

        // Public member variables
        .def_readwrite("success", &shasta::GTest::success,
             "A boolean indicating if the G-test was successful.")
        .def_readwrite("hypotheses", &shasta::GTest::hypotheses,
             "A list of all generated hypotheses.")

        // Public member methods
        .def("writeHtml", [](const shasta::GTest &g) {
            std::stringstream ss;
            g.writeHtml(ss);
            return ss.str();
        }, "Writes the G-test results as an HTML table and returns it as a string.");
}
