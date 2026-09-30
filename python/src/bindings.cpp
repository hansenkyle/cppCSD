// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

// nanobind bindings for ldcsd_core, built as ldcsd._core (see notes/python-bindings.md).
//
// Arrays cross the boundary in ldcsd.read()'s layout, not the C++ one: corner fields come out as
// [..., nx, 2 (up, down), 2 (L, R)], source as [G, M, nx, 2, 2], bc as [G, 2 (up, down), M].
// Everything returned is a read-only copy, so a Python-side edit can neither bypass the setters'
// validation nor silently vanish, and no view outlives the C++ data it points into.

#include <filesystem>
#include <format>
#include <initializer_list>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "file_manager.h"
#include "input_deck.h"
#include "logger.h"
#include "method.h"
#include "smm.h"
#include "source_iteration.h"
#include "terminal.h"
#include "version.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace {

// const scalar -> nanobind marks the numpy array read-only.
using Array = nb::ndarray<nb::numpy, const double>;
using RowMajor = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

// Hands data to numpy, which frees it when the array is collected. data is in C order for shape.
// Returned already converted, since the getters' default reference_internal policy can't apply to
// an array that owns its data.
nb::object toNumpy(std::vector<double> data, std::initializer_list<size_t> shape) {
  auto* owned = new std::vector<double>(std::move(data));
  nb::capsule owner(owned, [](void* p) noexcept { delete static_cast<std::vector<double>*>(p); });
  return nb::cast(Array(owned->data(), shape, owner));
}

nb::object vector(const Eigen::VectorXd& v) {
  return toNumpy({v.begin(), v.end()}, {static_cast<size_t>(v.size())});
}

nb::object matrix(const Eigen::MatrixXd& m) {
  const RowMajor row_major = m;
  return toNumpy({row_major.data(), row_major.data() + row_major.size()},
                 {static_cast<size_t>(m.rows()), static_cast<size_t>(m.cols())});
}

// [4nx x G] corner field -> [G, nx, 2, 2]. Column-major storage is already that C order: each
// column (group) is 4nx contiguous rows, row 4i + 2*(up, down) + (L, R).
nb::object byGroup(const Eigen::MatrixXd& corners) {
  return toNumpy(
      {corners.data(), corners.data() + corners.size()},
      {static_cast<size_t>(corners.cols()), static_cast<size_t>(corners.rows() / 4), 2, 2});
}

// G matrices [4nx x M] -> [G, M, nx, 2, 2], by the same column-major argument per group.
nb::object byGroupAngle(const std::vector<Eigen::MatrixXd>& per_group) {
  if (per_group.empty()) {
    return toNumpy({}, {0, 0, 0, 2, 2});
  }
  std::vector<double> data;
  for (const Eigen::MatrixXd& group : per_group) {
    data.insert(data.end(), group.data(), group.data() + group.size());
  }
  return toNumpy(std::move(data), {per_group.size(), static_cast<size_t>(per_group[0].cols()),
                                   static_cast<size_t>(per_group[0].rows() / 4), 2, 2});
}

using In3 = nb::ndarray<const double, nb::ndim<3>, nb::c_contig, nb::device::cpu>;
using In5 = nb::ndarray<const double, nb::ndim<5>, nb::c_contig, nb::device::cpu>;

// nanobind's array casters take only float64 numpy arrays; this lets every array argument also be
// a list, an int array, or non-contiguous, as numpy.asarray would accept it. expected names the
// argument and its shape for the error when the rank is wrong.
template <typename T> T asArray(nb::handle obj, const char* expected) {
  const nb::object array =
      nb::module_::import_("numpy").attr("ascontiguousarray")(obj, "dtype"_a = "float64");
  T value;
  if (!nb::try_cast(array, value)) {
    throw std::invalid_argument(
        std::format("expected {}, got shape {}", expected, nb::repr(array.attr("shape")).c_str()));
  }
  return value;
}

// [G, 2, M] -> bc.values [2G x M]: C order is already row 2g + (up, down).
Eigen::MatrixXd bcFromNumpy(const In3& bc) {
  if (bc.shape(1) != 2) {
    throw std::invalid_argument("bc must be shaped [G, 2 (up, down), M]");
  }
  return Eigen::Map<const RowMajor>(bc.data(), 2 * bc.shape(0), bc.shape(2));
}

// [G, M, nx, 2, 2] -> G matrices [4nx x M], inverting byGroupAngle.
std::vector<Eigen::MatrixXd> sourceFromNumpy(const In5& q) {
  if (q.shape(3) != 2 || q.shape(4) != 2) {
    throw std::invalid_argument("source must be shaped [G, M, nx, 2 (up, down), 2 (L, R)]");
  }
  const size_t M = q.shape(1);
  const size_t rows = 4 * q.shape(2);
  std::vector<Eigen::MatrixXd> values(q.shape(0));
  for (size_t g = 0; g < values.size(); ++g) {
    values[g] = Eigen::Map<const Eigen::MatrixXd>(q.data() + g * M * rows, rows, M);
  }
  return values;
}

// Binds what SourceIteration and SecondMoment share beyond Method: construction, write_h5 and
// the solution every method produces.
template <typename Solver> void bindSolver(nb::class_<Solver, Method>& cls) {
  cls.def(nb::init<InputDeck>(), "deck"_a)
      .def(
          "write_h5",
          [](const Solver& self, const std::filesystem::path& path,
             std::optional<std::string> timestamp) {
            self.writeH5(path, timestamp.value_or(make_timestamp()));
          },
          "path"_a, "timestamp"_a = nb::none(),
          "Writes the whole run to a new .h5 file (read it back with ldcsd.read). timestamp "
          "defaults to now.")
      .def_prop_ro(
          "scalar_flux", [](const Solver& self) { return byGroup(self.solution.scalar_flux); },
          "[G, nx, 2, 2]")
      .def_prop_ro(
          "angular_flux",
          [](const Solver& self) { return byGroupAngle(self.solution.angular_flux); },
          "[G, M, nx, 2, 2]");
}

} // namespace

NB_MODULE(_core, m) {
  m.doc() = "ldcsd C++ core: input deck builder and solvers";
  m.attr("__version__") = std::format("{}.{}.{}", version_major, version_minor, version_revision);

  m.def(
      "configure_logging",
      [](const std::filesystem::path& log_path, std::optional<std::filesystem::path> out_path) {
        Logger::configure(log_path);
        if (out_path) {
          Terminal::configure(*out_path);
        }
      },
      "log_path"_a, "out_path"_a = nb::none(),
      "Sends the C++ log to log_path (and terminal output to out_path). Silent until called.");

  // Plain data; checked when handed to InputDeck.set_materials, which also sets its name.
  nb::class_<Material>(m, "Material")
      .def(
          "__init__",
          [](Material* self, nb::handle total, nb::handle S, nb::handle S_b, nb::handle scatter) {
            new (self) Material();
            self->total = asArray<Eigen::VectorXd>(total, "total [G]");
            self->S = asArray<Eigen::VectorXd>(S, "S [G]");
            self->S_b = asArray<Eigen::VectorXd>(S_b, "S_b [G+1]");
            self->scatter = asArray<Eigen::MatrixXd>(scatter, "scatter [G, G]");
          },
          "total"_a, "S"_a, "S_b"_a, "scatter"_a,
          "total [G], S [G] (group-average stopping power), S_b [G+1] (at group boundaries), "
          "scatter [G, G] (from, to), all 0-indexed.")
      .def_ro("name", &Material::name)
      .def_prop_ro("total", [](const Material& self) { return vector(self.total); })
      .def_prop_ro("S", [](const Material& self) { return vector(self.S); })
      .def_prop_ro("S_b", [](const Material& self) { return vector(self.S_b); })
      .def_prop_ro("scatter", [](const Material& self) { return matrix(self.scatter); });

  nb::class_<InputDeck::Mesh>(m, "Mesh")
      .def_ro("n_x", &InputDeck::Mesh::n_x)
      .def_prop_ro("x_boundary",
                   [](const InputDeck::Mesh& self) { return vector(self.x_boundary); });

  nb::class_<InputDeck::Energy>(m, "Energy")
      .def_ro("G", &InputDeck::Energy::G)
      .def_prop_ro("E_boundary",
                   [](const InputDeck::Energy& self) { return vector(self.E_boundary); });

  nb::class_<InputDeck::Angle>(m, "Angle")
      .def_ro("M", &InputDeck::Angle::M)
      .def_prop_ro("mu", [](const InputDeck::Angle& self) { return vector(self.mu); })
      .def_prop_ro("w", [](const InputDeck::Angle& self) { return vector(self.w); });

  nb::class_<InputDeck>(m, "InputDeck")
      .def(nb::init<>())
      .def_static(
          "read",
          [](const std::filesystem::path& path) {
            InputDeck deck;
            deck.load(path);
            return deck;
          },
          "path"_a, "Reads and validates a YAML input deck; raises RuntimeError on any error.")
      .def("to_yaml", &InputDeck::write, "path"_a, "Writes the deck in the format read() reads.")
      .def("validate", &InputDeck::validate,
           "Runs the cross-struct checks (solvers also run it on construction).")
      .def(
          "set_mesh",
          [](InputDeck& self, nb::handle x_boundary) {
            self.set_mesh(asArray<Eigen::VectorXd>(x_boundary, "x_boundary [nx+1]"));
          },
          "x_boundary"_a, "Cell boundaries [nx+1], strictly ascending.")
      .def(
          "set_energy",
          [](InputDeck& self, nb::handle E_boundary) {
            self.set_energy(asArray<Eigen::VectorXd>(E_boundary, "E_boundary [G+1]"));
          },
          "E_boundary"_a, "Group boundaries [G+1], strictly descending.")
      .def(
          "set_angle",
          [](InputDeck& self, nb::handle mu, nb::handle w) {
            self.set_angle(asArray<Eigen::VectorXd>(mu, "mu [M]"),
                           asArray<Eigen::VectorXd>(w, "w [M]"));
          },
          "mu"_a, "w"_a, "Ordinates [M], strictly ascending, and weights [M] summing to 2.")
      .def(
          "set_materials",
          [](InputDeck& self, const nb::dict& materials, const std::vector<std::string>& regions) {
            // Iterated directly (not via std::map) to keep the dict's order.
            std::vector<Material> list;
            for (auto [name, material] : materials) {
              list.push_back(nb::cast<Material>(material));
              list.back().name = nb::cast<std::string>(name);
            }
            self.set_materials(std::move(list), regions);
          },
          "materials"_a, "regions"_a,
          "materials: {name: Material}; regions: one material name per cell.")
      .def(
          "set_bc",
          [](InputDeck& self, nb::handle bc) {
            self.set_bc(bcFromNumpy(asArray<In3>(bc, "bc [G, 2, M]")));
          },
          "bc"_a, "Incoming angular flux [G, 2 (up, down), M].")
      .def(
          "set_source",
          [](InputDeck& self, nb::handle q) {
            self.set_source(sourceFromNumpy(asArray<In5>(q, "source [G, M, nx, 2, 2]")));
          },
          "source"_a, "External source [G, M, nx, 2 (up, down), 2 (L, R)].")
      .def_prop_ro("mesh", [](const InputDeck& self) { return self.mesh; })
      .def_prop_ro("energy", [](const InputDeck& self) { return self.energy; })
      .def_prop_ro("angle", [](const InputDeck& self) { return self.angle; })
      .def_prop_ro("materials",
                   [](const InputDeck& self) {
                     nb::dict materials;
                     for (const Material& material : self.xs.material_list) {
                       materials[material.name.c_str()] = material;
                     }
                     return materials;
                   })
      .def_prop_ro("regions", [](const InputDeck& self) { return self.xs.material_names(); })
      .def_prop_ro(
          "bc",
          [](const InputDeck& self) {
            const RowMajor bc = self.bc.values;
            return toNumpy({bc.data(), bc.data() + bc.size()},
                           {static_cast<size_t>(bc.rows() / 2), 2, static_cast<size_t>(bc.cols())});
          },
          "[G, 2 (up, down), M]")
      .def_prop_ro(
          "source", [](const InputDeck& self) { return byGroupAngle(self.source.values); },
          "[G, M, nx, 2, 2]");

  nb::class_<Method>(m, "Method")
      .def_ro("name", &Method::name)
      .def("solve", &Method::solve, "epsilon"_a, "max_iterations"_a = 1000);

  nb::class_<SourceIteration, Method> source_iteration(m, "SourceIteration");
  bindSolver(source_iteration);

  nb::class_<SecondMoment, Method> second_moment(m, "SecondMoment");
  bindSolver(second_moment);
  second_moment
      .def_prop_ro(
          "current", [](const SecondMoment& self) { return byGroup(self.solution.current); },
          "[G, nx, 2, 2]")
      .def_prop_ro(
          "reconstructed_scalar",
          [](const SecondMoment& self) { return byGroup(self.solution.reconstructed_scalar); },
          "[G, nx, 2, 2]");
}
