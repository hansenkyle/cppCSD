// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

// nanobind bindings for ldcsd_core, built as ldcsd._core (see notes/python-bindings.md).

#include <format>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include "version.h"

namespace nb = nanobind;

NB_MODULE(_core, m) {
  m.doc() = "ldcsd C++ core: input deck builder and solvers";
  m.attr("__version__") = std::format("{}.{}.{}", version_major, version_minor, version_revision);
}
