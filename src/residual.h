#ifndef RESIDUAL_H
#define RESIDUAL_H

#include "cross_section.h"
#include "fe_space.h"
#include "mesh.h"

// Evaluates the discretized transport equation's residual for a candidate
// solution. Deliberately re-derives the governing equations independently of
// whatever matrix a Solver assembled, rather than reusing it -- the point is
// to catch bugs in that assembly, which a residual computed from the same
// matrix could never reveal. Takes plain problem data (Mesh, CrossSection,
// FESpace) and a candidate solution rather than a Solver, so it has no
// access to and no relationship with any Solver's internals. Templated on
// Scalar so the same formula can be evaluated in `double` for routine
// diagnostics or in an extended-precision type (e.g.
// boost::multiprecision::cpp_bin_float_50) to verify convergence to machine
// precision without the evaluation's own rounding error swamping the
// result.
//
// TODO: takes no solution parameter yet -- the corner-quad solution storage
// this will read from doesn't exist until something actually produces one
// (Solver::sweep isn't implemented). Add it then, rather than now.
template <typename Scalar>
Scalar residualNorm(const Mesh& mesh, const CrossSection& cross_section, const FESpace& fe_space) {
  return Scalar{};
}

#endif
