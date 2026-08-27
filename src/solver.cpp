#include "solver.h"

Solver::Solver(const Mesh& mesh, const CrossSection& cross_section)
    : mesh(mesh), cross_section(this->mesh, cross_section.total, cross_section.scattering,
                                cross_section.stop_power, cross_section.stop_power_boundary,
                                cross_section.material) {}

void Solver::sweep() {}
