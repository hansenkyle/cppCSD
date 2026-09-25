// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "method.h"
#include "output_block.h"
#include <filesystem>
#include <format>
#include <fstream>
#include <ranges>
#include <stdexcept>
#include <string>

namespace {

std::vector<int> intseq(int stop, int start = 1) {
  std::vector<int> result = {};
  for (int i = start; i < (stop + 1); i++) {
    result.push_back(i);
  }
  return result;
}

std::vector<std::string> int_label_seq(int max) {
  auto intview =
      std::views::iota(1, max + 1) | std::views::transform([](int x) { return std::to_string(x); });
  return std::vector<std::string>(intview.begin(), intview.end());
}

std::vector<std::string> evdoub_to_string(const Eigen::VectorXd& dvec,
                                          std::string_view fmt = "{:.2e}") {
  std::vector<std::string> result;
  for (Eigen::Index i = 0; i < dvec.size(); i++) {
    result.push_back(std::vformat(fmt, std::make_format_args(dvec(i))));
  }
  return result;
}

// Pairs each label with an empty one, for axes whose data has two entries (e.g. two bounds) per
// label.
std::vector<std::string> interleave(const std::vector<std::string>& original) {
  std::vector<std::string> doubled(2 * original.size());
  for (std::size_t i = 0; i < original.size(); i++) {
    doubled[2 * i] = original[i];
  }
  return doubled;
}
} // namespace

Eigen::MatrixXd MethodResult::spectrum() const {
  using Eigen::seqN;
  using Eigen::placeholders::all;
  int I = scalar_flux.rows() / 4;
  int G = scalar_flux.cols();

  // average over each spatial cell
  Eigen::MatrixXd xave =
      (scalar_flux(seqN(0, 2 * I, 2), all) + scalar_flux(seqN(1, 2 * I, 2), all)) / 2;
  // transform (2x by G) to (x by 2g)
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(I, 2 * G);
  for (int i = 0; i < I; i++) {
    for (int g = 0; g < G; g++) {
      result(i, seqN(g * 2, 2)) = xave(seqN(2 * i, 2), g).transpose();
    }
  }
  return result;
}

Eigen::MatrixXd MethodResult::multigroup() const {
  using Eigen::seqN;
  using Eigen::placeholders::all;

  int I = scalar_flux.rows() / 4;
  int G = scalar_flux.cols();

  Eigen::MatrixXd left = (scalar_flux(seqN(0, I, 4), all) + scalar_flux(seqN(2, I, 4), all)) / 2;
  Eigen::MatrixXd right = (scalar_flux(seqN(1, I, 4), all) + scalar_flux(seqN(3, I, 4), all)) / 2;

  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(2 * I, G);
  for (int g = 0; g < G; g++) {
    result(seqN(0, I, 2), g) = left(all, g);
    result(seqN(1, I, 2), g) = right(all, g);
  }

  return result.transpose();
}

Eigen::MatrixXd MethodResult::cell_average_scalar() const {
  using Eigen::seqN;
  using Eigen::placeholders::all;

  int I = scalar_flux.rows() / 4;
  int G = scalar_flux.cols();

  Eigen::MatrixXd leftsum = (scalar_flux(seqN(0, I, 4), all) + scalar_flux(seqN(2, I, 4), all));
  Eigen::MatrixXd rightsum = (scalar_flux(seqN(1, I, 4), all) + scalar_flux(seqN(3, I, 4), all));

  Eigen::MatrixXd result = (leftsum + rightsum) / 4;
  return result;
}

std::vector<Eigen::MatrixXd> MethodResult::cell_average_angular() const {
  using Eigen::seqN;
  using Eigen::placeholders::all;

  int M = angular_flux[0].cols();
  int I = angular_flux[0].rows() / 4;
  int G = angular_flux.size();

  std::vector<Eigen::MatrixXd> result(G);

  for (int g = 0; g < G; g++) {
    Eigen::MatrixXd leftsum =
        (angular_flux[g](seqN(0, I, 4), all) + angular_flux[g](seqN(2, I, 4), all));
    Eigen::MatrixXd rightsum =
        (angular_flux[g](seqN(1, I, 4), all) + angular_flux[g](seqN(3, I, 4), all));

    result[g] = (leftsum + rightsum) / 4;
  }

  return result;
}

void Method::appendToFile(const std::filesystem::path& file_path, const std::string& text) const {
  std::ofstream out(file_path, std::ios::app);
  if (!out.is_open()) {
    throw std::runtime_error("Solver: failed to open '" + file_path.string() + "' for writing");
  }
  out << text;
}

void Method::writeMetadata(const std::filesystem::path& file_path, std::string timestamp) const {
  KeyValueOutput metadata("run info");
  metadata.add("execution date/time", timestamp);
  metadata.add("n_groups", input_deck.energy.G);
  metadata.add("n_cells", input_deck.mesh.n_x);
  metadata.add("n_angles", input_deck.angle.M);
  metadata.add("method", name);
  appendToFile(file_path, metadata.render_txt() + "\n");
}

void Method::writeInputEcho(const std::filesystem::path& file_path) const {

  std::vector<int> x_index = intseq(input_deck.mesh.n_x);
  std::vector<int> e_index = intseq(input_deck.energy.G);
  std::vector<int> m_index = intseq(input_deck.angle.M);

  UnitGroup input_echo("input echo");

  HorizontalTable spatialdata;
  spatialdata.add_row("i", x_index);
  spatialdata.add_row("dx", input_deck.mesh.dx);
  spatialdata.add_row("material", input_deck.xs.material_names());
  spatialdata.add_row("cell bounds", input_deck.mesh.x_boundary);

  HorizontalTable energydata;
  energydata.add_row("g", e_index);
  energydata.add_row("dE", input_deck.energy.dE);
  energydata.add_row("group boundaries", input_deck.energy.E_boundary);

  VerticalTable quadrature("angular quadrature");
  quadrature.add_column("m", m_index);
  quadrature.add_column("mu", input_deck.angle.mu);
  quadrature.add_column("w", input_deck.angle.w);

  VerticalTable boundary("boundary conditions");
  boundary.add_column("m", m_index);
  for (int gplusone : e_index) {
    int g = gplusone - 1;
    boundary.add_column("g" + std::to_string(gplusone) + "_up",
                        input_deck.bc[g](0, Eigen::placeholders::all));
    boundary.add_column("g" + std::to_string(gplusone) + "_down",
                        input_deck.bc[g](1, Eigen::placeholders::all));
  }

  UnitGroup q0("source, zeroth moment", "q integrated over all angles, weight 1");
  for (int gplusone : e_index) {
    int g = gplusone - 1;
    HorizontalTable q0g("g=" + std::to_string(gplusone));
    q0g.add_row("i", x_index);
    q0g.add_row("left,up", input_deck.source.q0[g](Eigen::seqN(0, input_deck.mesh.n_x, 4)));
    q0g.add_row("right,up", input_deck.source.q0[g](Eigen::seqN(1, input_deck.mesh.n_x, 4)));
    q0g.add_row("left,down", input_deck.source.q0[g](Eigen::seqN(2, input_deck.mesh.n_x, 4)));
    q0g.add_row("right,down", input_deck.source.q0[g](Eigen::seqN(3, input_deck.mesh.n_x, 4)));

    q0.add(q0g, "{:.4e}");
  }

  UnitGroup q1("source, first moment", "q integrated over all angles, weight mu");
  for (int gplusone : e_index) {
    int g = gplusone - 1;
    HorizontalTable q1g("g=" + std::to_string(gplusone));
    q1g.add_row("i", x_index);
    q1g.add_row("left,up", input_deck.source.q1[g](Eigen::seqN(0, input_deck.mesh.n_x, 4)));
    q1g.add_row("right,up", input_deck.source.q1[g](Eigen::seqN(1, input_deck.mesh.n_x, 4)));
    q1g.add_row("left,down", input_deck.source.q1[g](Eigen::seqN(2, input_deck.mesh.n_x, 4)));
    q1g.add_row("right,down", input_deck.source.q1[g](Eigen::seqN(3, input_deck.mesh.n_x, 4)));

    q1.add(q1g, "{:.4e}");
  }

  UnitGroup materials("materials");

  for (auto m : input_deck.xs.material_list) {
    // make material's own block
    UnitGroup mat(m.name);
    // add total xs
    HorizontalTable xs_t_and_s("cross sections and stopping power");
    xs_t_and_s.add_row("g", intseq(input_deck.energy.G));
    xs_t_and_s.add_row("sigma_total", m.total);
    xs_t_and_s.add_row("grp. avg. stoping power", m.S);
    xs_t_and_s.add_row("grp. bound. stopping power", m.S_b);
    xs_t_and_s.add_row("in-group scattering xs", m.scatter.diagonal());
    xs_t_and_s.add_row("scattering / total",
                       m.scatter.diagonal().cwiseProduct(m.total.cwiseInverse()));

    HorizontalTable scatter("scattering matrix (from, to)");
    scatter.add_row("", intseq(input_deck.energy.G));
    for (int gplusone : intseq(input_deck.energy.G)) {
      int g = gplusone - 1;
      scatter.add_row(std::to_string(gplusone), m.scatter(g, Eigen::placeholders::all));
    }

    mat.add(xs_t_and_s, "{:.4e}");
    mat.add(scatter, "{:.4e}");

    // add to materials block
    materials.add(mat);
  }

  input_echo.add(spatialdata, "{:.2e}");
  input_echo.add(energydata, "{:.3e}");
  input_echo.add(quadrature, "{:.4e}");
  input_echo.add(boundary, "{:.4e}");
  input_echo.add(q0);
  input_echo.add(q1);
  input_echo.add(materials);

  appendToFile(file_path, input_echo.render_txt());
}

UnitGroup Method::solutionBlock(const MethodResult& solution) const {
  auto iseq = int_label_seq(input_deck.mesh.n_x);
  auto gseq = int_label_seq(input_deck.energy.G);
  auto mseq = int_label_seq(input_deck.angle.M);

  auto x_center = evdoub_to_string(input_deck.mesh.x_center);
  auto mu = evdoub_to_string(input_deck.angle.mu);

  std::vector<std::string> xb, eb;

  for (int i = 0; i < input_deck.mesh.n_x; i++) {
    auto s = evdoub_to_string(input_deck.mesh.x_boundary(Eigen::seqN(i, 2)));
    xb.insert(xb.end(), s.begin(), s.end());
  }

  for (int i = 0; i < input_deck.energy.G; i++) {
    auto s = evdoub_to_string(input_deck.energy.E_boundary(Eigen::seqN(i, 2)));
    eb.insert(eb.end(), s.begin(), s.end());
  }

  // cell-average scalar flux
  MatrixTable scalar("cell-average scalar flux", "averaged over each space-energy cell");
  scalar.set_data(solution.cell_average_scalar().transpose(), x_center, gseq);
  scalar.add_column_label(iseq);
  scalar.set_corner_grid({{"x_i"}, {"g \\ i"}});

  // cell-average angular flux
  UnitGroup angular("cell-average angular flux", "averaged over each space-energy cell");
  const auto cell_average_angular = solution.cell_average_angular();
  for (int g = 0; g < input_deck.energy.G; g++) {
    MatrixTable group("g = " + std::to_string(g + 1));
    group.set_data(cell_average_angular[g].transpose(), x_center, mu);
    group.add_column_label(iseq);
    group.add_row_label(mseq);
    group.set_corner_grid({{"", "x_i"}, {"mu", "m \\i"}});
    angular.add(group, "{:.4e}");
  }

  MatrixTable multigroup("multigroup scalar flux", "averaged over each energy group, not space");
  multigroup.set_data(solution.multigroup(), xb, gseq);
  multigroup.add_column_label(interleave(iseq));
  multigroup.set_corner_grid({{"x_boundary"}, {"g \\ i"}});

  MatrixTable spectrum("energy spectrum", "averaged over each spatial cell, not energy");
  spectrum.set_data(solution.spectrum(), eb, x_center);
  spectrum.add_row_label(iseq);
  spectrum.add_column_label(interleave(gseq));
  spectrum.set_corner_grid({{"", "E_bound"}, {"x_i", "i \\g"}});

  UnitGroup sol_block("solution");
  sol_block.add(scalar, "{:.4e}");
  sol_block.add(angular);
  sol_block.add(multigroup, "{:.4e}");
  sol_block.add(spectrum, "{:.4e}");
  return sol_block;
}

void Method::writeConvergence(const std::filesystem::path& results_path) const {
  auto gseq = int_label_seq(input_deck.energy.G);
  VerticalTable summary("iteration summary");
  summary.add_column("g", gseq);
  summary.add_column("# iterations", convergence_.iterations);

  UnitGroup convergence("per-group convergence history",
                        "delta = (phi_n - phi_n-1). absolute change.");
  for (int g = 0; g < input_deck.energy.G; g++) {
    VerticalTable group("g = " + gseq[g]);
    group.add_column("iteration", int_label_seq(convergence_.records[g].size()));

    std::vector<double> l2, li;
    for (const auto& record : convergence_.records[g]) {
      l2.push_back(record.norm2);
      li.push_back(record.norminf);
    }

    group.add_column("|delta|_2", l2);
    group.add_column("|delta|_infty", li);
    convergence.add(group, "{:.4e}");
  }

  UnitGroup result("convergence");
  result.add(summary);
  result.add(convergence);
  appendToFile(results_path, result.render_txt());
}
