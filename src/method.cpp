// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "method.h"
#include "output_block.h"
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {
std::string timestamp() {
  const std::time_t now = std::time(nullptr);
  const std::tm* tm = std::localtime(&now);
  std::ostringstream oss;
  oss << std::put_time(tm, "%Y-%m-%d %H:%M:%S");
  return oss.str();
}

std::vector<int> intseq(int stop, int start = 1) {
  std::vector<int> result = {};
  for (int i = start; i < (stop + 1); i++) {
    result.push_back(i);
  }
  return result;
}
} // namespace

void Method::appendToFile(const std::filesystem::path& file_path, const std::string& text) const {
  std::ofstream out(file_path, std::ios::app);
  if (!out.is_open()) {
    throw std::runtime_error("Solver: failed to open '" + file_path.string() + "' for writing");
  }
  out << text;
}

void Method::writeMetadata(const std::filesystem::path& file_path) const {
  KeyValueOutput metadata("run info");
  metadata.add("execution date/time", timestamp());
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