// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef CONVERGENCE_H
#define CONVERGENCE_H

#include <string>
#include <vector>

// One iteration's worth of convergence data, recorded by whichever
// iterative method produced it. Plain data: no method owns it, and nothing
// here knows what source iteration is -- a future method (SMM, Krylov,
// whatever) records into the same type.
//
// Norms are of the group's scalar flux vector, length 4 * n_x.
struct IterationRecord {
  int group = 0;             // energy group being solved
  int iteration = 0;         // 1-based, counted within the group
  double phi_norm = 0.0;     // ||phi_new||_2
  double abs_diff = 0.0;     // ||phi_new - phi_old||_2
  double rel_diff = 0.0;     // abs_diff / ||phi_new||_2 (0 if phi_norm is 0)
  double max_residual = 0.0; // max |Ax - b| over the group, or -1 if not evaluated
};

// What happened in one energy group, derived from its IterationRecords.
struct GroupSummary {
  int group = 0;
  int iterations = 0;
  double phi_norm = 0.0;
  double abs_diff = 0.0;
  double rel_diff = 0.0;
  double max_residual = 0.0;
  double seconds = 0.0;
  bool converged = false;
};

// The convergence record for a whole solve: every iteration of every group,
// in the order they happened, plus the per-group roll-up derived from them.
//
// Deliberately storage-only. It doesn't decide when a group has converged
// (the method does that and reports it via finishGroup) and it doesn't
// format anything (SolverFormatter does that) -- so the same history can be
// written to a file, printed to the terminal, or inspected from a test
// without any of those growing knowledge of the others.
class ConvergenceHistory {
public:
  // Appends one iteration. Fills in rel_diff from abs_diff/phi_norm.
  void record(IterationRecord record);

  // Closes out the group most recently recorded into: stamps whether it
  // converged and how long it took. Call once per group, after its last
  // record().
  void finishGroup(int group, bool converged, double seconds);

  const std::vector<IterationRecord>& records() const { return records_; }

  // One entry per group that was recorded into, in group order.
  const std::vector<GroupSummary>& groups() const { return groups_; }

  // Total iterations across all groups.
  int totalIterations() const;

  // Total wall time across all finished groups, in seconds.
  double totalSeconds() const;

  // True only if every finished group converged. An empty history is
  // trivially converged.
  bool allConverged() const;

  // Groups that hit the iteration cap, in group order. Empty when
  // allConverged().
  std::vector<int> unconvergedGroups() const;

  // Largest max_residual over every recorded iteration, ignoring
  // iterations where residuals weren't evaluated. Returns -1 if none were.
  double worstResidual() const;

  bool empty() const { return records_.empty(); }

private:
  std::vector<IterationRecord> records_;
  std::vector<GroupSummary> groups_;
};

#endif
