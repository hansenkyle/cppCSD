// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "convergence.h"

#include <algorithm>
#include <stdexcept>

namespace {

// Finds the summary for `group`, or nullptr. Group counts are small (tens),
// so a linear scan is cheaper than any index.
GroupSummary* find(std::vector<GroupSummary>& groups, int group) {
  for (GroupSummary& summary : groups) {
    if (summary.group == group) {
      return &summary;
    }
  }
  return nullptr;
}

} // namespace

void ConvergenceHistory::record(IterationRecord record) {
  record.rel_diff = record.phi_norm > 0.0 ? record.abs_diff / record.phi_norm : 0.0;
  records_.push_back(record);

  GroupSummary* summary = find(groups_, record.group);
  if (summary == nullptr) {
    groups_.push_back(GroupSummary{record.group});
    summary = &groups_.back();
  }

  // The roll-up always reflects the group's latest iteration, except
  // max_residual, which is the worst seen anywhere in the group -- a
  // residual that spiked mid-solve is worth surfacing even if the last
  // iteration looked clean.
  summary->iterations = record.iteration;
  summary->phi_norm = record.phi_norm;
  summary->abs_diff = record.abs_diff;
  summary->rel_diff = record.rel_diff;
  summary->max_residual = std::max(summary->max_residual, record.max_residual);
}

void ConvergenceHistory::finishGroup(int group, bool converged, double seconds) {
  GroupSummary* summary = find(groups_, group);
  if (summary == nullptr) {
    throw std::runtime_error("ConvergenceHistory::finishGroup: group " + std::to_string(group) +
                             " has no recorded iterations");
  }
  summary->converged = converged;
  summary->seconds = seconds;
}

int ConvergenceHistory::totalIterations() const { return static_cast<int>(records_.size()); }

double ConvergenceHistory::totalSeconds() const {
  double total = 0.0;
  for (const GroupSummary& summary : groups_) {
    total += summary.seconds;
  }
  return total;
}

bool ConvergenceHistory::allConverged() const {
  return std::all_of(groups_.begin(), groups_.end(),
                     [](const GroupSummary& summary) { return summary.converged; });
}

std::vector<int> ConvergenceHistory::unconvergedGroups() const {
  std::vector<int> groups;
  for (const GroupSummary& summary : groups_) {
    if (!summary.converged) {
      groups.push_back(summary.group);
    }
  }
  return groups;
}

double ConvergenceHistory::worstResidual() const {
  double worst = -1.0;
  for (const IterationRecord& record : records_) {
    worst = std::max(worst, record.max_residual);
  }
  return worst;
}
