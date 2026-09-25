// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "convergence.h"

#include <vector>

IterationRecord::IterationRecord(Norms solution, Norms delta) : solution(solution), delta(delta) {}

ConvergenceHistory::ConvergenceHistory(int G) {
  iterations = std::vector<int>(G);
  group_times = std::vector<double>(G);
  records = std::vector<std::vector<IterationRecord>>(G);
}

void ConvergenceHistory::log_group(int g, IterationRecord record) {
  iterations[g]++;
  records[g].push_back(record);
}

void ConvergenceHistory::time_group(int g, double time) { group_times[g] = time; }