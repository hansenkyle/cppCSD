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

struct IterationRecord {
public:
  IterationRecord();
  IterationRecord(double l2, double li);
  double norm2;
  double norminf;
};

class ConvergenceHistory {
public:
  ConvergenceHistory(int g);
  void log_group(int g, IterationRecord record);
  void time_group(int g, double time);
  std::vector<int> iterations;
  std::vector<double> group_times;
  std::vector<std::vector<IterationRecord>> records;
};

#endif
