// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "logger.h"

#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

const char* to_string(LogLevel level) {
  switch (level) {
  case LogLevel::Trace:
    return "TRACE";
  case LogLevel::Debug:
    return "DEBUG";
  case LogLevel::Info:
    return "INFO";
  case LogLevel::Warn:
    return "WARN";
  case LogLevel::Error:
    return "ERROR";
  }
  return "UNKNOWN";
}

namespace {

std::string timestamp() {
  const std::time_t now = std::time(nullptr);
  const std::tm* tm = std::localtime(&now);
  std::ostringstream oss;
  oss << std::put_time(tm, "%Y-%m-%d %H:%M:%S");
  return oss.str();
}

} // namespace

void Logger::configure(const std::filesystem::path& log_path,
                       const std::filesystem::path& out_path) {
  Logger& logger = instance();
  if (logger.stream_.is_open()) {
    logger.stream_.close();
  }
  logger.stream_.open(log_path);
  if (!logger.stream_.is_open()) {
    throw std::runtime_error("failed to open log file: " + log_path.string());
  }

  if (logger.out_stream_.is_open()) {
    logger.out_stream_.close();
  }
  logger.out_stream_.open(out_path);
  if (!logger.out_stream_.is_open()) {
    throw std::runtime_error("failed to open out file: " + out_path.string());
  }
}

void Logger::log(LogLevel level, const std::string& message, bool echo) {
  const std::string line = "[" + timestamp() + "] [" + to_string(level) + "] " + message;
  stream_ << line << "\n";
  stream_.flush();
  if (echo) {
    std::cout << line << "\n";
  }
}

void Logger::print(const std::string& message) {
  std::cout << message << "\n";
  out_stream_ << message << "\n";
  out_stream_.flush();
}
