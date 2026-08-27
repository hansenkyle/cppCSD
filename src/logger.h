// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef LOGGER_H
#define LOGGER_H

#include <filesystem>
#include <fstream>
#include <string>

enum class LogLevel { Trace, Debug, Info, Warn, Error };

const char* to_string(LogLevel level);

/// @class Logger
/// @brief Minimal, serial logger, exposed as global singleton. Any file that includes @ref logger.h
/// can use the logger after it's configured in main.
/// @details Call Logger::configure() once at startup to set log file location. Calls to log()
/// append a single line to the logfile with timestamp and level. Trace/debug-level calls are
/// optionally compiled, such that release builds skip these with no overhead.
class Logger {
public:
  static Logger& instance() {
    static Logger logger;
    return logger;
  }

  // Opens log_path for writing. Throws std::runtime_error if the file
  // can't be opened.
  static void configure(const std::filesystem::path& log_path);

  void log(LogLevel level, const std::string& message);

private:
  Logger() = default;

  std::ofstream stream_;
};

#ifdef LDCSD_ENABLE_DEBUG_LOGGING
#define LDCSD_LOG_TRACE(message) Logger::instance().log(LogLevel::Trace, message)
#define LDCSD_LOG_DEBUG(message) Logger::instance().log(LogLevel::Debug, message)
#else
#define LDCSD_LOG_TRACE(message) ((void)0)
#define LDCSD_LOG_DEBUG(message) ((void)0)
#endif

#define LDCSD_LOG_INFO(message) Logger::instance().log(LogLevel::Info, message)
#define LDCSD_LOG_WARN(message) Logger::instance().log(LogLevel::Warn, message)
#define LDCSD_LOG_ERROR(message) Logger::instance().log(LogLevel::Error, message)

#endif
