// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef LOGGER_H
#define LOGGER_H

#include <filesystem>
#include <format>
#include <fstream>
#include <string>

enum class Level { Trace, Debug, Info, Warn, Error };
enum class Channel { General, Iteration };

const char* to_string(Level level);

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

  // Opens log_path  for writing. Throws std::runtime_error if
  // it can't be opened.
  static void configure(const std::filesystem::path& log_path);

  // Appends a level-tagged, timestamped line to the log file. The
  // channel-less overload logs to Channel::General.
  void log(Level level, const std::string& message);
  void log(Channel channel, Level level, const std::string& message);

private:
  Logger() = default;

  std::ofstream stream_;
};

// Picks NAME based on how many variadic arguments were actually passed to
// the calling macro: 1 argument (just a message) selects the DEFAULT slot,
// 2 arguments (channel, message) selects the CH slot. Shared by every
// LDCSD_LOG_* macro below so each level only has to name its two shapes.
#define LDCSD_LOG_PICK_ARITY(_1, _2, NAME, ...) NAME

#ifdef LDCSD_ENABLE_DEBUG_LOGGING
#define LDCSD_LOG_TRACE(...)                                                                       \
  LDCSD_LOG_PICK_ARITY(__VA_ARGS__, LDCSD_LOG_TRACE_CH, LDCSD_LOG_TRACE_DEFAULT)(__VA_ARGS__)
#define LDCSD_LOG_DEBUG(...)                                                                       \
  LDCSD_LOG_PICK_ARITY(__VA_ARGS__, LDCSD_LOG_DEBUG_CH, LDCSD_LOG_DEBUG_DEFAULT)(__VA_ARGS__)
#else
#define LDCSD_LOG_TRACE(...) ((void)0)
#define LDCSD_LOG_DEBUG(...) ((void)0)
#endif
#define LDCSD_LOG_TRACE_DEFAULT(message) Logger::instance().log(Level::Trace, message)
#define LDCSD_LOG_TRACE_CH(channel, message) Logger::instance().log(channel, Level::Trace, message)
#define LDCSD_LOG_DEBUG_DEFAULT(message) Logger::instance().log(Level::Debug, message)
#define LDCSD_LOG_DEBUG_CH(channel, message) Logger::instance().log(channel, Level::Debug, message)

#define LDCSD_LOG_INFO(...)                                                                        \
  LDCSD_LOG_PICK_ARITY(__VA_ARGS__, LDCSD_LOG_INFO_CH, LDCSD_LOG_INFO_DEFAULT)(__VA_ARGS__)
#define LDCSD_LOG_INFO_DEFAULT(message) Logger::instance().log(Level::Info, message)
#define LDCSD_LOG_INFO_CH(channel, message) Logger::instance().log(channel, Level::Info, message)

#define LDCSD_LOG_WARN(...)                                                                        \
  LDCSD_LOG_PICK_ARITY(__VA_ARGS__, LDCSD_LOG_WARN_CH, LDCSD_LOG_WARN_DEFAULT)(__VA_ARGS__)
#define LDCSD_LOG_WARN_DEFAULT(message) Logger::instance().log(Level::Warn, message)
#define LDCSD_LOG_WARN_CH(channel, message) Logger::instance().log(channel, Level::Warn, message)

#define LDCSD_LOG_ERROR(...)                                                                       \
  LDCSD_LOG_PICK_ARITY(__VA_ARGS__, LDCSD_LOG_ERROR_CH, LDCSD_LOG_ERROR_DEFAULT)(__VA_ARGS__)
#define LDCSD_LOG_ERROR_DEFAULT(message) Logger::instance().log(Level::Error, message)
#define LDCSD_LOG_ERROR_CH(channel, message) Logger::instance().log(channel, Level::Error, message)

// Generic escape hatch when both channel and level are runtime values rather
// than known at the call site (e.g. forwarding).
#define LDCSD_LOG(channel, level, message) Logger::instance().log(channel, level, message)

#endif
