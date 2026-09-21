// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "logger.h"

#include <fstream>
#include <sstream>

#include <doctest.h>

namespace {

std::filesystem::path tempLogPath(const std::string& name) {
  return std::filesystem::temp_directory_path() / ("ldcsd_logger_test_" + name + ".log");
}

std::filesystem::path tempOutPath(const std::string& name) {
  return std::filesystem::temp_directory_path() / ("ldcsd_logger_test_" + name + ".out");
}

std::string readFile(const std::filesystem::path& path) {
  std::ifstream in(path);
  std::ostringstream oss;
  oss << in.rdbuf();
  return oss.str();
}

} // namespace

TEST_SUITE("Logger") {
  TEST_CASE("to_string names every level") {
    CHECK(std::string(to_string(Level::Trace)) == "TRACE");
    CHECK(std::string(to_string(Level::Debug)) == "DEBUG");
    CHECK(std::string(to_string(Level::Info)) == "INFO");
    CHECK(std::string(to_string(Level::Warn)) == "WARN");
    CHECK(std::string(to_string(Level::Error)) == "ERROR");
  }

  TEST_CASE("log() writes a timestamped, channel- and level-tagged line to the log file") {
    const auto path = tempLogPath("basic");
    Logger::configure(path);
    Logger::instance().log(Channel::General, Level::Info, "hello");

    const std::string contents = readFile(path);
    CHECK(contents.find("[INFO] [GENERAL] hello") != std::string::npos);
    // "[YYYY-MM-DD HH:MM:SS] " prefix.
    CHECK(contents.substr(0, 1) == "[");
    CHECK(contents.find("] [INFO] [GENERAL] hello") != std::string::npos);
  }

  TEST_CASE("each call to log() appends and flushes immediately") {
    const auto path = tempLogPath("appends");
    Logger::configure(path);

    Logger::instance().log(Channel::General, Level::Warn, "first");
    CHECK(readFile(path).find("[WARN] [GENERAL] first") != std::string::npos);

    Logger::instance().log(Channel::General, Level::Error, "second");
    const std::string contents = readFile(path);
    CHECK(contents.find("[WARN] [GENERAL] first") != std::string::npos);
    CHECK(contents.find("[ERROR] [GENERAL] second") != std::string::npos);
  }

  TEST_CASE("configure() throws if the log file can't be opened") {
    const std::filesystem::path bad_path =
        tempLogPath("missing-dir-does-not-exist") / "sub" / "log";
    CHECK_THROWS_AS(Logger::configure(bad_path), std::runtime_error);
  }

  TEST_CASE("LDCSD_LOG_INFO/WARN/ERROR always write") {
    const auto path = tempLogPath("macros-always-on");
    Logger::configure(path);

    LDCSD_LOG_INFO("info message");
    LDCSD_LOG_WARN("warn message");
    LDCSD_LOG_ERROR("error message");

    const std::string contents = readFile(path);
    CHECK(contents.find("[INFO] [GENERAL] info message") != std::string::npos);
    CHECK(contents.find("[WARN] [GENERAL] warn message") != std::string::npos);
    CHECK(contents.find("[ERROR] [GENERAL] error message") != std::string::npos);
  }

  TEST_CASE("LDCSD_LOG_TRACE/DEBUG are compiled out unless LDCSD_ENABLE_DEBUG_LOGGING is set") {
    const auto path = tempLogPath("macros-trace-debug");
    Logger::configure(path);

    LDCSD_LOG_TRACE("trace message");
    LDCSD_LOG_DEBUG("debug message");

    const std::string contents = readFile(path);
#ifdef LDCSD_ENABLE_DEBUG_LOGGING
    CHECK(contents.find("[TRACE] [GENERAL] trace message") != std::string::npos);
    CHECK(contents.find("[DEBUG] [GENERAL] debug message") != std::string::npos);
#else
    CHECK(contents.empty());
#endif
  }

  TEST_CASE("LDCSD_LOG_INFO(channel, message) tags the line with that channel") {
    const auto path = tempLogPath("macros-channel");
    Logger::configure(path);

    LDCSD_LOG_INFO("default-channel message");
    LDCSD_LOG_INFO(Channel::Iteration, "iteration message");

    const std::string contents = readFile(path);
    CHECK(contents.find("[INFO] [GENERAL] default-channel message") != std::string::npos);
    CHECK(contents.find("[INFO] [ITERATION] iteration message") != std::string::npos);
  }
}
