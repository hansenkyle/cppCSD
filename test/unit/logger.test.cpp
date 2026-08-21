#include "logger.h"

#include <fstream>
#include <sstream>

#include <doctest.h>

namespace {

std::filesystem::path tempLogPath(const std::string& name) {
  return std::filesystem::temp_directory_path() / ("ldcsd_logger_test_" + name + ".log");
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
    CHECK(std::string(to_string(LogLevel::Trace)) == "TRACE");
    CHECK(std::string(to_string(LogLevel::Debug)) == "DEBUG");
    CHECK(std::string(to_string(LogLevel::Info)) == "INFO");
    CHECK(std::string(to_string(LogLevel::Warn)) == "WARN");
    CHECK(std::string(to_string(LogLevel::Error)) == "ERROR");
  }

  TEST_CASE("log() writes a timestamped, leveled line to the log file") {
    const auto path = tempLogPath("basic");
    {
      Logger logger(path);
      logger.log(LogLevel::Info, "hello");
    }

    const std::string contents = readFile(path);
    CHECK(contents.find("[INFO] hello") != std::string::npos);
    // "[YYYY-MM-DD HH:MM:SS] " prefix.
    CHECK(contents.substr(0, 1) == "[");
    CHECK(contents.find("] [INFO] hello") != std::string::npos);
  }

  TEST_CASE("each call to log() appends and flushes immediately") {
    const auto path = tempLogPath("appends");
    Logger logger(path);

    logger.log(LogLevel::Warn, "first");
    CHECK(readFile(path).find("[WARN] first") != std::string::npos);

    logger.log(LogLevel::Error, "second");
    const std::string contents = readFile(path);
    CHECK(contents.find("[WARN] first") != std::string::npos);
    CHECK(contents.find("[ERROR] second") != std::string::npos);
  }

  TEST_CASE("constructor throws if the log file can't be opened") {
    const std::filesystem::path bad_path =
        tempLogPath("missing-dir-does-not-exist") / "sub" / "log";
    CHECK_THROWS_AS(Logger logger(bad_path), std::runtime_error);
  }

  TEST_CASE("LDCSD_LOG_INFO/WARN/ERROR always write") {
    const auto path = tempLogPath("macros-always-on");
    Logger logger(path);

    LDCSD_LOG_INFO(logger, "info message");
    LDCSD_LOG_WARN(logger, "warn message");
    LDCSD_LOG_ERROR(logger, "error message");

    const std::string contents = readFile(path);
    CHECK(contents.find("[INFO] info message") != std::string::npos);
    CHECK(contents.find("[WARN] warn message") != std::string::npos);
    CHECK(contents.find("[ERROR] error message") != std::string::npos);
  }

  TEST_CASE("LDCSD_LOG_TRACE/DEBUG are compiled out unless LDCSD_ENABLE_DEBUG_LOGGING is set") {
    const auto path = tempLogPath("macros-trace-debug");
    Logger logger(path);

    LDCSD_LOG_TRACE(logger, "trace message");
    LDCSD_LOG_DEBUG(logger, "debug message");

    const std::string contents = readFile(path);
#ifdef LDCSD_ENABLE_DEBUG_LOGGING
    CHECK(contents.find("[TRACE] trace message") != std::string::npos);
    CHECK(contents.find("[DEBUG] debug message") != std::string::npos);
#else
    CHECK(contents.empty());
#endif
  }
}
