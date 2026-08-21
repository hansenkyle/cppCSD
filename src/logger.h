#ifndef LOGGER_H
#define LOGGER_H

#include <filesystem>
#include <fstream>
#include <string>

enum class LogLevel { Trace, Debug, Info, Warn, Error };

const char *to_string(LogLevel level);

// Minimal, serial-only logger. Every call to log() appends one
// "[timestamp] [LEVEL] message" line to the log file and flushes
// immediately -- all levels go to the same file; a separate script is
// expected to split/filter it later. Levels aren't filtered at runtime:
// Trace/Debug calls are instead compiled out entirely unless
// LDCSD_ENABLE_DEBUG_LOGGING is defined (see the LDCSD_LOG_* macros below
// and the LDCSD_ENABLE_DEBUG_LOGGING CMake option).
class Logger {
public:
  explicit Logger(const std::filesystem::path &log_path);

  void log(LogLevel level, const std::string &message);

private:
  std::ofstream stream_;
};

#ifdef LDCSD_ENABLE_DEBUG_LOGGING
#define LDCSD_LOG_TRACE(logger, message) (logger).log(LogLevel::Trace, message)
#define LDCSD_LOG_DEBUG(logger, message) (logger).log(LogLevel::Debug, message)
#else
#define LDCSD_LOG_TRACE(logger, message) ((void) 0)
#define LDCSD_LOG_DEBUG(logger, message) ((void) 0)
#endif

#define LDCSD_LOG_INFO(logger, message) (logger).log(LogLevel::Info, message)
#define LDCSD_LOG_WARN(logger, message) (logger).log(LogLevel::Warn, message)
#define LDCSD_LOG_ERROR(logger, message) (logger).log(LogLevel::Error, message)

#endif
