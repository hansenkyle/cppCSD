#include "logger.h"

#include <ctime>
#include <iomanip>
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

Logger::Logger(const std::filesystem::path& log_path) : stream_(log_path) {
  if (!stream_.is_open()) {
    throw std::runtime_error("failed to open log file: " + log_path.string());
  }
}

void Logger::log(LogLevel level, const std::string& message) {
  stream_ << "[" << timestamp() << "] [" << to_string(level) << "] " << message << "\n";
  stream_.flush();
}
