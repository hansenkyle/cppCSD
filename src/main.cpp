#include "logger.h"

int main() {
  Logger::configure("ldcsd.log");
  LDCSD_LOG_INFO("ldcsd starting");
  return 0;
}
