#ifndef TERMINAL_H
#define TERMINAL_H

#include <filesystem>
#include <fstream>

#include <ctime>   // timestamp()
#include <sstream> // timestamp()

// Wrap std::out. print() prints to stdout and out file, configured at start of run.

class Terminal {
public:
  static Terminal& instance() {
    static Terminal terminal;
    return terminal;
  }
  // default constructor; default outfile path is current path
  Terminal() {
    // path = std::filesystem::current_path() / "out.txt";
    // outfile.open(path);
  }
  static void configure(std::filesystem::path newpath) {
    Terminal& terminal = instance();
    // close ofstream
    if (terminal.outfile.is_open()) {
      terminal.outfile.close();
    }
    // open ofstream with new path
    terminal.outfile.open(newpath, std::ios::app);
    if (!terminal.outfile.is_open()) {
      throw std::runtime_error("failed to open output file file: " + newpath.string());
    }
    // copy contents of old outfile to new one, if there was one (first-time
    // configure() has no prior file, and reading from an unopened ifstream's
    // rdbuf sets failbit on outfile, silently breaking all future writes)
    std::ifstream in(terminal.path);
    if (in.is_open()) {
      terminal.outfile << in.rdbuf();
      in.close();
      // delete old outfile
      std::filesystem::remove(terminal.path);
    }
    // update this->outfile
    terminal.path = newpath;
  }

  void Print(const std::string& message) {
    outfile << timestamp() << " : " << message << "\n";
    std::cout << message << "\n";
    outfile.flush();
  }
  void PrintSilent(const std::string& message) {
    outfile << timestamp() << " : " << message << "\n";
    outfile.flush();
  }

private:
  std::filesystem::path path;
  std::ofstream outfile;

  std::string timestamp() {
    const std::time_t now = std::time(nullptr);
    const std::tm* tm = std::localtime(&now);
    std::ostringstream oss;
    oss << std::put_time(tm, "%Y-%m-%d %H:%M:%S");
    return oss.str();
  }
};

#define PRINT(...) Terminal::instance().Print(__VA_ARGS__)
#define PRINT_SILENT(...) Terminal::instance().PrintSilent(__VA_ARGS__)

#endif // TERMINAL_H
