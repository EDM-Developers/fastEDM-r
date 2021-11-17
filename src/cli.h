#include "edm.h"
#include <fstream>
#include <iomanip>
#include <iostream>

#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
#endif
#include <fmt/format.h>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

class ConsoleIO : public IO
{
public:
  ConsoleIO() { this->verbosity = 1; }
  ConsoleIO(int v) { this->verbosity = v; }
  virtual void out(const char* s) const { std::cout << s; }
  virtual void error(const char* s) const { std::cerr << s; }
  virtual void flush() const { fflush(stdout); }
};

struct Inputs
{
  Options opts;
  ManifoldGenerator generator;
  int E;
  std::vector<bool> libraryRows, predictionRows;
};

Inputs parse_lowlevel_inputs_file(const json& j);
Inputs read_lowlevel_inputs_file(std::string fName);
void append_to_dumpfile(std::string fName, const json& taskGroup);
std::vector<bool> int_to_bool(std::vector<int> iv);
json run_tests(json testInputs, int nthreads, IO* io);