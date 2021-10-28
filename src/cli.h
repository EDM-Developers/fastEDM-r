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
  std::vector<bool> trainingRows, predictionRows;
};

/*! \brief Read in a dump file.
 *
 * Read in a dump file created with the Stata `saveinputs' option.
 *
 * \param fname dump filename
 * \param pointer to InputVars struct to store the read
 */
Inputs parse_lowlevel_inputs_file(const json& j)
{
  int E = j["E"];
  Options opts = j["opts"];
  ManifoldGenerator generator = j["generator"];
  std::vector<bool> trainingRows = j["trainingRows"], predictionRows = j["predictionRows"];

  return { opts, generator, E, trainingRows, predictionRows };
}

Inputs read_lowlevel_inputs_file(std::string fName)
{
  std::ifstream i(fName);
  json j;
  i >> j;

  return parse_lowlevel_inputs_file(j);
}

void append_to_dumpfile(std::string fName, const json& taskGroup)
{
  json allTaskGroups;

  std::ifstream i(fName);
  if (i.is_open()) {
    i >> allTaskGroups;
  }

  allTaskGroups.push_back(taskGroup);

  // Add "o << std::setw(4) << allTaskGroups" to pretty-print the saved JSON
  std::ofstream o(fName);
  o << allTaskGroups << std::endl;
}

std::vector<bool> int_to_bool(std::vector<int> iv)
{
  std::vector<bool> bv;
  std::copy(iv.begin(), iv.end(), std::back_inserter(bv));
  return bv;
}

json run_tests(json testInputs, int nthreads, IO* io)
{
  int rc = 0;
  json results;

  int numTaskGroups = testInputs.size();

  io->print(fmt::format("Number of tests in this JSON file is {}\n", numTaskGroups));

  for (int taskGroupNum = 0; taskGroupNum < numTaskGroups; taskGroupNum++) {
    json taskGroup = testInputs[taskGroupNum];

    Options opts = taskGroup["opts"];
    opts.nthreads = nthreads;

    io->print(fmt::format("[{}] Starting the Stata command: {}\n", taskGroupNum, opts.cmdLine));

    ManifoldGenerator generator = taskGroup["generator"];
    std::vector<int> Es = taskGroup["Es"];
    std::vector<int> libraries = taskGroup["libraries"];
    int k = taskGroup["k"];
    int numReps = taskGroup["numReps"];
    int crossfold = taskGroup["crossfold"];
    bool explore = taskGroup["explore"];
    bool full = taskGroup["full"];
    bool saveFinalPredictions = taskGroup["saveFinalPredictions"];
    bool saveFinalCoPredictions = taskGroup["saveFinalCoPredictions"];
    bool saveSMAPCoeffs = taskGroup["saveSMAPCoeffs"];
    bool copredictMode = taskGroup["copredictMode"];
    std::vector<bool> usable = int_to_bool(taskGroup["usable"]);
    std::string rngState = taskGroup["rngState"];

    std::vector<std::future<Prediction>> futures =
      launch_task_group(generator, opts, Es, libraries, k, numReps, crossfold, explore, full, saveFinalPredictions,
                        saveFinalCoPredictions, saveSMAPCoeffs, copredictMode, usable, rngState, io, nullptr, nullptr);

    // Collect the results of this task group before moving on to the next task group
    for (int f = 0; f < futures.size(); f++) {
      const Prediction pred = futures[f].get();
      io->print(io->get_and_clear_async_buffer());
      io->flush();

      results.push_back(pred);
      if (pred.rc > rc) {
        rc = pred.rc;
      }
    }
  }

  io->print(fmt::format("Return code is {}\n", rc));

  return results;
}
