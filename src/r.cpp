
#define STRICT_R_HEADERS

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppThread)]]

#include <Rcpp.h>
using namespace Rcpp;


#include <RcppThread.h>
#define RCPPTHREAD_OVERRIDE_THREAD 1

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

//#include "common.h"
//#include "edm.h"
#include "cli.h"
#include "cpu.h"

class RConsoleIO : public IO
{
public:
  RConsoleIO() { this->verbosity = 1; }
  RConsoleIO(int v) { this->verbosity = v; }
  virtual void out(const char* s) const { Rcpp::Rcout << s; }
  virtual void error(const char* s) const { Rcpp::Rcerr << s; }
  virtual void flush() const { R_FlushConsole(); }
};

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


#define rListGet(l, key, def) (l.containsElementNamed(key) ? l[key] : def)

Options parse_options(const List &l) {
  Options opts;
  
  opts.nthreads = rListGet(l, "nthreads", 1);

  opts.copredict = rListGet(l, "copredict", false);
  opts.forceCompute = rListGet(l, "forceCompute", true);
  opts.savePrediction = rListGet(l, "savePrediction", true);
  opts.saveSMAPCoeffs = rListGet(l, "saveSMAPCoeffs", false);
  
  opts.k = rListGet(l, "k", 1);
  opts.missingdistance = rListGet(l, "missingdistance", 0);
  opts.dtWeight = rListGet(l, "dtWeight", 0);
  // opts.panelMode = rListGet(l, "panelMode", false);
  opts.idw = rListGet(l, "idw", 0);
  if (l.containsElementNamed("thetas")) {
    opts.thetas = Rcpp::as<std::vector<double>>(l["thetas"]);
  } else {
    opts.thetas = { 1.0 };
  }

  std::string alg = rListGet(l, "algorithm", "simplex");
  if (alg == "simplex") {
    opts.algorithm = Algorithm::Simplex;
  } else if (alg == "smap") {
    opts.algorithm = Algorithm::SMap;
  } else {
    return {}; // TODO
  }

  opts.calcRhoMAE = rListGet(l, "calcRhoMAE", true);
  
  //opts.taskNum = 1;
  //opts.numTasks = 1;
  opts.aspectRatio = 0;
  opts.distance = Distance::Euclidean;
  opts.metrics = {};
  opts.cmdLine = "";
  opts.saveKUsed = false;

  return opts;
}

bool rcpp_keep_going() {
  Rcpp::checkUserInterrupt();
  return !RcppThread::isInterrupted();
}

// [[Rcpp::export]]
std::string run_command(DataFrame df, NumericVector es, NumericVector libs, List ropts,
                        int tau, int p, int numReps=1, int crossfold=0,
                        bool full=false, bool dtMode=false, bool allowMissing=false,
                        int verbosity = 1)
{
  RConsoleIO io(verbosity);
  
  Rcout << "Verbosity set to " << verbosity << "\n";
  
  std::vector<int> Es = Rcpp::as<std::vector<int>>(es);
  std::vector<int> libraries = Rcpp::as<std::vector<int>>(libs);
  
  Options opts = parse_options(ropts);

  Rcout << "Num threads used is " << opts.nthreads << "\n";
  Rcout << "CPU has " << num_logical_cores() << " logical cores and " << num_physical_cores() << " physical cores\n";
  
  Rcout << "k is " << opts.k << "\n";
  
  
  //io->print(fmt::format("[{}] Starting the Stata command: {}\n", taskGroupNum, opts.cmdLine));

  std::vector<double> t = Rcpp::as<std::vector<double>>(df["t"]);
  std::vector<double> x = Rcpp::as<std::vector<double>>(df["x"]);

  //Rcout << "t is " << Rcpp::as<NumericVector>(df["t"]) << "\n";
  //Rcout << "x is " << Rcpp::as<NumericVector>(df["x"]) << "\n";

  bool explore;
  std::vector<double> xmap;
  if (df.containsElementNamed("y")) {
    xmap = Rcpp::as<std::vector<double>>(df["y"]);
    explore = false;
  } else {
    xmap = x;
    explore = true;
  }
  Rcout << "explore mode? " << explore << "\n";
  
  
  std::vector<int> panelIDs;
  if (df.containsElementNamed("id")) {
    panelIDs = Rcpp::as<std::vector<int>>(df["id"]);
    opts.panelMode = true;
  } else {
    opts.panelMode = false;
  }
  

  std::vector<double> co_x = {};

  std::vector<std::vector<double>> extras = {};
  int numExtrasLagged = 0;

  bool reldt = false;

  const ManifoldGenerator generator(t, x, tau, p, xmap, co_x, panelIDs, extras,
                                    numExtrasLagged, dtMode, reldt, allowMissing);

  int maxE = Es[Es.size()-1];
  
  /*
  Manifold M_all = generator.create_manifold(maxE, {}, false, false);
  
  Rcout << "M_all size is (" << M_all.nobs() << ", " << M_all.E_actual() << ")\n";
  
  for (int i = 0; i < M_all.nobs(); i++) {
    for (int j = 0; j < M_all.E_actual(); j++){ 
      Rcout << M_all(i, j) << " ";
    }
    Rcout << "\n";
  }
  */
  
  std::vector<bool> usable = generator.generate_usable(maxE);

  int numUsable = std::accumulate(usable.begin(), usable.end(), 0);
  if (numUsable == 0) {
    Rcout << "Num usable is 0!\n";
    return "";
  }

  int k = opts.k;

  // int numReps = taskGroup["numReps"];
  // int crossfold = taskGroup["crossfold"];
  // bool full = taskGroup["full"];

  bool saveFinalPredictions = opts.savePrediction;
  bool saveSMAPCoeffs = opts.saveSMAPCoeffs;


  bool copredictMode = false; // taskGroup["copredictMode"];

  
  bool saveFinalCoPredictions = false;
  std::vector<bool> coTrainingRows = {}; // = int_to_bool(taskGroup["coTrainingRows"]);
  std::vector<bool> coPredictionRows = {}; // int_to_bool(taskGroup["coPredictionRows"]);

  std::string rngState = ""; // taskGroup["rngState"];

  //opts.numTasks = numReps * crossfold * libraries.size() * Es.size();
  
  io.print("Starting the command!\n");
  io.flush();
  
  bool shuffle = false;
  
  std::vector<std::future<PredictionResult>> futures = launch_task_group(
    generator, opts, Es, libraries, k, numReps, crossfold, explore, full, shuffle,
    saveFinalPredictions, saveFinalCoPredictions, saveSMAPCoeffs,
    copredictMode, usable, rngState, &io, rcpp_keep_going, nullptr);
  
  Rcout << "Waiting for " << futures.size() << " results to come back\n";
  io.flush();
  
  json results;
  json summaryTable;

  // Collect the results of this task group before moving on to the next task group
  int rc = 0;
  
  // futures.size() iterations in loop, update progress every 1 sec
  RcppThread::ProgressBar bar(futures.size(), 1);
  //RcppThread::parallelFor(0, 100, [&] (int i) {
  //  std::this_thread::sleep_for(std::chrono::milliseconds(500));
  //  bar++;
  //});
  
  for (int f = 0; f < futures.size(); f++) {
    const PredictionResult pred = futures[f].get();
    bar++;
    io.print(io.get_and_clear_async_buffer());
    io.flush();

    summaryTable.push_back(pred.stats);
    if (pred.rc > rc) {
      rc = pred.rc;
    }
    
    std::vector<double> predictionsVec, coeffsVec;
    if (pred.predictions != nullptr) {
      predictionsVec = std::vector<double>(pred.predictions.get(),
                                           pred.predictions.get() + pred.numThetas * pred.numPredictions);
      results["predictions"] = predictionsVec;
    }
    if (pred.coeffs != nullptr) {
      coeffsVec = std::vector<double>(pred.coeffs.get(),
                                      pred.coeffs.get() + pred.numPredictions * pred.numCoeffCols);
      results["coeffs"] = coeffsVec;
    }
  }
  
  results["summaryTable"] = summaryTable;
  results["rc"] = rc;
  
  Rcout << "Got them all\n rc";

  Rcout << fmt::format("Return code is {}\n", rc);
  //io->print(fmt::format("Return code is {}\n", rc));

  return results.dump();
}

/*
// [[Rcpp::export]]
void rcpp_launch_command() {
  
  int E = 2;
  int tau = 1;
  int p = 1;
  
  std::vector<double> t = { 1, 2, 3, 4 };
  std::vector<double> x = { 11, 12, 13, 14 };

  ManifoldGenerator generator(t, x, tau, p);
  
  std::vector<bool> usable = generator.generate_usable(E);
  
  Options opts;
  
  opts.copredict = false;
  opts.forceCompute = true;
  opts.savePrediction = true;
  opts.saveSMAPCoeffs = false;
  
  opts.k = 1;
  opts.nthreads = 1;
  opts.missingdistance = 0;
  opts.dtWeight = 0;
  opts.panelMode = false;
  opts.idw = 0;
  opts.thetas = { 1.0 };
  opts.algorithm = Algorithm::Simplex;
  opts.taskNum = 1;
  opts.numTasks = 1;
  opts.calcRhoMAE = false;
  opts.aspectRatio = 0;
  opts.distance = Distance::Euclidean;
  opts.metrics = {};
  opts.cmdLine = "edm explore x";
  opts.saveKUsed = false;

  const std::vector<int> Es = { 2 };
  const std::vector<int> libraries = { 2 };
  
  int k = 1; // TODO: Why is this in opt & here also?
  
  int numReps = 1;
  int crossfold = 0;
  bool explore = true;
  bool full = true;
  bool saveFinalPredictions = true; // TODO: Same as k
  bool saveSMAPCoeffs = false;
  bool copredictMode = false;
  
  const std::vector<bool> coTrainingRows = {};
  const std::vector<bool> coPredictionRows = {};
  const std::string rngState = "";
  
  IO* io = nullptr;
  //bool keep_going() = nullptr;
  //void all_tasks_finished() = nullptr;
  
  std::vector<std::future<Prediction>> res = launch_task_group(
      generator, opts, Es, libraries, k, numReps, crossfold, explore, full,
      saveFinalPredictions, saveSMAPCoeffs, copredictMode, usable,
      coTrainingRows, coPredictionRows, rngState,
      io, nullptr, nullptr);
  
  
}
*/

// [[Rcpp::export]]
std::string run_json_test(std::string fnameIn) {
  
  //std::string fnameIn("easy.json");
  
  int nthreads = 1;

  int verbosity = 1;
  ConsoleIO io(verbosity);
  
  std::ifstream i(fnameIn);
  json testInputs;
  i >> testInputs;
  
  json results = run_tests(testInputs, nthreads, &io);
  
  return results.dump();
  
  /*  
  std::vector<double> rhos;
  
  for (int r = 0; r < results.size(); r++) {
    json pJS = results[r];
    //Rcout << pJS << "\n";
    Prediction p = pJS;

    for (int s = 0; s < p.stats.size(); s++) {
      rhos.push_back(p.stats[s].rho);
    }
    
    // Crashes, as p.numPredictions not necessarily the size of p.ystar, bug?
    / *    
    if(p.numPredictions > 0) {
      std::vector<double> preds(p.numPredictions);
      
      for (int i = 0; i < p.numPredictions; i++) {
        preds[i] = p.ystar[i];
      }
      
      return preds;
    }
     * /
  }
  
  return rhos;
  */
}

/**  R
run_json_test("easy.json")

*/
