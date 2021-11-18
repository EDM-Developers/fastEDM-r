
#define STRICT_R_HEADERS

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppThread)]]

#include <Rcpp.h>
using namespace Rcpp;

#include <RcppThread.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

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

bool rcpp_keep_going() {
  // The following two calls cause all kinds of crashes.
  //Rcpp::checkUserInterrupt();
  //return !RcppThread::isInterrupted();
  return 1;
}

void replace_nan(std::vector<double>& v) {
  for (int i = 0; i < v.size(); i++) {
    if (!std::isnormal(v[i])) {
      v[i] = MISSING_D;
    }
  }
}

// [[Rcpp::export]]
List run_command(DataFrame df, IntegerVector es,
                        int tau, NumericVector thetas, Nullable<IntegerVector> libs,
                        int k = 0, std::string algorithm = "simplex", int numReps = 1,
                        int p = 1, int crossfold = 0, bool full = false, bool shuffle = false,
                        bool saveFinalPredictions=false, bool saveSMAPCoeffs=false,
                        bool dt = false, bool allowMissing = false, int numThreads = 1, int verbosity = 1)
{
  RConsoleIO io(verbosity);
  
  io.print(fmt::format("Verbosity set to {}\n", verbosity));

  Options opts;

  opts.nthreads = numThreads;
  opts.copredict = false;
  opts.forceCompute = true;
  opts.savePrediction = saveFinalPredictions;
  opts.saveSMAPCoeffs = saveSMAPCoeffs;

  opts.k = k;
  opts.missingdistance = 0.0;
  opts.dtWeight = 0.0;
  opts.panelMode = false;
  opts.idw = 0;

  if (thetas.size() > 0) {
      opts.thetas = Rcpp::as<std::vector<double>>(thetas);
  } else {
      opts.thetas = { 1.0 };
  }


  std::vector<int> Es = Rcpp::as<std::vector<int>>(es);

  std::vector<int> libraries;
  if (libs.isNotNull()) {
    libraries = Rcpp::as<std::vector<int>>(libs);
  }
  if (algorithm == "simplex") {
      opts.algorithm = Algorithm::Simplex;
  } else if (algorithm == "smap") {
      opts.algorithm = Algorithm::SMap;
  } else {
      return {}; // TODO
  }

  opts.calcRhoMAE = true; // TODO: When is this off?

  opts.aspectRatio = 0;
  opts.distance = Distance::Euclidean;
  opts.metrics = {};
  opts.cmdLine = "";
  opts.saveKUsed = true; // TODO: Check again

  io.print(fmt::format("Num threads used is {}\n", opts.nthreads));
  io.print(fmt::format("CPU has {} logical cores and {} physical cores\n", num_logical_cores(), num_physical_cores()));
  io.print(fmt::format("k is {}\n", k));

  std::vector<double> t = Rcpp::as<std::vector<double>>(df["t"]);
  std::vector<double> x = Rcpp::as<std::vector<double>>(df["x"]);
  
  replace_nan(t);
  replace_nan(x);

  //Rcout << "t is " << Rcpp::as<NumericVector>(df["t"]) << "\n";
  //Rcout << "x is " << Rcpp::as<NumericVector>(df["x"]) << "\n";

  bool explore;
  std::vector<double> xmap;
  if (df.containsElementNamed("y")) {
    xmap = Rcpp::as<std::vector<double>>(df["y"]);
    replace_nan(xmap);
    explore = false;
  } else {
    xmap = x;
    explore = true;
  }
  
  io.print(fmt::format("explore mode is {}\n", explore));
  
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
                                    numExtrasLagged, dt, reldt, allowMissing);

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
    io.print("Num usable is 0!\n");
    return "";
  }
  
  if (libraries.size() == 0) {
    std::vector<bool> usable = generator.generate_usable(maxE);
    int numUsable = std::accumulate(usable.begin(), usable.end(), 0);
    libraries = { numUsable };
  }

  bool copredictMode = false; // taskGroup["copredictMode"];
  bool saveFinalCoPredictions = false;
  std::vector<bool> coTrainingRows = {}; // = int_to_bool(taskGroup["coTrainingRows"]);
  std::vector<bool> coPredictionRows = {}; // int_to_bool(taskGroup["coPredictionRows"]);

  std::string rngState = ""; // taskGroup["rngState"];

  //opts.numTasks = numReps * crossfold * libraries.size() * Es.size();
  
  io.print("Starting the command!\n");
  io.flush();
  
  std::vector<std::future<PredictionResult>> futures = launch_task_group(
    generator, opts, Es, libraries, k, numReps, crossfold, explore, full, shuffle,
    saveFinalPredictions, saveFinalCoPredictions, saveSMAPCoeffs,
    copredictMode, usable, rngState, &io, rcpp_keep_going, nullptr);
  
  io.print(fmt::format("Waiting for {} results to come back\n", futures.size()));
  io.flush();

  int rc = 0;
  
  RcppThread::ProgressBar bar(futures.size(), 1);
  
  int kMin, kMax;
  
  std::vector<double> predictionsVec, coeffsVec;
  DataFrame summary;
  
  { 
    IntegerVector Es, libraries;
    NumericVector thetas, rhos, maes;
    
    auto Rint = [](double v) { return (v != MISSING_D) ? v : NA_INTEGER; };
    auto Rdouble = [](double v) { return (v != MISSING_D) ? v : NA_REAL; };
    
    for (int f = 0; f < futures.size(); f++) {
      const PredictionResult pred = futures[f].get();
      bar++;
      
      if (f == 0 || pred.kMin < kMin) {
        kMin = pred.kMin;
      }
      if (f == 0 || pred.kMax > kMax) {
        kMax = pred.kMax;
      }
      
      //io.print(io.get_and_clear_async_buffer());
      //io.flush();
      
      for (int t = 0; t < pred.stats.size(); t++) {
        Es.push_back(Rint(pred.stats[t].E));
        thetas.push_back(Rdouble(pred.stats[t].theta));
        libraries.push_back(Rint(pred.stats[t].library));
        rhos.push_back(Rdouble(pred.stats[t].rho));
        maes.push_back(Rdouble(pred.stats[t].mae));
      }
  
      if (pred.rc > rc) {
        rc = pred.rc;
      }
      
      if (pred.predictions != nullptr) {
        predictionsVec = std::vector<double>(pred.predictions.get(),
                                             pred.predictions.get() + pred.numThetas * pred.numPredictions);
      }
      if (pred.coeffs != nullptr) {
        coeffsVec = std::vector<double>(pred.coeffs.get(),
                                        pred.coeffs.get() + pred.numPredictions * pred.numCoeffCols);
      }
    }
    
    summary = Rcpp::DataFrame::create(_["E"] = Es, _["library"] = libraries,
                                      _["theta"] = thetas, _["rho"] = rhos,
                                      _["mae"] = maes);
  }
  
  //Rcpp::NumericVector summaryTable = Rcpp::wrap(summary);
  //summaryTable.attr("dim") = Rcpp::Dimension(summary.size() / 5, 5);
  
  //results["summaryTable"] = summaryTable;
  
  //results["rc"] = rc;
  
  io.print(fmt::format("k value was between {} and {}\n", kMin, kMax));
  
  io.print(fmt::format("Return code is {}\n", rc));

  return List::create(_["summary"] = summary,
                      _["predictions"] = predictionsVec,
                      _["coeffs"] = coeffsVec,
                      _["rc"] = rc,
                      _["kMin"] = kMin,
                      _["kMax"] = kMax);
}

// [[Rcpp::export]]
std::string run_json_test(std::string fnameIn) {
  
  int nthreads = 1;
  int verbosity = 1;
  
  ConsoleIO io(verbosity);
  
  std::ifstream i(fnameIn);
  json testInputs;
  i >> testInputs;
  
  json results = run_tests(testInputs, nthreads, &io);
  
  return results.dump();
}

/**  R
run_json_test("easy.json")
*/
