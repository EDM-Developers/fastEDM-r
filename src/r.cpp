
#define STRICT_R_HEADERS

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
using namespace Rcpp;

#include <RcppEigen.h>

#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
#endif
#include <fmt/format.h>

#include "edm.h"
#include "cpu.h"
#include "stats.h"

class RConsoleIO : public IO
{
public:
  RConsoleIO() { this->verbosity = 1; }
  RConsoleIO(int v) { this->verbosity = v; }
  virtual void out(const char* s) const { Rcpp::Rcout << s; }
  virtual void error(const char* s) const { Rcpp::Rcerr << s; }
  virtual void flush() const { R_FlushConsole(); }
};

bool rcpp_keep_going()
{
  // The following two calls cause all kinds of crashes.
  // Rcpp::checkUserInterrupt();
  // return !RcppThread::isInterrupted();
  return 1;
}

void replace_nan(std::vector<double>& v)
{
  for (int i = 0; i < v.size(); i++) {
    if (!std::isnormal(v[i])) {
      v[i] = MISSING_D;
    }
  }
}

Rcpp::NumericMatrix to_R_matrix(const double* v, int r, int c,
                                std::vector<bool> filter = {},
                                bool rowMajor = false)
{
  Rcpp::NumericMatrix mat(r, c);

  int obsNum = 0;

  if (rowMajor) {
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < c; j++) {
        if (filter.size() > 0 && !filter[i]) {
          mat(i, j) = NA_REAL;
          continue;
        }
        mat(i, j) = v[obsNum] != MISSING_D ? v[obsNum] : NA_REAL;
        obsNum += 1;
      }
    }
  } else {
    for (int j = 0; j < c; j++) {
      for (int i = 0; i < r; i++) {
        if (filter.size() > 0 && !filter[i]) {
          mat(i, j) = NA_REAL;
          continue;
        }

        mat(i, j) = v[obsNum] != MISSING_D ? v[obsNum] : NA_REAL;
        obsNum += 1;
      }
    }
  }

  return mat;
}

// [[Rcpp::export]]
List run_command(DataFrame df, IntegerVector es, int tau, NumericVector thetas, Nullable<IntegerVector> libs, int k = 0,
                 std::string algorithm = "simplex", int numReps = 1, int p = 1, int crossfold = 0, bool full = false,
                 bool shuffle = false, bool saveFinalPredictions = false, bool saveFinalCoPredictions = false,
                 bool saveManifolds = false, bool saveSMAPCoeffs = false, bool dt = false, bool reldt = false,
                 double dtWeight = 0.0, Nullable<List> extras = R_NilValue, bool allowMissing = false,
                 double missingDistance = 0.0, int numThreads = 1, int verbosity = 1)
{
  RConsoleIO io(verbosity);

  io.print(fmt::format("Verbosity set to {}\n", verbosity));

  Options opts;

  opts.nthreads = numThreads;
  opts.copredict = false;
  opts.forceCompute = true;
  opts.saveManifolds = saveManifolds;
  opts.saveSMAPCoeffs = saveSMAPCoeffs;

  opts.k = k;
  opts.missingdistance = missingDistance;
  
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

  // Rcout << "t is " << Rcpp::as<NumericVector>(df["t"]) << "\n";
  // Rcout << "x is " << Rcpp::as<NumericVector>(df["x"]) << "\n";

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

  io.print(fmt::format("panel data mode is {}\n", opts.panelMode));

  std::vector<double> co_x;

  if (df.containsElementNamed("co_x")) {
    co_x = Rcpp::as<std::vector<double>>(df["co_x"]);
    replace_nan(co_x);
  }

  std::vector<std::vector<double>> extrasVecs;

  if (extras.isNotNull()) {
    List extrasList = Rcpp::as<List>(extras);

    io.print(fmt::format("Num extras is {}\n", extrasList.size()));
    for (int e = 0; e < extrasList.size(); e++) {
      extrasVecs.emplace_back(Rcpp::as<std::vector<double>>(extrasList[e]));
      replace_nan(extrasVecs[extrasVecs.size() - 1]);
      opts.metrics.push_back(Metric::Diff); // TODO: Handle factor extras
    }
  }

  // for (int e = 0; e < extrasVecs.size(); e++) {
  //   io.print(fmt::format("Extra[{}] len {}: ", e, extrasVecs[e].size()));
  //   for (int i = 0; i < extrasVecs[e].size(); i++) {
  //     io.print(fmt::format("{} ", extrasVecs[e][i]));
  //   }
  //   io.print("\n");
  // }

  int numExtrasLagged = 0;

  ManifoldGenerator generator(t, x, tau, p, xmap, co_x, panelIDs, extrasVecs,
                              numExtrasLagged, dt, reldt, allowMissing,
                              dtWeight);

  // Manifold M_all = generator.create_manifold(maxE, {}, true);
  //
  // Rcout << "M_all size is (" << M_all.numPoints() << ", " << M_all.E_actual() << ")\n";
  //
  // for (int i = 0; i < M_all.numPoints(); i++) {
  //   for (int j = 0; j < M_all.E_actual(); j++){
  //     Rcout << M_all(i, j) << " ";
  //   }
  //   Rcout << "\n";
  // }

  if (allowMissing && opts.missingdistance == 0) {
    opts.missingdistance = default_missing_distance(x);
  }

  int maxE = Es[Es.size() - 1];
  std::vector<bool> usable = generator.generate_usable(maxE);

  int numUsable = std::accumulate(usable.begin(), usable.end(), 0);
  if (numUsable == 0) {
    io.print("Num usable is 0!\n");
    return "";
  }

  if (libraries.size() == 0) {
    libraries = { numUsable };
  }

  bool copredictMode = co_x.size() > 0;
  std::string rngState = ""; // taskGroup["rngState"];

  // opts.numTasks = numReps * crossfold * libraries.size() * Es.size();

  io.print("Starting the command!\n");
  io.flush();

  auto genPtr = std::shared_ptr<ManifoldGenerator>(&generator, [](ManifoldGenerator*) {});
  opts.lowMemoryMode = false;

  std::vector<std::future<PredictionResult>> futures = launch_task_group(
    genPtr, opts, Es, libraries, k, numReps, crossfold, explore, full, shuffle, saveFinalPredictions,
    saveFinalCoPredictions, saveSMAPCoeffs, copredictMode, usable, rngState, &io, rcpp_keep_going, nullptr);

  io.print(fmt::format("Waiting for {} results to come back\n", futures.size()));
  io.flush();

  int rc = 0;

  // RcppThread::ProgressBar bar(futures.size(), 1);

  int kMin, kMax;

  NumericMatrix predictions, coPredictions, coeffs;
  DataFrame summary, co_summary;
  std::vector<NumericMatrix> Ms, Mps;

  {
    IntegerVector Es, libraries;
    NumericVector thetas, rhos, maes;

    IntegerVector co_Es, co_libraries;
    NumericVector co_thetas, co_rhos, co_maes;

    auto Rint = [](double v) { return (v != MISSING_D) ? v : NA_INTEGER; };
    auto Rdouble = [](double v) { return (v != MISSING_D) ? v : NA_REAL; };

    for (int f = 0; f < futures.size(); f++) {
      const PredictionResult pred = futures[f].get();
      // bar++;

      if (f == 0 || pred.kMin < kMin) {
        kMin = pred.kMin;
      }
      if (f == 0 || pred.kMax > kMax) {
        kMax = pred.kMax;
      }

      // io.print(io.get_and_clear_async_buffer());
      // io.flush();

      if (!pred.copredict) {
        for (int t = 0; t < pred.stats.size(); t++) {
          Es.push_back(Rint(pred.stats[t].E));
          thetas.push_back(Rdouble(pred.stats[t].theta));
          libraries.push_back(Rint(pred.stats[t].library));
          rhos.push_back(Rdouble(pred.stats[t].rho));
          maes.push_back(Rdouble(pred.stats[t].mae));
        }
      } else {
        for (int t = 0; t < pred.stats.size(); t++) {
          co_Es.push_back(Rint(pred.stats[t].E));
          co_thetas.push_back(Rdouble(pred.stats[t].theta));
          co_libraries.push_back(Rint(pred.stats[t].library));
          co_rhos.push_back(Rdouble(pred.stats[t].rho));
          co_maes.push_back(Rdouble(pred.stats[t].mae));
        }
      }

      if (pred.rc > rc) {
        rc = pred.rc;
      }

      if (pred.predictions != nullptr) {
        if (!pred.copredict) {
          predictions = to_R_matrix(pred.predictions.get(),
                                    pred.predictionRows.size(),
                                    pred.numThetas,
                                    pred.predictionRows);
        } else {
          coPredictions = to_R_matrix(pred.predictions.get(),
                                      pred.predictionRows.size(),
                                      pred.numThetas,
                                      pred.predictionRows);
        }
      }
      if (pred.coeffs != nullptr) {
        coeffs = to_R_matrix(pred.coeffs.get(), pred.predictionRows.size(),
                             pred.numCoeffCols, pred.predictionRows);
      }
      
      if (saveManifolds) {
        Ms.push_back(to_R_matrix(pred.M->data(), pred.M->numPoints(), pred.M->E_actual(), {}, true));
        Mps.push_back(to_R_matrix(pred.Mp->data(), pred.Mp->numPoints(), pred.Mp->E_actual(), {}, true));
      }
    }

    summary = Rcpp::DataFrame::create(_["E"] = Es, _["library"] = libraries, _["theta"] = thetas, _["rho"] = rhos,
                                      _["mae"] = maes);

    if (copredictMode) {
      co_summary = Rcpp::DataFrame::create(_["E"] = co_Es, _["library"] = co_libraries, _["theta"] = co_thetas,
                                           _["rho"] = co_rhos, _["mae"] = co_maes);
    }
  }

  io.print(fmt::format("k value was between {} and {}\n", kMin, kMax));
  io.print(fmt::format("Return code is {}\n", rc));

  List res;
  res["summary"] = summary;
  res["rc"] = rc;
  res["kMin"] = kMin;
  res["kMax"] = kMax;
  
  if (copredictMode) {
    res["co_summary"] = co_summary;
  }
  
  if (saveFinalPredictions) {
    res["predictions"] = predictions;
  }
  
  if (saveFinalCoPredictions) {
    res["copredictions"] = coPredictions;
  }
  
  if (saveManifolds) {
    res["Ms"] = Ms;
    res["Mps"] = Mps;
  }
  
  if (saveSMAPCoeffs) {
    res["coeffs"] = coeffs;
  }
  
  if (allowMissing) {
    res["missingdistance"] = opts.missingdistance;
  }
  
  if (dt || reldt) {
    res["dtWeight"] = generator.dtWeight();
  }
  
  return res;
}