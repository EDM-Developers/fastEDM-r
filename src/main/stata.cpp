#pragma warning(disable : 4018)

#include "cpu.h"
#include "edm.h"
#include "stplugin.h"

#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
#endif
#include <fmt/format.h>

#include <algorithm>
#include <future>
#include <numeric> // for std::accumulate
#include <stdexcept>
#include <string>
#include <vector>

#include "cli.h" // to save the inputs to a local file for debugging
#include "stats.h"

// These are all the variables in the edm.ado script we modify in the plugin.
// These definitions also suppress the "C++ doesn't permit string literals as char*" warnings.
char* FINISHED_SCALAR = (char*)"plugin_finished";
char* MISSING_DISTANCE_USED = (char*)"_missing_dist_used";
char* DTW_USED = (char*)"_dtw_used";
char* NUM_USABLE = (char*)"_num_usable";
char* K_MIN = (char*)"_k_min";
char* K_MAX = (char*)"_k_max";

// These are all the variables in the edm.ado script we will read from in the plugin.
char* RNG_STATE = (char*)"_rngstate";
char* NUM_NEIGHBOURS = (char*)"_k";
char* SAVE_DT = (char*)"_dtsave";
char* SAVE_MANIFOLD = (char*)"_savemanifold";
char* SAVE_PREDICTION = (char*)"_predictionsave";
char* SAVE_COPREDICTION = (char*)"_copredictionsave";
char* SAVE_SMAP = (char*)"_savesmap";
char* SAVE_INPUTS = (char*)"_saveinputs";
char* NUM_REPS = (char*)"_round";

class StataIO : public IO
{
public:
  virtual void out(const char* s) const { SF_display((char*)s); }
  virtual void error(const char* s) const { SF_error((char*)s); }
  virtual void flush() const { _stata_->spoutflush(); }
};

// Global state, needed to persist between multiple edm calls
StataIO io;
std::atomic<bool> breakButtonPressed = false;
std::atomic<bool> allTasksFinished = false;

bool keep_going()
{
  return !breakButtonPressed;
}

void all_tasks_finished()
{
  allTasksFinished = true;
}

void print_error(std::string command, ST_retcode rc)
{
  // Don't print header if rc=SUCCESS or rc=1 (when Break button pressed)
  if (rc > 1 && io.verbosity > 1) {
    io.error((char*)fmt::format("Error in edm '{}': ", command).c_str());
  }
  switch (rc) {
    case TOO_FEW_VARIABLES:
      io.error("Too few arguments\n");
      break;
    case TOO_MANY_VARIABLES:
      io.error("Too many arguments\n");
      break;
    case NOT_IMPLEMENTED:
      io.error("Method is not yet implemented\n");
      break;
    case INSUFFICIENT_UNIQUE:
      io.error("Insufficient number of unique observations, consider "
               "tweaking the values of E, k or use -force- option\n");
      break;
    case INVALID_ALGORITHM:
      io.error("Invalid algorithm argument\n");
      break;
    case INVALID_DISTANCE:
      io.error("Invalid distance argument\n");
      break;
    case INVALID_METRICS:
      io.error("Invalid metrics argument\n");
      break;
    case UNKNOWN_ERROR:
      io.error("Unknown error\n");
      break;
  }
}

/*
 * Count the number of rows that aren't being filtered out
 * by Stata's 'if' or 'in' expressions.
 */
int num_if_in_rows()
{
  int num = 0;
  for (ST_int i = SF_in1(); i <= SF_in2(); i++) {
    if (SF_ifobs(i)) {
      num += 1;
    }
  }
  return num;
}

/*
 * Read in columns from Stata (i.e. what Stata calls variables).
 * Starting from column number 'j0', read in 'numCols' of columns
 */
template<typename T>
std::vector<T> stata_columns(ST_int j0, int numCols = 1)
{
  // Allocate space for the matrix of data from Stata
  int numRows = num_if_in_rows();
  std::vector<T> M(numRows * numCols);
  int ind = 0; // Flattened index of M matrix
  int r = 0;   // Count each row that isn't filtered by Stata 'if'
  for (ST_int i = SF_in1(); i <= SF_in2(); i++) {
    if (SF_ifobs(i)) { // Skip rows according to Stata's 'if'
      for (ST_int j = j0; j < j0 + numCols; j++) {
        ST_double value;
        ST_retcode rc = SF_vdata(j, i, &value);
        if (rc) {
          throw std::runtime_error(fmt::format("Cannot read Stata's variable {}", j));
        }
        if (SF_is_missing(value)) {
          if (std::is_floating_point<T>::value) {
            value = MISSING_D;
          } else {
            value = 0;
          }
        }
        M[ind] = (T)value;
        ind += 1;
      }
      r += 1;
    }
  }

  return M;
}

/*
 * Write data to a column number 'j' in Stata (i.e. to a Stata 'variable').
 * If supplied, we consider each row 'i' only if filter[i] == true.
 */
void write_stata_column(ST_double* data, int len, ST_int j, const std::vector<bool>& filter = {})
{
  bool useEveryRow = (filter.size() == 0);
  int obs = 0;
  int r = 0; // Count each row that isn't filtered by Stata 'if'
  for (ST_int i = SF_in1(); i <= SF_in2(); i++) {
    if (SF_ifobs(i)) { // Skip rows according to Stata's 'if'
      if (useEveryRow || filter[r]) {
        // Convert MISSING_D back to Stata's missing value
        ST_double value = (data[obs] == MISSING_D) ? SV_missval : data[obs];
        ST_retcode rc = SF_vstore(j, i, value);
        if (rc) {
          throw std::runtime_error(fmt::format("Cannot write to Stata's variable {}", j));
        }
        obs += 1;
      }
      r += 1;
      if (obs >= len) {
        break;
      }
    }
  }
}

/*
 * Write data to columns ('variables') in Stata, starting from column number 'j0'.
 * If supplied, we consider each row 'i' only if filter[i] == true.
 */
void write_stata_columns(double* matrix, int matrixNumRows, int matrixNumCols, ST_int j0,
                         const std::vector<bool>& filter = {})
{
  bool useEveryRow = (filter.size() == 0);
  int obs = 0;
  int r = 0; // Count each row that isn't filtered by Stata 'if'
  for (ST_int i = SF_in1(); i <= SF_in2(); i++) {
    if (SF_ifobs(i)) { // Skip rows according to Stata's 'if'
      if (useEveryRow || filter[r]) {
        for (ST_int j = j0; j < j0 + matrixNumCols; j++) {
          // Convert MISSING_D back to Stata's missing value
          ST_double value = matrix[obs * matrixNumCols + (j - j0)];
          if (value == MISSING_D) {
            value = SV_missval;
          }
          ST_retcode rc = SF_vstore(j, i, value);
          if (rc) {
            throw std::runtime_error(fmt::format("Cannot write to Stata's variable {}", j));
          }
        }
        obs += 1;
      }
      r += 1;
      if (obs >= matrixNumRows) {
        break;
      }
    }
  }
}

std::vector<std::string> split_string(std::string list)
{
  std::vector<std::string> splitList;

  size_t found = list.find(' ');
  while (found != std::string::npos) {
    std::string part = list.substr(0, found);
    splitList.push_back(part);
    list = list.substr(found + 1);
    found = list.find(' ');
  }

  if (!list.empty()) {
    splitList.push_back(list);
  }

  return splitList;
}

template<typename T>
std::vector<T> numlist_to_vector(std::string list)
{
  std::vector<T> numList;

  size_t found = list.find(' ');
  while (found != std::string::npos) {
    std::string theta = list.substr(0, found);
    numList.push_back(atof(theta.c_str()));
    list = list.substr(found + 1);
    found = list.find(' ');
  }

  if (!list.empty()) {
    numList.push_back((T)atof(list.c_str()));
  }

  return numList;
}

template<typename T>
std::vector<T> stata_numlist(std::string macro)
{
  char buffer[1000];
  SF_macro_use((char*)("_" + macro).c_str(), buffer, 1000);

  std::string list(buffer);
  return numlist_to_vector<T>(list);
}

template<typename T>
void print_vector(std::string name, std::vector<T> vec)
{
  if (io.verbosity > 1) {
    io.print(fmt::format("{} [{}]:\n", name, vec.size()));
    for (int i = 0; i < vec.size(); i++) {
      if (i == 10) {
        io.print("... ");
        continue;
      }
      if (i > 10 && i < vec.size() - 10) {
        continue;
      }
      io.print(fmt::format("{} ", vec[i]));
    }
    io.print("\n");
  }
}

std::vector<int> bool_to_int(std::vector<bool> bv)
{
  std::vector<int> iv;
  std::copy(bv.begin(), bv.end(), std::back_inserter(iv));
  return iv;
}

std::vector<double> bool_to_double(std::vector<bool> bv)
{
  std::vector<double> dv;
  std::copy(bv.begin(), bv.end(), std::back_inserter(dv));
  return dv;
}

std::vector<std::future<PredictionResult>> futures;

// In case we have some remnants of previous runs still
// in the system (e.g. after a 'break'), clear our past results.
void reset_global_state()
{
  io.get_and_clear_async_buffer();
  breakButtonPressed = false;
  allTasksFinished = false;
}

ST_retcode launch_edm_tasks(int argc, char* argv[])
{
  if (argc < 28) {
    return TOO_FEW_VARIABLES;
  }
  if (argc > 28) {
    return TOO_MANY_VARIABLES;
  }

  reset_global_state();

  Options opts;
  opts.copredict = false;

  opts.calcRhoMAE = true;
  int numExtras = atoi(argv[0]);
  bool dtMode = atoi(argv[1]);
  opts.dtWeight = dtMode ? atof(argv[2]) : 0.0;
  std::string alg = std::string(argv[3]);
  if (alg.empty() || alg == "simplex") {
    opts.algorithm = Algorithm::Simplex;
  } else if (alg == "smap") {
    opts.algorithm = Algorithm::SMap;
  } else if (alg == "llr") {
    return NOT_IMPLEMENTED;
  } else {
    return INVALID_ALGORITHM;
  }
  opts.forceCompute = (std::string(argv[4]) == "force");
  opts.missingdistance = atof(argv[5]);
  char* reqThreads = argv[6];
  opts.nthreads = atoi(reqThreads);
  io.verbosity = atoi(argv[7]);
  opts.numTasks = atoi(argv[8]);
  bool explore = atoi(argv[9]);
  bool full = atoi(argv[10]);
  bool shuffle = atoi(argv[11]);
  int crossfold = atoi(argv[12]);
  int tau = atoi(argv[13]);
  int maxE = atoi(argv[14]);
  bool allowMissing = atoi(argv[15]);
  opts.thetas = numlist_to_vector<double>(std::string(argv[16]));
  opts.aspectRatio = atof(argv[17]);
  std::string distance(argv[18]);
  std::string requestedMetrics(argv[19]);
  bool copredictMode = atoi(argv[20]);
  opts.cmdLine = argv[21];
  int numExtrasLagged = atoi(argv[22]);
  opts.idw = atof(argv[23]);
  opts.panelMode = atoi(argv[24]);
  bool reldt = atoi(argv[25]);
  bool wassDT = atoi(argv[26]);
  int p = atoi(argv[27]);

  auto extrasFactorVariables = stata_numlist<bool>("z_factor_var");

  if (distance == "l1" || distance == "L1" || distance == "mae" || distance == "MAE") {
    opts.distance = Distance::MeanAbsoluteError;
  } else if (distance == "l2" || distance == "L2" || distance == "euclidean" || distance == "Euclidean") {
    opts.distance = Distance::Euclidean;
  } else if (distance == "wasserstein" || distance == "Wasserstein") {
    opts.distance = Distance::Wasserstein;
  } else {
    return INVALID_DISTANCE;
  }

  // Default number of threads is the number of physical cores available
  ST_int npcores = (ST_int)num_physical_cores();
  if (opts.nthreads <= 0) {
    opts.nthreads = npcores;
  }

  // Restrict going over the number of logical cores available
  ST_int nlcores = (ST_int)num_logical_cores();
  if (opts.nthreads > nlcores) {
    io.print(fmt::format("Restricting to {} threads (recommend {} threads)\n", nlcores, npcores));
    opts.nthreads = nlcores;
  }

  // Read in some macros from Stata
  char buffer[200];

  // What is k?
  if (SF_macro_use(NUM_NEIGHBOURS, buffer, 200)) {
    io.print("Got an error rc from macro_use!\n");
  }
  int k = atoi(buffer);

  // Are we saving the manifold?
  if (SF_macro_use(SAVE_MANIFOLD, buffer, 200)) {
    io.print("Got an error rc from macro_use!\n");
  }
  bool saveManifold = !(std::string(buffer).empty());

  // Are we saving the dt variable?
  if (SF_macro_use(SAVE_DT, buffer, 200)) {
    io.print("Got an error rc from macro_use!\n");
  }
  bool saveDT = !(std::string(buffer).empty());

  // Are we saving the predictions?
  if (SF_macro_use(SAVE_PREDICTION, buffer, 200)) {
    io.print("Got an error rc from macro_use!\n");
  }
  bool saveFinalPredictions = !(std::string(buffer).empty());

  // Are we saving the copredictions?
  if (SF_macro_use(SAVE_COPREDICTION, buffer, 200)) {
    io.print("Got an error rc from macro_use!\n");
  }
  bool saveFinalCoPredictions = !(std::string(buffer).empty());

  // Are we saving the S-map coefficients (only in xmap mode)?
  bool saveSMAPCoeffs;
  if (explore) {
    saveSMAPCoeffs = false;
  } else {
    if (SF_macro_use(SAVE_SMAP, buffer, 200)) {
      io.print("Got an error rc from macro_use!\n");
    }
    saveSMAPCoeffs = !(std::string(buffer).empty());
  }

  // Are we saving the inputs to a JSON file?
  if (SF_macro_use(SAVE_INPUTS, buffer, 200)) {
    io.print("Got an error rc from macro_use!\n");
  }
  std::string saveInputsFilename(buffer);

  // Number of replications
  if (SF_macro_use(NUM_REPS, buffer, 200)) {
    io.print("Got an error rc from macro_use!\n");
  }
  int numReps = atoi(buffer);

  std::vector<int> Es = stata_numlist<int>("e");

  // Read in the main data from Stata
  std::vector<ST_double> x = stata_columns<ST_double>(2);

  // Read in the time series which forms the basis of the y targets.
  // In 'explore' mode, it is just 'x' again so nothing is needed.
  std::vector<ST_double> xmap;
  if (!explore) {
    xmap = stata_columns<ST_double>(3);
  }

  // Read in the extras
  std::vector<std::vector<ST_double>> extras(numExtras);

  for (int z = 0; z < numExtras; z++) {
    extras[z] = stata_columns<ST_double>(4 + z);
  }

  std::vector<Metric> metrics;
  if (requestedMetrics == "auto" || requestedMetrics.empty()) {
    for (bool isFactorVariable : extrasFactorVariables) {
      metrics.push_back(isFactorVariable ? Metric::CheckSame : Metric::Diff);
    }
  } else {
    for (std::string& metric : split_string(requestedMetrics)) {
      if (metric == "same" || metric == "indicator" || metric == "onehot") {
        metrics.push_back(Metric::CheckSame);
      } else {
        metrics.push_back(Metric::Diff);
      }
    }

    // If the user supplied fewer than the required number of metrics,
    // just repeat the last one to pad out the list.
    while (metrics.size() < numExtras) {
      metrics.push_back(metrics.back());
    }
  }
  opts.metrics = metrics;

  std::vector<ST_double> t = stata_columns<ST_double>(1);
  print_vector<ST_double>("t", t);

  std::vector<ST_double> co_x;
  if (copredictMode) {
    co_x = stata_columns<ST_double>(3 + numExtras + 1);
  }

  std::vector<int> panelIDs;
  if (opts.panelMode) {
    panelIDs = stata_columns<int>(3 + numExtras + copredictMode + 1);
  }

  if (dtMode || (opts.distance == Distance::Wasserstein && wassDT)) {
    if (wassDT && !dtMode) {
      opts.dtWeight = 1.0;
      reldt = true;
    }
  }

  if (allowMissing && opts.missingdistance == 0) {
    opts.missingdistance = default_missing_distance(x);
  }
  SF_macro_save(MISSING_DISTANCE_USED, (char*)fmt::format("{}", opts.missingdistance).c_str());

  if (dtMode && opts.dtWeight == 0.0) {
    // If we have to set the default 'dt' weight, then make a manifold with dtweight of 1 then
    // we can rescale this by the appropriate variances in the future.
    // TODO: How would we change this for 'reldt'?
    const ManifoldGenerator dtgenerator(t, x, tau, p, xmap, co_x, panelIDs, extras, numExtrasLagged, dtMode, false,
                                        allowMissing);
    double DT_WEIGHT = 1.0;
    Manifold manifold = dtgenerator.create_manifold(maxE, {}, true, DT_WEIGHT);
    std::vector<double> dts(manifold.numPoints());
    for (int i = 0; i < dts.size(); i++) {
      dts[i] = manifold.dt(i, 1);
    }

    opts.dtWeight = default_dt_weight(dts, x);
    if (opts.dtWeight < 0) {
      dtMode = false;
    }
  }
  SF_macro_save(DTW_USED, (char*)fmt::format("{}", opts.dtWeight).c_str());

  const ManifoldGenerator generator(t, x, tau, p, xmap, co_x, panelIDs, extras, numExtrasLagged, dtMode, reldt,
                                    allowMissing);

  // Save some variables back to Stata, like the manifold (if savemanifold) or the dt (if dtsave).
  // Start by generating the 'usable' variable.
  std::vector<bool> usable = generator.generate_usable(maxE);

  int numUsable = std::accumulate(usable.begin(), usable.end(), 0);
  SF_macro_save(NUM_USABLE, (char*)fmt::format("{}", numUsable).c_str());

  // Save the dt column (before scaling by 'dtWeight') back to Stata if requested.
  if (saveDT) {
    double SAVE_DT_WEIGHT = 1.0;
    Manifold manifold = generator.create_manifold(maxE, {}, true, SAVE_DT_WEIGHT);

    std::vector<double> dts(manifold.numPoints());
    for (int i = 0; i < dts.size(); i++) {
      dts[i] = manifold.dt(i, 1);
    }
    ST_int startCol = 3 + numExtras + copredictMode + opts.panelMode + 1;
    write_stata_column(dts.data(), (int)dts.size(), startCol);
  }

  // Save the manifold back to Stata if requested.
  if (saveManifold) {
    Manifold manifold = generator.create_manifold(maxE, {}, true, opts.dtWeight);

    ST_int startCol = 3 + numExtras + copredictMode + opts.panelMode + saveDT + 1;
    write_stata_columns(manifold.data(), manifold.numPoints(), manifold.E_actual(), startCol);
  }

  if (numUsable == 0) {
    SF_scal_save(FINISHED_SCALAR, 1.0);
    return SUCCESS; // Let Stata give the error here.
  }

  // If we need to create a randomised library/prediction split, then sync the state of the
  // Mersenne Twister in Stata to that in the splitter instance.
  bool requiresRandomNumbers = numReps > 0 || (shuffle && !full);

  std::string rngState;

  if (requiresRandomNumbers) {
    char buffer[5200]; // Need at least 5011 + 1 bytes.
    if (SF_macro_use(RNG_STATE, buffer, 5200)) {
      io.print("Got an error rc from macro_use!\n");
    }

    rngState = std::string(buffer);

    // If the rng state is too small, it must be malformed.
    if (rngState.size() < 3 + 312 * 16 + 4 + 8) {
      rngState.clear();
    } else {
      // However, even if the right length, it may be all zeros, which indicates
      // an invalid state which should also be ignored.
      bool allZeros = true;
      for (int i = 0; i < 312; i++) {
        if (std::stoull(rngState.substr(3 + i * 16, 16), nullptr, 16) != 0) {
          allZeros = false;
          break;
        }
      }

      if (allZeros) {
        rngState.clear();
      }
    }
  }

  std::vector<int> libraries;
  if (!explore) {
    libraries = stata_numlist<int>("l_ori"); // The 'library' macro gets overwritten
    if (libraries.empty()) {
      libraries.push_back(numUsable);
    }
  }

  // If requested, save the inputs to a local file for testing
  if (!saveInputsFilename.empty()) {
    if (io.verbosity > 1) {
      io.print(fmt::format("Saving inputs to '{}.json'\n", saveInputsFilename));
      io.flush();
    }

    json taskGroup;
    taskGroup["generator"] = generator;
    taskGroup["opts"] = opts;
    taskGroup["Es"] = Es;
    taskGroup["libraries"] = libraries;
    taskGroup["k"] = k;
    taskGroup["numReps"] = numReps;
    taskGroup["crossfold"] = crossfold;
    taskGroup["explore"] = explore;
    taskGroup["full"] = full;
    taskGroup["shuffle"] = shuffle;
    taskGroup["saveFinalPredictions"] = saveFinalPredictions;
    taskGroup["saveFinalCoPredictions"] = saveFinalCoPredictions;
    taskGroup["saveSMAPCoeffs"] = saveSMAPCoeffs;
    taskGroup["copredictMode"] = copredictMode;
    taskGroup["usable"] = bool_to_int(usable);
    taskGroup["rngState"] = rngState;

    append_to_dumpfile(saveInputsFilename + ".json", taskGroup);

    // If we just want to save the input file and not actually run the command,
    // then uncomment the following two lines to end early.
    // SF_scal_save(FINISHED_SCALAR, 1.0);
    // return SUCCESS; // Let Stata give the error here.
  }

  futures = launch_task_group(generator, opts, Es, libraries, k, numReps, crossfold, explore, full, shuffle,
                              saveFinalPredictions, saveFinalCoPredictions, saveSMAPCoeffs, copredictMode, usable,
                              rngState, &io, keep_going, all_tasks_finished);

  return SUCCESS;
}

void save_to_stata_matrix(std::string matName, ST_int r, ST_int c, double v, ST_retcode& rc)
{
  if (v == MISSING_D) {
    v = SV_missval;
  }

  if (SF_mat_store((char*)matName.c_str(), r, c, v)) {
    io.error(fmt::format("Error: failed to save {} to matrix '{}[{},{}]'\n", v, matName, r, c).c_str());
    rc = CANNOT_SAVE_RESULTS;
  }
}

ST_retcode save_all_task_results_to_stata(int argc, char* argv[])
{
  if (argc < 3) {
    return TOO_FEW_VARIABLES;
  }
  if (argc > 3) {
    return TOO_MANY_VARIABLES;
  }

  std::string resultMatrix = argv[0];
  bool savePredictMode = atoi(argv[1]);
  bool saveCoPredictMode = atoi(argv[2]);

  ST_retcode rc = 0;
  int kMin, kMax;

  int numCoeffColsSaved = 0;

  for (int i = 0; i < futures.size(); i++) {

    // If there are no errors, store the predictions and S-map coefficients to Stata variables.
    const PredictionResult pred = futures[i].get();

    if (i == 0 || pred.kMin < kMin) {
      kMin = pred.kMin;
    }
    if (i == 0 || pred.kMax > kMax) {
      kMax = pred.kMax;
    }

    if (pred.rc == SUCCESS) {
      // Save the rho/MAE results
      std::string matName = pred.copredict ? "co_" + resultMatrix : resultMatrix;

      for (int t = 0; t < pred.stats.size(); t++) {
        ST_int matRow = pred.configNum + t + 1;
        if (pred.explore) {
          save_to_stata_matrix(matName, matRow, 1, pred.stats[t].E, rc);
          save_to_stata_matrix(matName, matRow, 2, pred.stats[t].theta, rc);
        } else {
          save_to_stata_matrix(matName, matRow, 2, pred.stats[t].library, rc);
        }
        save_to_stata_matrix(matName, matRow, 3, pred.stats[t].rho, rc);
        save_to_stata_matrix(matName, matRow, 4, pred.stats[t].mae, rc);
      }

      if (pred.predictions != nullptr) {
        if (savePredictMode && !pred.copredict) {
          write_stata_column(pred.predictions.get(), pred.numPredictions, 1, pred.predictionRows);
        }
        if (saveCoPredictMode && pred.copredict) {
          write_stata_column(pred.predictions.get(), pred.numPredictions, savePredictMode + 1, pred.predictionRows);
        }
      }

      if (pred.coeffs != nullptr) {
        write_stata_columns(pred.coeffs.get(), pred.numPredictions, pred.numCoeffCols,
                            savePredictMode + saveCoPredictMode + numCoeffColsSaved + 1, pred.predictionRows);
        numCoeffColsSaved += pred.numCoeffCols;
      }
    }

    if (pred.rc > rc) {
      rc = pred.rc;
    }
  }

  SF_macro_save(K_MIN, (char*)fmt::format("{}", kMin).c_str());
  SF_macro_save(K_MAX, (char*)fmt::format("{}", kMax).c_str());

  return rc;
}

STDLL stata_call(int argc, char* argv[])
{
  try {
    ST_retcode rc = UNKNOWN_ERROR + 1;
    std::string command(argv[0]);

    if (command == "launch_edm_tasks") {
      rc = launch_edm_tasks(argc - 1, argv + 1);
    } else if (command == "report_progress") {
      if (!breakButtonPressed) {
        io.print(io.get_and_clear_async_buffer());
      }

      bool breakHit = (argc == 2) && atoi(argv[1]);
      if (breakHit) {
        breakButtonPressed = true;
        io.out("Aborting edm run (this may take a few seconds).\n");
      }

      if (allTasksFinished) {
        SF_scal_save(FINISHED_SCALAR, 1.0);
      }

      rc = SUCCESS;
    } else if (command == "collect_results") {
      io.print(io.get_and_clear_async_buffer());
      rc = save_all_task_results_to_stata(argc - 1, argv + 1);
    } else {
      rc = UNKNOWN_ERROR;
    }

    print_error(command, rc);
    return rc;
  } catch (const std::exception& e) {
    io.error(e.what());
    io.error("\n");
  } catch (...) {
    io.error("Unknown error in edm plugin\n");
  }
  return UNKNOWN_ERROR;
}
