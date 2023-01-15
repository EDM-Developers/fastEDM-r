#pragma once

#define SUCCESS 0
#define BREAK_HIT 1
#define TOO_FEW_VARIABLES 102
#define TOO_MANY_VARIABLES 103
#define INVALID_ALGORITHM 400
#define INVALID_DISTANCE 401
#define INVALID_METRICS 402
#define INSUFFICIENT_UNIQUE 503
#define NOT_IMPLEMENTED 908
#define CANNOT_SAVE_RESULTS 1000
#define UNKNOWN_ERROR 8000

#define WITH_GPU_PROFILING 0

typedef int retcode;

#include <future>
#include <map>
#include <memory> // For unique_ptr
#include <queue>
#include <string>
#include <vector>

#include "manifold.h"
#if defined(WITH_ARRAYFIRE)
#include <arrayfire.h>
#endif

enum class Algorithm
{
  Simplex,
  SMap,
};

enum class Distance
{
  MeanAbsoluteError,
  Euclidean,
  Wasserstein
};

enum class Metric
{
  Diff,
  CheckSame
};

struct DistanceIndexPairs
{
  std::vector<int> inds;
  std::vector<double> dists;
};

#if defined(WITH_ARRAYFIRE)
struct DistanceIndexPairsOnGPU
{
  af::array valids;
  af::array dists;
};
#endif

struct Options
{
  bool explore;
  bool copredict;
  bool forceCompute;
  bool saveTargets;
  bool savePrediction;
  bool saveSMAPCoeffs;
  int k, nthreads, library;
  double missingdistance;
  bool panelMode;
  double idw;
  std::map<std::pair<int, int>, float> idWeights;
  std::vector<double> thetas;
  Algorithm algorithm;
  int taskNum, numTasks, configNum;
  bool calcRhoMAE;
  double aspectRatio;
  Distance distance;
  std::vector<Metric> metrics;
  std::string cmdLine;
  bool saveKUsed;
  bool saveManifolds;
  bool lowMemoryMode;
  bool useOnlyPastToPredictFuture;
  bool hasMissing;
  bool hasCategorical;
};

struct PredictionStats
{
  int library, E;
  double theta, mae, rho;
};

struct PredictionResult
{
  retcode rc;
  size_t numThetas, numPredictions, numCoeffCols;
  std::unique_ptr<double[]> targets;
  std::unique_ptr<double[]> predictions;
  std::unique_ptr<double[]> coeffs;
  std::unique_ptr<Manifold> M, Mp;
  std::vector<PredictionStats> stats;
  std::vector<bool> predictionRows;
  int kMin, kMax;
  std::string cmdLine;
  bool explore, copredict;
  int configNum;
};

class IO
{
public:
  int verbosity = 0;

  virtual void print(std::string s)
  {
    if (verbosity > 0) {
      out(s.c_str());
      flush();
    }
  }

  virtual void print_async(std::string s)
  {
    if (verbosity > 0) {
      std::lock_guard<std::mutex> guard(bufferMutex);
      buffer += s;
    }
  }

  virtual std::string get_and_clear_async_buffer()
  {
    std::lock_guard<std::mutex> guard(bufferMutex);
    std::string ret = buffer;
    buffer.clear();
    return ret;
  }

  virtual void progress_bar(double progress)
  {
    std::lock_guard<std::mutex> guard(bufferMutex);

    if (progress == 0.0) {
      finished = false;
      buffer += "Percent complete: 0";
      nextMessage = 1.0 / 40;
      dots = 0;
      tens = 0;
      return;
    }

    while (progress >= nextMessage && nextMessage < 1.0) {
      if (dots < 3) {
        buffer += ".";
        dots += 1;
      } else {
        tens += 1;
        buffer += std::to_string(tens * 10);
        dots = 0;
      }
      nextMessage += 1.0 / 40;
    }

    if (progress >= 1.0 && !finished) {
      buffer += "\n";
      finished = true;
    }
  }

  // Actual implementation of IO functions are in the subclasses
  virtual void out(const char*) const = 0;
  virtual void error(const char*) const = 0;
  virtual void flush() const = 0;

private:
  std::string buffer = "";
  std::mutex bufferMutex;

  bool finished;
  int dots, tens;
  double nextMessage;
};

#ifdef JSON

void to_json(json& j, const Algorithm& a);
void from_json(const json& j, Algorithm& a);

void to_json(json& j, const Distance& d);
void from_json(const json& j, Distance& d);

void to_json(json& j, const Metric& m);
void from_json(const json& j, Metric& m);

void to_json(json& j, const Options& o);
void from_json(const json& j, Options& o);

void to_json(json& j, const PredictionStats& s);
void from_json(const json& j, PredictionStats& s);

void to_json(json& j, const PredictionResult& p);
void from_json(const json& j, PredictionResult& p);

#endif
