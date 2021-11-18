#include <benchmark/benchmark.h>

#ifdef _WIN32
#pragma comment(lib, "Shlwapi.lib")
#ifdef _DEBUG
#pragma comment(lib, "benchmarkd.lib")
#else
#pragma comment(lib, "benchmark.lib")
#endif
#endif

#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
#endif
#include <fmt/format.h>

#include "cli.h"
#include "common.h"
#include "cpu.h"
#include "distances.h"
#include "edm.h"

#define EIGEN_NO_DEBUG
#define EIGEN_DONT_PARALLELIZE
#include <Eigen/SVD>

// Function declarations for 'private' functions not listed in the relevant header files.
std::vector<int> potential_neighbour_indices(int Mp_i, const Options& opts, const Manifold& M, const Manifold& Mp);
DistanceIndexPairs kNearestNeighbours(const DistanceIndexPairs& potentialNeighbours, int k);

// Compiler flags tried on Windows: "/GL" and "/GL /LTCG", both slightly worse. "/O2" is the default.

std::vector<std::string> lowLevelInputDumps = { "logmapsmall.json", "logmaplarge.json", "affectsmall.json",
                                                "affectbige.json" };

ConsoleIO io(0);

static void bm_basic_distances(benchmark::State& state)
{
  std::string filename = lowLevelInputDumps[state.range(0)];
  state.SetLabel(filename);

  Inputs vars = read_lowlevel_inputs_file(filename);
  Manifold M = vars.generator.create_manifold(vars.E, vars.libraryRows, false, false);
  Manifold Mp = vars.generator.create_manifold(vars.E, vars.predictionRows, false, true);

  int Mp_i = 0;

  for (auto _ : state) {
    std::vector<int> tryInds = potential_neighbour_indices(Mp_i, vars.opts, M, Mp);
    DistanceIndexPairs potentialNN = lp_distances(Mp_i, vars.opts, M, Mp, tryInds);
    Mp_i = (Mp_i + 1) % Mp.numPoints();
  }
}

BENCHMARK(bm_basic_distances)->DenseRange(0, lowLevelInputDumps.size() - 1)->Unit(benchmark::kMicrosecond);

static void bm_wasserstein_distances(benchmark::State& state)
{
  std::string filename = lowLevelInputDumps[state.range(0)];
  state.SetLabel(filename);

  Inputs vars = read_lowlevel_inputs_file(filename);
  Manifold M = vars.generator.create_manifold(vars.E, vars.libraryRows, false, false);
  Manifold Mp = vars.generator.create_manifold(vars.E, vars.predictionRows, false, true);

  int Mp_i = 0;

  for (auto _ : state) {
    std::vector<int> tryInds = potential_neighbour_indices(Mp_i, vars.opts, M, Mp);
    DistanceIndexPairs potentialNN = wasserstein_distances(Mp_i, vars.opts, M, Mp, tryInds);
    Mp_i = (Mp_i + 1) % Mp.numPoints();
  }
}

BENCHMARK(bm_wasserstein_distances)->DenseRange(0, lowLevelInputDumps.size() - 1)->Unit(benchmark::kMillisecond);

static void bm_nearest_neighbours(benchmark::State& state)
{
  std::string filename = lowLevelInputDumps[state.range(0)];
  state.SetLabel(filename);

  Inputs vars = read_lowlevel_inputs_file(filename);
  Manifold M = vars.generator.create_manifold(vars.E, vars.libraryRows, false, false);
  Manifold Mp = vars.generator.create_manifold(vars.E, vars.predictionRows, false, true);

  int Mp_i = 0;
  Options opts = vars.opts;

  std::vector<int> tryInds = potential_neighbour_indices(Mp_i, vars.opts, M, Mp);
  DistanceIndexPairs potentialNN;

  if (opts.distance == Distance::Wasserstein) {
    potentialNN = wasserstein_distances(Mp_i, vars.opts, M, Mp, tryInds);
  } else {
    potentialNN = lp_distances(Mp_i, vars.opts, M, Mp, tryInds);
  }

  int numValidDistances = potentialNN.inds.size();
  int k = opts.k;
  if (k > numValidDistances) {
    k = numValidDistances;
  }

  // If k < 0 then we normally skip this step, so let's set it one less than that.
  if (k < 0 || k == numValidDistances) {
    k = numValidDistances - 1;
  }

  for (auto _ : state) {
    DistanceIndexPairs kNNs = kNearestNeighbours(potentialNN, k);
  }
}

BENCHMARK(bm_nearest_neighbours)->DenseRange(0, lowLevelInputDumps.size() - 1)->Unit(benchmark::kMicrosecond);

static void bm_simplex(benchmark::State& state)
{
  std::string filename = lowLevelInputDumps[state.range(0)];
  state.SetLabel(filename);

  Inputs vars = read_lowlevel_inputs_file(filename);
  Manifold M = vars.generator.create_manifold(vars.E, vars.libraryRows, false, false);
  Manifold Mp = vars.generator.create_manifold(vars.E, vars.predictionRows, false, true);

  int Mp_i = 0;
  Options opts = vars.opts;

  std::vector<int> tryInds = potential_neighbour_indices(Mp_i, vars.opts, M, Mp);
  DistanceIndexPairs potentialNN;

  if (opts.distance == Distance::Wasserstein) {
    potentialNN = wasserstein_distances(Mp_i, vars.opts, M, Mp, tryInds);
  } else {
    potentialNN = lp_distances(Mp_i, vars.opts, M, Mp, tryInds);
  }

  int numValidDistances = potentialNN.inds.size();
  int k = opts.k;
  if (k > numValidDistances) {
    k = numValidDistances;
  }

  // If we asked for all of the neighbours to be considered (e.g. with k = -1), return this index vector directly.
  DistanceIndexPairs kNNs;
  if (k < 0 || k == potentialNN.inds.size()) {
    kNNs = potentialNN;
  } else {
    kNNs = kNearestNeighbours(potentialNN, k);
  }

  auto dists = kNNs.dists;
  auto kNNInds = kNNs.inds;

  for (auto _ : state) {
    for (int t = 0; t < opts.thetas.size(); t++) {
      int k = kNNInds.size();

      // Find the smallest distance (closest neighbour) among the supplied neighbours.
      double minDist = *std::min_element(dists.begin(), dists.end());

      // Calculate our weighting of each neighbour, and the total sum of these weights.
      std::vector<double> w(k);
      double sumw = 0.0;
      const double theta = opts.thetas[t];

      for (int j = 0; j < k; j++) {
        w[j] = exp(-theta * (dists[j] / minDist));
        sumw = sumw + w[j];
      }

      // Make the simplex projection/prediction.
      double r = 0.0;
      for (int j = 0; j < k; j++) {
        r = r + M.target(kNNInds[j]) * (w[j] / sumw);
      }
    }
  }
}

BENCHMARK(bm_simplex)->DenseRange(0, lowLevelInputDumps.size() - 1);

static void bm_smap(benchmark::State& state)
{
  std::string filename = lowLevelInputDumps[state.range(0)];
  state.SetLabel(filename);

  Inputs vars = read_lowlevel_inputs_file(filename);
  Manifold M = vars.generator.create_manifold(vars.E, vars.libraryRows, false, false);
  Manifold Mp = vars.generator.create_manifold(vars.E, vars.predictionRows, false, true);

  int Mp_i = 0;
  Options opts = vars.opts;

  std::vector<int> tryInds = potential_neighbour_indices(Mp_i, vars.opts, M, Mp);
  DistanceIndexPairs potentialNN;

  if (opts.distance == Distance::Wasserstein) {
    potentialNN = wasserstein_distances(Mp_i, vars.opts, M, Mp, tryInds);
  } else {
    potentialNN = lp_distances(Mp_i, vars.opts, M, Mp, tryInds);
  }

  int numValidDistances = potentialNN.inds.size();
  int k = opts.k;
  if (k > numValidDistances) {
    k = numValidDistances;
  }

  // If we asked for all of the neighbours to be considered (e.g. with k = -1), return this index vector directly.
  DistanceIndexPairs kNNs;
  if (k < 0 || k == potentialNN.inds.size()) {
    kNNs = potentialNN;
  } else {
    kNNs = kNearestNeighbours(potentialNN, k);
  }

  auto dists = kNNs.dists;
  auto kNNInds = kNNs.inds;

  for (auto _ : state) {
    for (int t = 0; t < opts.thetas.size(); t++) {
      int k = kNNInds.size();

      // Pull out the nearest neighbours from the manifold, and
      // simultaneously prepend a column of ones in front of the manifold data.
      Eigen::MatrixXd X_ls_cj(k, M.E_actual() + 1);
      X_ls_cj << Eigen::VectorXd::Ones(k), M.map()(kNNInds, Eigen::all);

      // Calculate the weight for each neighbour
      Eigen::Map<const Eigen::VectorXd> distsMap(&(dists[0]), dists.size());
      Eigen::VectorXd w = Eigen::exp(-opts.thetas[t] * (distsMap.array() / distsMap.mean()));

      // Scale everything by our weights vector
      X_ls_cj.array().colwise() *= w.array();
      Eigen::VectorXd y_ls = M.targetsMap()(kNNInds).array() * w.array();

      // The old way to solve this system:
      // Eigen::BDCSVD<Eigen::MatrixXd> svd(X_ls_cj, Eigen::ComputeThinU | Eigen::ComputeThinV);
      //  Eigen::VectorXd ics = svd.solve(y_ls);

      // The pseudo-inverse of X can be calculated as (X^T * X)^(-1) * X^T
      // see https://scicomp.stackexchange.com/a/33375
      const int svdOpts =
        Eigen::ComputeThinU | Eigen::ComputeThinV; // 'ComputeFull*' would probably work identically here.
      Eigen::BDCSVD<Eigen::MatrixXd> svd(X_ls_cj.transpose() * X_ls_cj, svdOpts);
      Eigen::VectorXd ics = svd.solve(X_ls_cj.transpose() * y_ls);

      double r = ics(0);
      for (int j = 0; j < M.E_actual(); j++) {
        if (Mp(Mp_i, j) != MISSING_D) {
          r += Mp(Mp_i, j) * ics(j + 1);
        }
      }
    }
  }
}

BENCHMARK(bm_smap)->DenseRange(0, lowLevelInputDumps.size() - 1);

// Larger tests which indicate typical use-cases of the plugin
std::vector<std::string> tests = {
  "chicago.json",             // ~ 1 min on 80 threads
  "affect.json",              // ~ 10 mins on 80 threads
  "chicago-wasserstein.json", // ~ 20 mins on 80 threads
};

static void bm_run_tests(benchmark::State& state)
{
  std::string filename = tests[state.range(0)];
  state.SetLabel(filename);

  std::ifstream i(filename);
  json testInputs;
  i >> testInputs;

  bool verbose = false;
  int nthreads = 80;

  ConsoleIO io(verbose);

  for (auto _ : state) {
    auto results = run_tests(testInputs, nthreads, &io);
  }
}

BENCHMARK(bm_run_tests)->DenseRange(0, tests.size() - 1)->Unit(benchmark::kMillisecond);
