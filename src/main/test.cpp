#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one cpp file
#include <catch2/catch.hpp>

#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
#endif
#include <fmt/format.h>

#include "edm.h"
#include "library_prediction_split.h"
#include "manifold.h"

const double NA = MISSING_D;

// Function declarations for 'private' functions not listed in the relevant header files.
std::vector<int> potential_neighbour_indices(int Mp_i, const Options& opts, const Manifold& M, const Manifold& Mp);

std::unique_ptr<double[]> wasserstein_cost_matrix(const Manifold& M, const Manifold& Mp, int i, int j,
                                                  const Options& opts, int& len_i, int& len_j);

DistanceIndexPairs wasserstein_distances(int Mp_i, const Options& opts, const Manifold& M, const Manifold& Mp,
                                         std::vector<int> inds);

#if defined(WITH_ARRAYFIRE)
DistanceIndexPairs afWassersteinDistances(int Mp_i, const Options& opts, const Manifold& hostM, const Manifold& hostMp,
                                          const ManifoldOnGPU& M, const ManifoldOnGPU& Mp, const std::vector<int>& inds,
                                          const af::array& metricOpts);
#endif

void print_raw_matrix(const double* M, int rows, int cols)
{

  auto stringVersion = [](double v) { return (v == NA) ? std::string(" . ") : fmt::format("{:.1f}", v); };

  std::cout << "\n";
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      std::cout << stringVersion(M[i * cols + j]) << " (" << i * cols + j << ") "
                << "\t";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}

void print_eig_matrix(const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>& M)
{

  auto stringVersion = [](double v) { return (v == NA) ? std::string(" . ") : fmt::format("{:.1f}", v); };

  std::cout << "\n";
  for (int i = 0; i < M.rows(); i++) {
    for (int j = 0; j < M.cols(); j++) {
      std::cout << stringVersion(M(i, j)) << "\t";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}

void print_manifold(const Manifold& M)
{
  auto stringVersion = [](double v) { return (v == NA) ? std::string(" . ") : fmt::format("{:.1f}", v); };

  std::cout << "\n";
  for (int i = 0; i < M.numPoints(); i++) {
    for (int j = 0; j < M.E_actual(); j++) {
      std::cout << stringVersion(M(i, j)) << "\t";
    }
    std::cout << "\t|\t" << stringVersion(M.target(i)) << "\n";
  }
  std::cout << "\n";
}

template<typename T>
void require_vectors_match(const std::vector<T>& u, const std::vector<T>& v)
{
  REQUIRE(u.size() == v.size());

  for (int i = 0; i < u.size(); i++) {
    CAPTURE(i);
    REQUIRE(u[i] == v[i]);
  }
}

void require_manifolds_match(const Manifold& M, const std::vector<std::vector<double>>& M_true,
                             const std::vector<double>& y_true)
{
  REQUIRE(M.numPoints() == M_true.size());
  REQUIRE(M.numTargets() == y_true.size()); // Shouldn't this always be the same as M.numPoints()?
  REQUIRE(M.E_actual() == M_true[0].size());

  for (int i = 0; i < M.numPoints(); i++) {
    CAPTURE(i);
    for (int j = 0; j < M.E_actual(); j++) {
      CAPTURE(j);
      REQUIRE(M(i, j) == M_true[i][j]);
    }
    REQUIRE(M.target(i) == y_true[i]);
  }
}

void check_usable_matches_prediction_set(const std::vector<bool>& usable, const Manifold& Mp)
{
  int numUsable = std::accumulate(usable.begin(), usable.end(), 0);
  REQUIRE(numUsable == Mp.numPoints());
}

TEST_CASE("Basic manifold creation", "[basicManifold]")
{
  int E = 2;
  int tau = 1;
  int p = 1;

  std::vector<double> t = { 1, 2, 3, 4 };
  std::vector<double> x = { 11, 12, 13, 14 };

  SECTION("Basic manifold, no extras or dt")
  {
    ManifoldGenerator generator(t, x, tau, p);

    std::vector<bool> usable = generator.generate_usable(E);
    std::vector<bool> usableTrue = { false, true, true, true };
    require_vectors_match<bool>(usable, usableTrue);

    Manifold M = generator.create_manifold(E, usable, false);
    std::vector<std::vector<double>> M_true = { { 12.0, 11.0 }, { 13.0, 12.0 } };
    std::vector<double> y_true = { 13.0, 14.0 };
    require_manifolds_match(M, M_true, y_true);

    Manifold Mp = generator.create_manifold(E, usable, true);
    check_usable_matches_prediction_set(usable, Mp);
    std::vector<std::vector<double>> Mp_true = { { 12.0, 11.0 }, { 13.0, 12.0 }, { 14.0, 13.0 } };
    std::vector<double> yp_true = { 13.0, 14.0, NA };
    require_manifolds_match(Mp, Mp_true, yp_true);
  }

  SECTION("Manifold with dt (not allowing missing)")
  {
    // TODO: This test is a bit fake; edm.ado would not allow dt to be applied when there's no gaps in the time
    // variable.
    bool dtMode = true;
    double dtWeight = 1.0;
    bool allowMissing = false;
    ManifoldGenerator generator(t, x, tau, p, {}, {}, {}, {}, 0, dtMode, false, allowMissing);

    std::vector<bool> usable = generator.generate_usable(E);
    std::vector<bool> usableTrue = { false, true, true, false };
    require_vectors_match<bool>(usable, usableTrue);

    Manifold M = generator.create_manifold(E, usable, false, dtWeight);
    std::vector<std::vector<double>> M_true = { { 12.0, 11.0, 1.0, 1.0 }, { 13.0, 12.0, 1.0, 1.0 } };
    std::vector<double> y_true = { 13.0, 14.0 };
    require_manifolds_match(M, M_true, y_true);

    // Here, there's no difference between library and prediction sets
    Manifold Mp = generator.create_manifold(E, usable, true, dtWeight);
    check_usable_matches_prediction_set(usable, Mp);
    require_manifolds_match(Mp, M_true, y_true);
  }
}

// These tests are used as examples in the Julia docs
TEST_CASE("Missing data manifold creation (tau = 1)", "[missingDataManifold]")
{
  int E = 2;
  int tau = 1;
  int p = 1;

  std::vector<double> t = { 1.0, 2.5, 3.0, 4.5, 5.0, 6.0 };
  std::vector<double> x = { 11, 12, NA, 14, 15, 16 };

  SECTION("Default")
  {
    ManifoldGenerator generator(t, x, tau, p);

    REQUIRE(generator.calculate_time_increment() == 0.5);

    std::vector<int> obsNums = { 0, 3, 4, 7, 8, 10 };
    for (int i = 0; i < obsNums.size(); i++) {
      CAPTURE(i);
      REQUIRE(generator.get_observation_num(i) == obsNums[i]);
    }

    std::vector<bool> usable = generator.generate_usable(E);
    std::vector<bool> usableTrue = { false, false, false, false, true, false };
    require_vectors_match<bool>(usable, usableTrue);

    Manifold M = generator.create_manifold(E, usable, false);
    REQUIRE(M.numPoints() == 0);
    REQUIRE(M.numTargets() == 0);
    REQUIRE(M.E_actual() == 2);

    Manifold Mp = generator.create_manifold(E, usable, true);
    check_usable_matches_prediction_set(usable, Mp);
    std::vector<std::vector<double>> Mp_true = { { 15.0, 14.0 } };
    std::vector<double> yp_true = { NA };
    require_manifolds_match(Mp, Mp_true, yp_true);
  }

  SECTION("dt")
  {
    bool dtMode = true;
    double dtWeight = 1.0;
    bool allowMissing = false;
    ManifoldGenerator generator(t, x, tau, p, {}, {}, {}, {}, 0, dtMode, false, allowMissing);

    std::vector<int> obsNums = { 0, 1, -1, 2, 3, 4 };
    for (int i = 0; i < obsNums.size(); i++) {
      CAPTURE(i);
      REQUIRE(generator.get_observation_num(i) == obsNums[i]);
    }

    std::vector<bool> usable = generator.generate_usable(E);
    std::vector<bool> usableTrue = { false, true, false, true, true, false };
    require_vectors_match<bool>(usable, usableTrue);

    Manifold M = generator.create_manifold(E, usable, false, dtWeight);
    std::vector<std::vector<double>> M_true = { { 12.0, 11.0, 2.0, 1.5 },
                                                { 14.0, 12.0, 0.5, 2.0 },
                                                { 15.0, 14.0, 1.0, 0.5 } };
    std::vector<double> y_true = { 14.0, 15.0, 16.0 };
    require_manifolds_match(M, M_true, y_true);

    // Here, there's no difference between library and prediction sets
    Manifold Mp = generator.create_manifold(E, usable, true, dtWeight);
    check_usable_matches_prediction_set(usable, Mp);
    require_manifolds_match(Mp, M_true, y_true);
  }

  SECTION("dt and allowingmissing")
  {
    bool dtMode = true;
    double dtWeight = 1.0;
    bool allowMissing = true;
    ManifoldGenerator generator(t, x, tau, p, {}, {}, {}, {}, 0, dtMode, false, allowMissing);

    std::vector<int> obsNums = { 0, 1, 2, 3, 4, 5 };
    for (int i = 0; i < obsNums.size(); i++) {
      CAPTURE(i);
      REQUIRE(generator.get_observation_num(i) == obsNums[i]);
    }

    std::vector<bool> usable = generator.generate_usable(E);
    std::vector<bool> usableTrue = { true, true, true, true, true, true };
    require_vectors_match<bool>(usable, usableTrue);

    Manifold M = generator.create_manifold(E, usable, false, dtWeight);
    std::vector<std::vector<double>> M_true = {
      { 11.0, NA, 1.5, NA }, { NA, 12.0, 1.5, 0.5 }, { 14.0, NA, 0.5, 1.5 }, { 15.0, 14.0, 1.0, 0.5 }
    };
    std::vector<double> y_true = { 12.0, 14.0, 15.0, 16.0 };
    require_manifolds_match(M, M_true, y_true);

    Manifold Mp = generator.create_manifold(E, usable, true, dtWeight);
    check_usable_matches_prediction_set(usable, Mp);
    std::vector<std::vector<double>> Mp_true = { { 11.0, NA, 1.5, NA },    { 12.0, 11.0, 0.5, 1.5 },
                                                 { NA, 12.0, 1.5, 0.5 },   { 14.0, NA, 0.5, 1.5 },
                                                 { 15.0, 14.0, 1.0, 0.5 }, { 16.0, 15.0, NA, 1.0 } };
    std::vector<double> yp_true = { 12.0, NA, 14.0, 15.0, 16.0, NA };
    require_manifolds_match(Mp, Mp_true, yp_true);
  }

  SECTION("reldt")
  {
    bool dtMode = true, reldtMode = true;
    double dtWeight = 1.0;
    bool allowMissing = false;
    ManifoldGenerator generator(t, x, tau, p, {}, {}, {}, {}, 0, dtMode, reldtMode, allowMissing);

    std::vector<int> obsNums = { 0, 1, -1, 2, 3, 4 };
    for (int i = 0; i < obsNums.size(); i++) {
      CAPTURE(i);
      REQUIRE(generator.get_observation_num(i) == obsNums[i]);
    }

    std::vector<bool> usable = generator.generate_usable(E);
    std::vector<bool> usableTrue = { false, true, false, true, true, false };
    require_vectors_match<bool>(usable, usableTrue);

    Manifold M = generator.create_manifold(E, usable, false, dtWeight);
    std::vector<std::vector<double>> M_true = { { 12.0, 11.0, 2.0, 3.5 },
                                                { 14.0, 12.0, 0.5, 2.5 },
                                                { 15.0, 14.0, 1.0, 1.5 } };
    std::vector<double> y_true = { 14.0, 15.0, 16.0 };
    require_manifolds_match(M, M_true, y_true);

    // Here, there's no difference between library and prediction sets
    Manifold Mp = generator.create_manifold(E, usable, true, dtWeight);
    check_usable_matches_prediction_set(usable, Mp);
    require_manifolds_match(Mp, M_true, y_true);
  }
}

TEST_CASE("Missing data dt manifold creation (tau = 2)", "[missingDataManifold2]")
{
  int E = 2;
  int tau = 2;
  int p = 1;

  std::vector<double> t = { 1.0, 2.5, 3.0, 4.5, 5.0, 6.0 };
  std::vector<double> x = { 11, 12, NA, 14, 15, 16 };

  SECTION("Allowing missing values")
  {
    bool dtMode = true;
    double dtWeight = 1.0;
    bool allowMissing = true;
    ManifoldGenerator generator(t, x, tau, p, {}, {}, {}, {}, 0, dtMode, false, allowMissing);

    std::vector<bool> usable = generator.generate_usable(E);
    std::vector<bool> usableTrue = { true, true, true, true, true, true };
    require_vectors_match<bool>(usable, usableTrue);

    Manifold M = generator.create_manifold(E, usable, false, dtWeight);
    std::vector<std::vector<double>> M_true = {
      { 11.0, NA, 1.5, NA },
      { NA, 11.0, 1.5, 2.0 },
      { 14.0, 12.0, 0.5, 2.0 },
      { 15.0, NA, 1.0, 2.0 },
    };
    std::vector<double> y_true = { 12.0, 14.0, 15.0, 16.0 };
    require_manifolds_match(M, M_true, y_true);

    Manifold Mp = generator.create_manifold(E, usable, true, dtWeight);
    check_usable_matches_prediction_set(usable, Mp);
    std::vector<std::vector<double>> Mp_true = {
      { 11.0, NA, 1.5, NA },    { 12.0, NA, 0.5, NA },  { NA, 11.0, 1.5, 2.0 },
      { 14.0, 12.0, 0.5, 2.0 }, { 15.0, NA, 1.0, 2.0 }, { 16.0, 14.0, NA, 1.5 },
    };
    std::vector<double> yp_true = { 12.0, NA, 14.0, 15.0, 16.0, NA };
    require_manifolds_match(Mp, Mp_true, yp_true);
  }

  SECTION("Not allowing missing values")
  {
    bool dtMode = true;
    double dtWeight = 1.0;
    bool allowMissing = false;
    ManifoldGenerator generator(t, x, tau, p, {}, {}, {}, {}, 0, dtMode, false, allowMissing);

    std::vector<bool> usable = generator.generate_usable(E);
    std::vector<bool> usableTrue = { false, false, false, true, true, false };
    require_vectors_match<bool>(usable, usableTrue);

    Manifold M = generator.create_manifold(E, usable, false, dtWeight);
    std::vector<std::vector<double>> M_true = { { 14.0, 11.0, 0.5, 3.5 }, { 15.0, 12.0, 1.0, 2.5 } };
    std::vector<double> y_true = { 15.0, 16.0 };
    require_manifolds_match(M, M_true, y_true);

    // Here, there's no difference between library and prediction sets
    Manifold Mp = generator.create_manifold(E, usable, true, dtWeight);
    check_usable_matches_prediction_set(usable, Mp);
    require_manifolds_match(Mp, M_true, y_true);
  }
}

TEST_CASE("Check negative times work", "[negativeTimes]")
{
  std::vector<double> t = { -9.0, -7.5, -7.0, -6.5, -5.0, -4.0 };
  std::vector<double> x = { 11, 12, NA, 14, 15, 16 };

  int tau = 1;
  int p = 1;

  ManifoldGenerator generator(t, x, tau, p);

  REQUIRE(generator.calculate_time_increment() == 0.5);

  std::vector<int> obsNums = { 0, 3, 4, 5, 8, 10 };
  for (int i = 0; i < obsNums.size(); i++) {
    CAPTURE(i);
    REQUIRE(generator.get_observation_num(i) == obsNums[i]);
  }
}

TEST_CASE("Library prediction splits", "[splitting]")
{

  std::vector<bool> usable = { true, true, true, false, false, true, true };

  SECTION("Explore")
  {
    LibraryPredictionSetSplitter splitter(true, false, false, 0, usable);

    splitter.update_library_prediction_split();
    std::vector<Set> splitTrue = { Set::Library, Set::Library,    Set::Prediction, Set::Neither,
                                   Set::Neither, Set::Prediction, Set::Prediction };
    require_vectors_match<Set>(splitter.setMemberships(), splitTrue);
  }

  SECTION("Explore and full")
  {
    LibraryPredictionSetSplitter splitter(true, true, false, 0, usable);

    splitter.update_library_prediction_split();
    std::vector<Set> splitTrue = { Set::Both, Set::Both, Set::Both, Set::Neither, Set::Neither, Set::Both, Set::Both };
    require_vectors_match<Set>(splitter.setMemberships(), splitTrue);
  }

  SECTION("Crossfold mode (explore)")
  {
    int crossfold = 3;

    LibraryPredictionSetSplitter splitter(true, false, false, crossfold, usable);

    std::vector<std::vector<Set>> splitsTrue = {
      { Set::Prediction, Set::Prediction, Set::Library, Set::Neither, Set::Neither, Set::Library, Set::Library },
      { Set::Library, Set::Library, Set::Prediction, Set::Neither, Set::Neither, Set::Prediction, Set::Library },
      { Set::Library, Set::Library, Set::Library, Set::Neither, Set::Neither, Set::Library, Set::Prediction },
    };

    for (int iter = 1; iter <= crossfold; iter++) {
      splitter.update_library_prediction_split(-1, iter);
      require_vectors_match<Set>(splitter.setMemberships(), splitsTrue[iter - 1]);
    }
  }

  SECTION("Xmap")
  {
    LibraryPredictionSetSplitter splitter(false, false, false, 0, usable);

    splitter.update_library_prediction_split(1);
    std::vector<Set> splitTrue1 = { Set::Both,    Set::Prediction, Set::Prediction, Set::Neither,
                                    Set::Neither, Set::Prediction, Set::Prediction };
    require_vectors_match<Set>(splitter.setMemberships(), splitTrue1);

    splitter.update_library_prediction_split(4);
    std::vector<Set> splitTrue4 = { Set::Both,    Set::Both, Set::Both,      Set::Neither,
                                    Set::Neither, Set::Both, Set::Prediction };
    require_vectors_match<Set>(splitter.setMemberships(), splitTrue4);
  }
}

TEST_CASE("Wasserstein distance", "[wasserstein]")
{
  int E = 5;
  int tau = 1;
  int p = 0;

  std::vector<double> t = { 0, 1, 2, 3, 4 };
  std::vector<double> x = { 1, 2, NA, NA, 5 };
  std::vector<double> co_x = { 1, NA, NA, 4, 5 };
  auto xmap = co_x;

  bool dt = true, reldt = true, allowMissing = true;
  ManifoldGenerator generator(t, x, tau, p, xmap, co_x, {}, {}, 0, dt, reldt, allowMissing);

  std::vector<bool> usable = generator.generate_usable(E);
  std::vector<bool> usableTrue = { true, true, true, true, true };
  require_vectors_match<bool>(usable, usableTrue);

  double dtWeight = 1.0;
  bool copredict = true;
  Manifold M = generator.create_manifold(E, usable, false, dtWeight, copredict);
  Manifold Mp = generator.create_manifold(E, usable, true, dtWeight, copredict);
  check_usable_matches_prediction_set(usable, Mp);

  std::vector<std::vector<double>> M_true = {
    { 1, NA, NA, NA, NA, 0, NA, NA, NA, NA },
    { NA, NA, 2, 1, NA, 0, 1, 2, 3, NA },
    { 5, NA, NA, 2, 1, 0, 1, 2, 3, 4 },
  };
  std::vector<double> y_true = { 1, 4, 5 };
  require_manifolds_match(M, M_true, y_true);

  std::vector<std::vector<double>> Mp_true = {
    { 1, NA, NA, NA, NA, 0, NA, NA, NA, NA }, { NA, 1, NA, NA, NA, 0, 1, NA, NA, NA },
    { NA, NA, 1, NA, NA, 0, 1, 2, NA, NA },   { 4, NA, NA, 1, NA, 0, 1, 2, 3, NA },
    { 5, 4, NA, NA, 1, 0, 1, 2, 3, 4 },
  };
  std::vector<double> yp_true = { 1, NA, NA, 4, 5 };
  require_manifolds_match(Mp, Mp_true, yp_true);

  Options opts;

  opts.missingdistance = 0;
  opts.aspectRatio = 1.0;
  opts.panelMode = false;

  opts.metrics = {};
  for (int i = 0; i < M.E_actual(); i++) {
    opts.metrics.push_back(Metric::Diff);
  }

  SECTION("Cost matrix")
  {
    int i = 2, j = 4;

    auto M_i = M.laggedObsMap(i);
    auto Mp_j = Mp.laggedObsMap(j);

    REQUIRE(M_i.rows() == 2);
    REQUIRE(Mp_j.rows() == 2);

    REQUIRE(M_i.cols() == 5);
    REQUIRE(Mp_j.cols() == 5);

    int len_i, len_j;
    std::unique_ptr<double[]> C = wasserstein_cost_matrix(M, Mp, i, j, opts, len_i, len_j);

    REQUIRE(len_i == 3);
    REQUIRE(len_j == 3);

    std::vector<double> C_true = { 0, 2, 8, 6, 4, 2, 8, 6, 0 };

    for (int ind = 0; ind < C_true.size(); ind++) {
      REQUIRE(C[ind] == C_true[ind]);
    }
  }

  SECTION("Multiple Wasserstein distances")
  {

    int Mp_j = 2;

    std::vector<int> tryInds = potential_neighbour_indices(Mp_j, opts, M, Mp);
    for (int i = 0; i < M.numPoints(); i++) {
      // std::cout << fmt::format("potential_neighbour_indices[{}] = {}\n", i, tryInds[i]);
    }

    for (int i = 0; i < M.numPoints(); i++) {
      // std::cout << fmt::format("Cost matrix M_i={}\n", i);

      auto M_i_map = M.laggedObsMap(i);
      auto Mp_j_map = Mp.laggedObsMap(Mp_j);

      // std::cout << "M_i_map={}\n";
      // print_eig_matrix(M_i_map);
      //
      // std::cout << "Mp_j_map={}\n";
      // print_eig_matrix(Mp_j_map);

      int len_i, len_j;
      std::unique_ptr<double[]> C = wasserstein_cost_matrix(M, Mp, i, Mp_j, opts, len_i, len_j);
      // print_raw_matrix(C.get(), len_i, len_j);
      // std::cout << fmt::format("len_i = {} len_j = {}\n", len_i, len_j) << std::endl;
    }

    DistanceIndexPairs wDistPairCPU = wasserstein_distances(Mp_j, opts, M, Mp, tryInds);
    // print_raw_matrix(wDistPairCPU.dists.data(), 1, wDistPairCPU.dists.size());
    assert(wDistPairCPU.dists.size() == tryInds.size());

#if defined(WITH_ARRAYFIRE)
    // Char is the internal representation of bool in ArrayFire
    std::vector<char> mopts;
    for (int j = 0; j < M.E_actual(); j++) {
      mopts.push_back(opts.metrics[j] == Metric::Diff);
    }

    af::array metricOpts(M.E_actual(), mopts.data());

    const ManifoldOnGPU gpuM = M.toGPU(false);
    const ManifoldOnGPU gpuMp = Mp.toGPU(false);

    // The next line currently crashes:
    DistanceIndexPairs wDistPairGPU = afWassersteinDistances(Mp_j, opts, M, Mp, gpuM, gpuMp, tryInds, metricOpts);
    assert(wDistPairGPU.dists.size() == tryInds.size());
#endif
  }
}

TEST_CASE("Coprediction and usable/library set/prediction sets", "[copredSets]")
{
  int E = 2;
  int tau = 1;
  int p = 1;

  std::vector<double> t = { 0, 1, 2, 3, 4, 5 };
  std::vector<double> x = { 1, NA, 3, 4, 5, 6 };
  std::vector<double> co_x = { 11, 12, 13, NA, 15, 16 };

  ManifoldGenerator generator(t, x, tau, p, {}, co_x);

  std::vector<bool> usable = generator.generate_usable(E);
  std::vector<bool> usableTrue = { false, false, false, true, true, true };
  require_vectors_match<bool>(usable, usableTrue);

  std::vector<bool> cousable = generator.generate_usable(E, true);
  std::vector<bool> cousableTrue = { false, true, true, false, false, true };
  require_vectors_match<bool>(cousable, cousableTrue);

  SECTION("Standard predictions")
  {
    Manifold M = generator.create_manifold(E, usable, false);
    Manifold Mp = generator.create_manifold(E, usable, true);
    check_usable_matches_prediction_set(usable, Mp);

    std::vector<std::vector<double>> M_true = { { 4, 3 }, { 5, 4 } };
    std::vector<double> y_true = { 5, 6 };
    require_manifolds_match(M, M_true, y_true);

    std::vector<std::vector<double>> Mp_true = { { 4, 3 }, { 5, 4 }, { 6, 5 } };
    std::vector<double> yp_true = { 5, 6, NA };
    require_manifolds_match(Mp, Mp_true, yp_true);
  }

  SECTION("Copredictions")
  {
    bool copredict = true;
    Manifold M_co = generator.create_manifold(E, usable, false, -1, copredict);
    Manifold Mp_co = generator.create_manifold(E, cousable, true, -1, copredict);

    std::vector<std::vector<double>> M_co_true = { { 4, 3 }, { 5, 4 } };
    std::vector<double> y_co_true = { 5, 6 };
    require_manifolds_match(M_co, M_co_true, y_co_true);

    std::vector<std::vector<double>> Mp_co_true = { { 12, 11 }, { 13, 12 }, { 16, 15 } };
    std::vector<double> yp_co_true = { 13, NA, NA };
    require_manifolds_match(Mp_co, Mp_co_true, yp_co_true);
  }
}