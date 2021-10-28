#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one cpp file
#include <catch2/catch.hpp>

#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
#endif
#include <fmt/format.h>

#include "edm.h"
#include "manifold.h"

const double NA = MISSING_D;

// Add in some function declarations for 'private' functions not listed
// in the relevant header files.

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
  for (int i = 0; i < M.nobs(); i++) {
    for (int j = 0; j < M.E_actual(); j++) {
      std::cout << stringVersion(M(i, j)) << "\t";
    }
    std::cout << "\t|\t" << stringVersion(M.y(i)) << "\n";
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
  REQUIRE(M.nobs() == M_true.size());  // TODO: Rename this to numPoints
  REQUIRE(M.ySize() == y_true.size()); // Shouldn't this always be the same as M.nobs()?
  REQUIRE(M.E_actual() == M_true[0].size());

  for (int i = 0; i < M.nobs(); i++) {
    CAPTURE(i);
    for (int j = 0; j < M.E_actual(); j++) {
      CAPTURE(j);
      REQUIRE(M(i, j) == M_true[i][j]);
    }
    REQUIRE(M.y(i) == y_true[i]);
  }
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
    std::vector<bool> usableTrue = { false, true, true, false };
    require_vectors_match<bool>(usable, usableTrue);

    Manifold M = generator.create_manifold(E, usable, false, false, false);
    std::vector<std::vector<double>> M_true = { { 12.0, 11.0 }, { 13.0, 12.0 } };
    std::vector<double> y_true = { 13.0, 14.0 };
    require_manifolds_match(M, M_true, y_true);
  }

  SECTION("Manifold with dt (not allowing missing)")
  {
    // TODO: This test is a bit fake; edm.ado would not allow dt to be
    // applied when there's no gaps in the time variable.
    bool dtMode = true;
    double dtWeight = 1.0;
    bool allowMissing = false;
    ManifoldGenerator generator(t, x, tau, p, {}, {}, {}, {}, 0, dtMode, true, false, allowMissing);

    std::vector<bool> usable = generator.generate_usable(E);
    std::vector<bool> usableTrue = { false, true, true, false };
    require_vectors_match<bool>(usable, usableTrue);

    Manifold M = generator.create_manifold(E, usable, false, false, dtWeight, false);
    std::vector<std::vector<double>> M_true = { { 12.0, 11.0, 1.0, 1.0 }, { 13.0, 12.0, 1.0, 1.0 } };
    std::vector<double> y_true = { 13.0, 14.0 };
    require_manifolds_match(M, M_true, y_true);
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
    std::vector<bool> usableTrue = { false, false, false, false, false, false };
    require_vectors_match<bool>(usable, usableTrue);

    Manifold M = generator.create_manifold(E, usable, false, false, false);
    REQUIRE(M.nobs() == 0);
    REQUIRE(M.ySize() == 0);
    REQUIRE(M.E_actual() == 2);
  }

  SECTION("dt")
  {
    bool dtMode = true;
    double dtWeight = 1.0;
    bool allowMissing = false;
    ManifoldGenerator generator(t, x, tau, p, {}, {}, {}, {}, 0, dtMode, true, false, allowMissing);

    std::vector<int> obsNums = { 0, 1, -1, 2, 3, 4 };
    for (int i = 0; i < obsNums.size(); i++) {
      CAPTURE(i);
      REQUIRE(generator.get_observation_num(i) == obsNums[i]);
    }

    std::vector<bool> usable = generator.generate_usable(E);
    std::vector<bool> usableTrue = { false, true, false, true, true, false };
    require_vectors_match<bool>(usable, usableTrue);

    Manifold M = generator.create_manifold(E, usable, false, false, dtWeight, false);
    std::vector<std::vector<double>> M_true = { { 12.0, 11.0, 2.0, 1.5 },
                                                { 14.0, 12.0, 0.5, 2.0 },
                                                { 15.0, 14.0, 1.0, 0.5 } };
    std::vector<double> y_true = { 14.0, 15.0, 16.0 };
    require_manifolds_match(M, M_true, y_true);
  }

  SECTION("dt and allowingmissing")
  {
    bool dtMode = true;
    double dtWeight = 1.0;
    bool allowMissing = true;
    ManifoldGenerator generator(t, x, tau, p, {}, {}, {}, {}, 0, dtMode, true, false, allowMissing);

    std::vector<int> obsNums = { 0, 1, 2, 3, 4, 5 };
    for (int i = 0; i < obsNums.size(); i++) {
      CAPTURE(i);
      REQUIRE(generator.get_observation_num(i) == obsNums[i]);
    }

    std::vector<bool> usable = generator.generate_usable(E);
    std::vector<bool> usableTrue = { true, false, true, true, true, false };
    require_vectors_match<bool>(usable, usableTrue);

    Manifold M = generator.create_manifold(E, usable, false, false, dtWeight, false);
    std::vector<std::vector<double>> M_true = {
      { 11.0, NA, 1.5, NA }, { NA, 12.0, 1.5, 0.5 }, { 14.0, NA, 0.5, 1.5 }, { 15.0, 14.0, 1.0, 0.5 }
    };
    std::vector<double> y_true = { 12.0, 14.0, 15.0, 16.0 };
    require_manifolds_match(M, M_true, y_true);
  }

  SECTION("reldt")
  {
    bool dtMode = true, reldtMode = true;
    double dtWeight = 1.0;
    bool allowMissing = false;
    ManifoldGenerator generator(t, x, tau, p, {}, {}, {}, {}, 0, dtMode, true, reldtMode, allowMissing);

    std::vector<int> obsNums = { 0, 1, -1, 2, 3, 4 };
    for (int i = 0; i < obsNums.size(); i++) {
      CAPTURE(i);
      REQUIRE(generator.get_observation_num(i) == obsNums[i]);
    }

    std::vector<bool> usable = generator.generate_usable(E);
    std::vector<bool> usableTrue = { false, true, false, true, true, false };
    require_vectors_match<bool>(usable, usableTrue);

    Manifold M = generator.create_manifold(E, usable, false, false, dtWeight, false);
    std::vector<std::vector<double>> M_true = { { 12.0, 11.0, 2.0, 3.5 },
                                                { 14.0, 12.0, 0.5, 2.5 },
                                                { 15.0, 14.0, 1.0, 1.5 } };
    std::vector<double> y_true = { 14.0, 15.0, 16.0 };
    require_manifolds_match(M, M_true, y_true);
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
    ManifoldGenerator generator(t, x, tau, p, {}, {}, {}, {}, 0, dtMode, true, false, allowMissing);

    std::vector<bool> usable = generator.generate_usable(E);
    std::vector<bool> usableTrue = { true, false, true, true, true, false };
    require_vectors_match<bool>(usable, usableTrue);

    Manifold M = generator.create_manifold(E, usable, false, false, dtWeight, false);
    std::vector<std::vector<double>> M_true = {
      { 11.0, NA, 1.5, NA },
      { NA, 11.0, 1.5, 2.0 },
      { 14.0, 12.0, 0.5, 2.0 },
      { 15.0, NA, 1.0, 2.0 },
    };
    std::vector<double> y_true = { 12.0, 14.0, 15.0, 16.0 };
    require_manifolds_match(M, M_true, y_true);
  }

  SECTION("Not allowing missing values")
  {
    bool dtMode = true;
    double dtWeight = 1.0;
    bool allowMissing = false;
    ManifoldGenerator generator(t, x, tau, p, {}, {}, {}, {}, 0, dtMode, true, false, allowMissing);

    std::vector<bool> usable = generator.generate_usable(E);
    std::vector<bool> usableTrue = { false, false, false, true, true, false };
    require_vectors_match<bool>(usable, usableTrue);

    Manifold M = generator.create_manifold(E, usable, false, false, dtWeight, false);
    std::vector<std::vector<double>> M_true = { { 14.0, 11.0, 0.5, 3.5 }, { 15.0, 12.0, 1.0, 2.5 } };
    std::vector<double> y_true = { 15.0, 16.0 };
    require_manifolds_match(M, M_true, y_true);
  }
}

TEST_CASE("Check negative times work", "[negativeTimes]")
{
  std::vector<double> t = { -9.0, -7.5, -7.0, -6.5, -5.0, -4.0 };
  std::vector<double> x = { 11, 12, NA, 14, 15, 16 };

  int tau = 1;
  int p = 1;
  int E = 2;

  std::vector<std::vector<double>> extras;
  int numExtrasLagged = 0;

  ManifoldGenerator generator(t, x, tau, p);

  REQUIRE(generator.calculate_time_increment() == 0.5);

  std::vector<int> obsNums = { 0, 3, 4, 5, 8, 10 };
  for (int i = 0; i < obsNums.size(); i++) {
    CAPTURE(i);
    REQUIRE(generator.get_observation_num(i) == obsNums[i]);
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

  bool dt = true, dt0 = true, reldt = true, allowMissing = true;
  ManifoldGenerator generator(t, x, tau, p, xmap, co_x, {}, {}, 0, dt, dt0, reldt, allowMissing);

  std::vector<bool> usable = generator.generate_usable(E);
  std::vector<bool> usableTrue = { true, false, false, true, true };
  require_vectors_match<bool>(usable, usableTrue);

  double dtWeight = 1.0;
  bool copredict = true;
  Manifold M = generator.create_manifold(E, usable, copredict, false, dtWeight);
  Manifold Mp = generator.create_manifold(E, usable, copredict, true, dtWeight);

  std::vector<std::vector<double>> M_true = {
    { 1, NA, NA, NA, NA, 0, NA, NA, NA, NA },
    { NA, NA, 2, 1, NA, 0, 1, 2, 3, NA },
    { 5, NA, NA, 2, 1, 0, 1, 2, 3, 4 },
  };
  std::vector<double> y_true = { 1, 4, 5 };
  require_manifolds_match(M, M_true, y_true);

  std::vector<std::vector<double>> Mp_true = {
    { 1, NA, NA, NA, NA, 0, NA, NA, NA, NA },
    { 4, NA, NA, 1, NA, 0, 1, 2, 3, NA },
    { 5, 4, NA, NA, 1, 0, 1, 2, 3, 4 },
  };
  std::vector<double> yp_true = y_true;
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
    int i = 2, j = 2;

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

    for (int i = 0; i < C_true.size(); i++) {
      REQUIRE(C[i] == C_true[i]);
    }
  }

  SECTION("Multiple Wasserstein distances")
  {

    int Mp_j = 2;

    std::vector<int> tryInds = potential_neighbour_indices(Mp_j, opts, M, Mp);
    for (int i = 0; i < M.nobs(); i++) {
      std::cout << fmt::format("potential_neighbour_indices[{}] = {}\n", i, tryInds[i]);
    }

    for (int i = 0; i < M.nobs(); i++) {
      std::cout << fmt::format("Cost matrix M_i={}\n", i);

      auto M_i_map = M.laggedObsMap(i);
      auto Mp_j_map = Mp.laggedObsMap(Mp_j);

      std::cout << "M_i_map={}\n";
      print_eig_matrix(M_i_map);

      std::cout << "Mp_j_map={}\n";
      print_eig_matrix(Mp_j_map);

      int len_i, len_j;
      std::unique_ptr<double[]> C = wasserstein_cost_matrix(M, Mp, i, Mp_j, opts, len_i, len_j);
      print_raw_matrix(C.get(), len_i, len_j);

      std::cout << fmt::format("len_i = {} len_j = {}\n", len_i, len_j);

      std::cout << std::endl;
    }

    DistanceIndexPairs wDistPairCPU = wasserstein_distances(Mp_j, opts, M, Mp, tryInds);
    print_raw_matrix(wDistPairCPU.dists.data(), 1, wDistPairCPU.dists.size());
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

TEST_CASE("Coprediction and usable/training set/prediction sets", "[copredSets]")
{
  int E = 2;
  int tau = 1;
  int p = 1;

  std::vector<double> t = { 0, 1, 2, 3, 4, 5 };
  std::vector<double> x = { 1, NA, 3, 4, 5, 6 };
  std::vector<double> co_x = { 11, 12, 13, NA, 15, 16 };

  ManifoldGenerator generator(t, x, tau, p, {}, co_x);

  std::vector<bool> usable = generator.generate_usable(E);
  std::vector<bool> usableTrue = { false, false, false, true, true, false };
  require_vectors_match<bool>(usable, usableTrue);

  std::vector<bool> cousable = generator.generate_usable(E, true);
  std::vector<bool> cousableTrue = { false, true, false, false, false, false };
  require_vectors_match<bool>(cousable, cousableTrue);

  SECTION("Standard predictions")
  {
    bool copredict = false;
    Manifold M = generator.create_manifold(E, usable, copredict, false);
    Manifold Mp = generator.create_manifold(E, usable, copredict, true);

    std::vector<std::vector<double>> M_true = { { 4, 3 }, { 5, 4 } };
    std::vector<double> y_true = { 5, 6 };
    require_manifolds_match(M, M_true, y_true);

    std::vector<std::vector<double>> Mp_true = { { 4, 3 }, { 5, 4 } };
    std::vector<double> yp_true = { 5, 6 };
    require_manifolds_match(Mp, Mp_true, yp_true);
  }

  SECTION("Copredictions")
  {
    bool copredict = true;
    Manifold M_co = generator.create_manifold(E, usable, copredict, false);
    Manifold Mp_co = generator.create_manifold(E, cousable, copredict, true);

    std::vector<std::vector<double>> M_co_true = { { 4, 3 }, { 5, 4 } };
    std::vector<double> y_co_true = { 5, 6 };
    require_manifolds_match(M_co, M_co_true, y_co_true);

    std::vector<std::vector<double>> Mp_co_true = { { 12, 11 } };
    std::vector<double> yp_co_true = { 13 };
    require_manifolds_match(Mp_co, Mp_co_true, yp_co_true);
  }
}