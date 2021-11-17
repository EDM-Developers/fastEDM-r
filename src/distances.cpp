
#include "distances.h"
#include "EMD_wrapper.h"

#define EIGEN_NO_DEBUG
#define EIGEN_DONT_PARALLELIZE
#include <Eigen/Dense>

#include <cmath> // for std::isnormal
#include <limits>
#include <vector>

#if defined(WITH_ARRAYFIRE)
#include <cuda_runtime.h>
#if WITH_GPU_PROFILING
#include <nvtx3/nvToolsExt.h>
#endif
#define AF_DEFINE_CUDA_TYPES
#include <af/cuda.h>
#include <arrayfire.h>

#include "lp_distance.cuh"

using af::array;
#endif

std::vector<Metric> expand_metrics(const ManifoldGenerator& generator, int E, Distance distance,
                                   const std::vector<Metric>& metrics)
{
  // Expand the 'metrics' vector now that we know the value of E.
  std::vector<Metric> expandedMetrics;

  // For the Wasserstein distance, it's more convenient to have one 'metric' for each variable (before taking lags).
  // However, for the L^1 / L^2 distances, it's more convenient to have one 'metric' for each individual
  // point of each observations, so metrics.size() == M.E_actual().
  if (distance == Distance::Wasserstein) {
    // Add a metric for the main variable and for the dt variable.
    // These are always treated as a continuous values (though perhaps in the future this will change).
    expandedMetrics.push_back(Metric::Diff);
    if (generator.E_dt(E) > 0) {
      expandedMetrics.push_back(Metric::Diff);
    }

    // Add in the metrics for the 'extra' variables as they were supplied to us.
    for (int k = 0; k < generator.numExtras(); k++) {
      expandedMetrics.push_back(metrics[k]);
    }
  } else {
    // Add metrics for the main variable and the dt variable and their lags.
    // These are always treated as a continuous values (though perhaps in the future this will change).
    for (int lagNum = 0; lagNum < E + generator.E_dt(E); lagNum++) {
      expandedMetrics.push_back(Metric::Diff);
    }

    // The user specified how to treat the extra variables.
    for (int k = 0; k < generator.numExtras(); k++) {
      int numLags = (k < generator.numExtrasLagged()) ? E : 1;
      for (int lagNum = 0; lagNum < numLags; lagNum++) {
        expandedMetrics.push_back(metrics[k]);
      }
    }
  }

  return expandedMetrics;
}

DistanceIndexPairs lp_distances(int Mp_i, const Options& opts, const Manifold& M, const Manifold& Mp,
                                std::vector<int> inpInds)
{
  std::vector<int> inds;
  std::vector<double> dists;

  // Compare every observation in the M manifold to the
  // Mp_i'th observation in the Mp manifold.
  for (int i : inpInds) {
    // Calculate the distance between M[i] and Mp[Mp_i]
    double dist_i = 0.0;

    // If we have panel data and the M[i] / Mp[Mp_j] observations come from different panels
    // then add the user-supplied penalty/distance for the mismatch.
    if (opts.panelMode && opts.idw > 0) {
      dist_i += opts.idw * (M.panel(i) != Mp.panel(Mp_i));
    }

    for (int j = 0; j < M.E_actual(); j++) {
      // Get the sub-distance between M[i,j] and Mp[Mp_i, j]
      double dist_ij;

      // If either of these values is missing, the distance from
      // M[i,j] to Mp[Mp_i, j] is opts.missingdistance.
      // However, if the user doesn't specify this, then the entire
      // M[i] to Mp[Mp_i] distance is set as missing.
      if ((M(i, j) == MISSING_D) || (Mp(Mp_i, j) == MISSING_D)) {
        if (opts.missingdistance == 0) {
          dist_i = MISSING_D;
          break;
        } else {
          dist_ij = opts.missingdistance;
        }
      } else { // Neither M[i,j] nor Mp[Mp_i, j] is missing.
        // How do we compare them? Do we treat them like continuous values and subtract them,
        // or treat them like unordered categorical variables and just check if they're the same?
        if (opts.metrics[j] == Metric::Diff) {
          dist_ij = M(i, j) - Mp(Mp_i, j);
        } else { // Metric::CheckSame
          dist_ij = (M(i, j) != Mp(Mp_i, j));
        }
      }

      if (opts.distance == Distance::MeanAbsoluteError) {
        dist_i += abs(dist_ij) / M.E_actual();
      } else { // Distance::Euclidean
        dist_i += dist_ij * dist_ij;
      }
    }

    if (dist_i != 0 && dist_i != MISSING_D) {
      if (opts.distance == Distance::MeanAbsoluteError) {
        dists.push_back(dist_i);
      } else { // Distance::Euclidean
        dists.push_back(sqrt(dist_i));
      }
      inds.push_back(i);
    }
  }

  return { inds, dists };
}

// This function compares the M(i,.) multivariate time series to the Mp(j,.) multivariate time series.
// The M(i,.) observation has data for E consecutive time points (e.g. time(i), time(i+1), ..., time(i+E-1)) and
// the Mp(j,.) observation corresponds to E consecutive time points (e.g. time(j), time(j+1), ..., time(j+E-1)).
// At each time instant we observe n >= 1 pieces of data.
// These may either be continuous data or unordered categorical data.
//
// The Wasserstein distance (using the 'curve-matching' strategy) is equivalent to the (minimum) cost of turning
// the first time series into the second time series. In a simple example, say E = 2 and n = 1, and
//         M(i,.) = [ 1, 2 ] and Mp(j,.) = [ 2, 2 ].
// To turn M(i,.) into Mp(j,.) the first element needs to be increased by 1, so the overall cost is
//         Wasserstein( M(i,.), Mp(j,.) ) = 1.
// The distance can also reorder the points, so for example say
//         M(i,.) = [ 1, 100 ] and Mp(j,.) = [ 100, 1 ].
// If we just change the 1 to 100 and the 100 to 1 then the cost of each is 99 + 99 = 198.
// However, Wasserstein can instead reorder these points at a cost of
//         Wasserstein( M(i,.), Mp(j,.) ) = 2 * gamma * (time(1)-time(2))
// so if the observations occur on a regular grid so time(i) = i then the distance will just be 2 * gamma.
//
// The return value of this function is a matrix which shows the pairwise costs associated to each
// potential Wasserstein solution. E.g. the (n,m) element of the returned matrix shows the cost
// of turning the individual point M(i, n) into Mp(j, m).
//
// When there are missing values in one or other observation, we can either ignore this time period
// and compute the Wasserstein for the mismatched regime where M(i,.) is of size len_i and Mp(j,.) is
// of size len_j, where len_i != len_j is possible. Alternatively, we can fill in the affected elements
// of the cost matrix with some user-supplied 'missingDistance' value and then len_i == len_j is upheld.
std::unique_ptr<double[]> wasserstein_cost_matrix(const Manifold& M, const Manifold& Mp, int i, int j,
                                                  const Options& opts, int& len_i, int& len_j)
{
  // The M(i,.) observation will be stored as one flat vector of length M.E_actual():
  // - the first M.E() observations will the lagged version of the main time series
  // - the next M.E() observations will be the lagged 'dt' time series (if it is included, i.e., if M.E_dt() > 0)
  // - the next n * M.E() observations will be the n lagged extra variables,
  //   so in total that is M.E_lagged_extras() = n * M.E() observations
  // - the remaining M.E_actual() - M.E() - M.E_dt() - M.E_lagged_extras() are the unlagged extras and the distance
  //   between those two vectors forms a kind of minimum distance which is added to the time-series curve matching
  //   Wasserstein distance.

  bool skipMissing = (opts.missingdistance == 0);

  auto M_i = M.laggedObsMap(i);
  auto Mp_j = Mp.laggedObsMap(j);

  auto M_i_missing = (M_i.array() == M.missing()).colwise().any();
  auto Mp_j_missing = (Mp_j.array() == Mp.missing()).colwise().any();

  if (skipMissing) {
    // N.B. Can't .sum() a vector of bools to count them in Eigen.
    len_i = M.E() - M_i_missing.count();
    len_j = Mp.E() - Mp_j_missing.count();
  } else {
    len_i = M.E();
    len_j = Mp.E();
  }

  double gamma = 1.0;
  if (M.E_dt() > 0) {
    // Imagine the M_i time series as a plot, and calculate the
    // aspect ratio of this plot, so we can rescale the time variable
    // to get the user-supplied aspect ratio.
    double minData = std::numeric_limits<double>::max();
    double maxData = std::numeric_limits<double>::min();
    double maxTime = 0.0;
    for (int t = 0; t < M_i.cols(); t++) {
      if (M_i(0, t) != MISSING_D) {
        if (M_i(0, t) < minData) {
          minData = M_i(0, t);
        }
        if (M_i(0, t) > maxData) {
          maxData = M_i(0, t);
        }
      }
      if (M_i(1, t) != MISSING_D && M_i(1, t) > maxTime) {
        maxTime = M_i(1, t);
      }
    }

    double epsilon = 1e-6; // Some small number in case the following ratio gets wildly large/small
    gamma = opts.aspectRatio * (maxData - minData + epsilon) / (maxTime + epsilon);
  }

  int timeSeriesDim = M_i.rows();

  double unlaggedDist = 0.0;
  int numUnlaggedExtras = M.E_extras() - M.E_lagged_extras();
  for (int e = 0; e < numUnlaggedExtras; e++) {
    double x = M.unlagged_extras(i, e), y = Mp.unlagged_extras(j, e);
    bool eitherMissing = (x == M.missing()) || (y == M.missing());

    if (eitherMissing) {
      unlaggedDist += opts.missingdistance;
    } else {
      if (opts.metrics[timeSeriesDim + e] == Metric::Diff) {
        unlaggedDist += abs(x - y);
      } else {
        unlaggedDist += x != y;
      }
    }
  }

  // If we have panel data and the M[i] / Mp[j] observations come from different panels
  // then add the user-supplied penalty/distance for the mismatch.
  if (opts.panelMode && opts.idw > 0) {
    unlaggedDist += opts.idw * (M.panel(i) != Mp.panel(j));
  }

  auto flatCostMatrix = std::make_unique<double[]>(len_i * len_j);
  std::fill_n(flatCostMatrix.get(), len_i * len_j, unlaggedDist);
  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> costMatrix(flatCostMatrix.get(),
                                                                                                len_i, len_j);
  for (int k = 0; k < timeSeriesDim; k++) {
    int n = 0;
    for (int nn = 0; nn < M_i.cols(); nn++) {
      if (skipMissing && M_i_missing[nn]) {
        continue;
      }

      int m = 0;

      for (int mm = 0; mm < Mp_j.cols(); mm++) {
        if (skipMissing && Mp_j_missing[mm]) {
          continue;
        }
        double dist;
        bool eitherMissing = M_i_missing[nn] || Mp_j_missing[mm];

        if (eitherMissing) {
          dist = opts.missingdistance;
        } else {
          if (opts.metrics[k] == Metric::Diff) {
            dist = abs(M_i(k, nn) - Mp_j(k, mm));
          } else {
            dist = M_i(k, nn) != Mp_j(k, mm);
          }
        }

        // For the time data, we add in the 'gamma' scaling factor calculated earlier
        if ((M.E_dt() > 0) && (k == 1)) {
          dist *= gamma;
        }

        costMatrix(n, m) += dist;

        m += 1;
      }

      n += 1;
    }
  }
  return flatCostMatrix;
}

// TODO: Subtract the D(x,x) and D(y,y) parts from this.
double approx_wasserstein(double* C, int len_i, int len_j, double eps, double stopErr)
{
  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> costMatrix(C, len_i, len_j);

  double r = 1.0 / len_i;
  double c = 1.0 / len_j;

  Eigen::MatrixXd K = Eigen::exp(-costMatrix.array() / eps);
  Eigen::MatrixXd Kp = len_i * K.array();

  Eigen::VectorXd u = Eigen::VectorXd::Ones(len_i) / len_i;
  Eigen::VectorXd v = Eigen::VectorXd::Ones(len_j) / len_j;

  int maxIter = 10000;
  for (int iter = 0; iter < maxIter; iter++) {

    v = c / (K.transpose() * u).array();
    u = 1.0 / (Kp * v).array();

    if (iter % 10 == 0) {
      // Compute right marginal (diag(u) K diag(v))^T1
      Eigen::VectorXd tempColSums = (u.asDiagonal() * K * v.asDiagonal()).colwise().sum();
      double LInfErr = (tempColSums.array() - c).abs().maxCoeff();
      if (LInfErr < stopErr) {
        break;
      }
    }
  }

  Eigen::MatrixXd transportPlan = u.asDiagonal() * K * v.asDiagonal();
  double dist = (transportPlan.array() * costMatrix.array()).sum();
  return dist;
}

double wasserstein(double* C, int len_i, int len_j)
{
  // Create vectors which are just 1/len_i and 1/len_j of length len_i and len_j.
  auto w_1 = std::make_unique<double[]>(len_i);
  std::fill_n(w_1.get(), len_i, 1.0 / len_i);
  auto w_2 = std::make_unique<double[]>(len_j);
  std::fill_n(w_2.get(), len_j, 1.0 / len_j);

  int maxIter = 10000;
  double cost;
  EMD_wrap(len_i, len_j, w_1.get(), w_2.get(), C, &cost, maxIter);
  return cost;
}

DistanceIndexPairs wasserstein_distances(int Mp_i, const Options& opts, const Manifold& M, const Manifold& Mp,
                                         std::vector<int> inpInds)
{
  std::vector<int> inds;
  std::vector<double> dists;

  // Compare every observation in the M manifold to the
  // Mp_i'th observation in the Mp manifold.
  for (int i : inpInds) {
    int len_i, len_j;
    auto C = wasserstein_cost_matrix(M, Mp, i, Mp_i, opts, len_i, len_j);

    if (len_i > 0 && len_j > 0) {
      double dist_i = wasserstein(C.get(), len_i, len_j);

      // Alternatively, the approximate version based on Sinkhorn's algorithm can be called with something like:
      // double dist_i = approx_wasserstein(C.get(), len_i, len_j, 0.1, 0.1)
      // In that case, the "std::isnormal" is really needed on the next line, as some
      // instability gives us some 'nan' distances using that method.

      if (dist_i != 0 && std::isnormal(dist_i)) {
        dists.push_back(dist_i);
        inds.push_back(i);
      }
    }
  }

  return { inds, dists };
}

/////////////////////////////////////////////////////////////// ArrayFire PORTED versions BEGIN HERE

#if defined(WITH_ARRAYFIRE)
DistanceIndexPairsOnGPU afLPDistances(const int npreds, const Options& opts, const ManifoldOnGPU& M,
                                      const ManifoldOnGPU& Mp, const af::array& metricOpts)
{
  constexpr bool useCustomKernel = true;

#if WITH_GPU_PROFILING
  auto range = nvtxRangeStartA(__FUNCTION__);
#endif

  const af_dtype cType = M.mdata.type();

  if (useCustomKernel) {
    af::array valids(M.numPoints, npreds, b8);
    af::array dists(M.numPoints, npreds, cType);
    if (cType == f64) {
      cuLPDistances(valids.device<char>(), dists.device<double>(), npreds, opts.distance == Distance::MeanAbsoluteError,
                    opts.panelMode, opts.idw, opts.missingdistance, M.E_actual, M.numPoints, M.mdata.device<double>(),
                    M.panel.device<int>(), Mp.mdata.device<double>(), Mp.panel.device<int>(), metricOpts.device<char>(),
                    afcu::getStream(0));
    } else if (cType == f32) {
      cuLPDistances(valids.device<char>(), dists.device<float>(), npreds, opts.distance == Distance::MeanAbsoluteError,
                    opts.panelMode, opts.idw, opts.missingdistance, M.E_actual, M.numPoints, M.mdata.device<float>(),
                    M.panel.device<int>(), Mp.mdata.device<float>(), Mp.panel.device<int>(), metricOpts.device<char>(),
                    afcu::getStream(0));
    }
    valids.unlock();
    dists.unlock();
    M.mdata.unlock();
    M.panel.unlock();
    Mp.mdata.unlock();
    Mp.panel.unlock();
    metricOpts.unlock();

#if WITH_GPU_PROFILING
    nvtxRangeEnd(range);
#endif
    return { valids, dists };
  } else {
    using af::array;
    using af::moddims;
    using af::select;
    using af::seq;
    using af::span;
    using af::sum;
    using af::tile;

    // Mp_i goes from 0 to npreds - 1
    // All coloumns of Manifold M are considered valid for batch operation

    const bool imdoZero = opts.missingdistance == 0;
    const int numLibraryPoints = M.numPoints;
    const int eacts = M.E_actual;

    array anyMCols, distsMat;
    {
      array predsM = tile(M.mdata, 1, 1, npreds);
      array predsMp = tile(moddims(Mp.mdata(span, seq(npreds)), eacts, 1, npreds), 1, numLibraryPoints);
      array diffMMp = predsM - predsMp;
      array compMMp = (predsM != predsMp).as(cType);
      array distMMp = select(tile(metricOpts, 1, numLibraryPoints, npreds), diffMMp, compMMp);
      array missing = predsM == MISSING_D || predsMp == MISSING_D;

      distsMat = (imdoZero ? distMMp : select(missing, opts.missingdistance, distMMp));
      anyMCols = anyTrue(missing, 0);

      if (opts.distance == Distance::MeanAbsoluteError) {
        distsMat = af::abs(distsMat) / double(eacts);
      } else {
        distsMat = distsMat * distsMat;
      }
    }
    if (opts.panelMode && opts.idw > 0) {
      array npPanelMp = tile(Mp.panel(seq(npreds)).T(), numLibraryPoints);
      array npPanelM = tile(M.panel, 1, npreds);
      array penalty = (opts.idw * (npPanelM != npPanelMp));
      array penalties = tile(moddims(penalty, 1, numLibraryPoints, npreds), eacts);

      distsMat += penalties;
    }
    array accDists = sum(distsMat, 0);
    array distances = select(anyMCols * imdoZero, double(MISSING_D), accDists);
    array valids = (distances != 0.0 && distances != double(MISSING_D));
    array dists = (opts.distance == Distance::MeanAbsoluteError ? distances : af::sqrt(distances));

    valids = moddims(valids, numLibraryPoints, npreds);
    dists = moddims(dists, numLibraryPoints, npreds);

#if WITH_GPU_PROFILING
    nvtxRangeEnd(range);
#endif
    return { valids, dists };
  }
}

array afWassersteinCostMatrix(const bool& skipMissing, const Options& opts, const array& metricOpts,
                              const ManifoldOnGPU& M, const array& M_i, const array& M_i_missing, const array& x,
                              const int& len_i, const array& Mp_j, const array& Mp_j_missing, const array& y,
                              const int& len_j, const bool arePanelIdsSame)
{
  using af::array;
  using af::constant;
  using af::dim4;
  using af::moddims;
  using af::seq;
  using af::span;
  using af::sum;
  using af::tile;
  using af::where;

  const af_dtype cType = M.mdata.type();

  double gamma = 1.0;
  if (M.E_dt > 0) {
    array firstColumn = M_i(span, 0);
    array nonMissings = firstColumn != MISSING_D;
    array validIndexs = where(nonMissings);
    array validValues = firstColumn(validIndexs);
    double minData = af::min<double>(validValues);
    double maxData = af::max<double>(validValues);
    double maxTime = af::max<double>(M_i(span, 1));

    // Some small number in case the following ratio gets wildly large/small
    constexpr double epsilon = 1e-6;

    gamma = opts.aspectRatio * (maxData - minData + epsilon) / (maxTime + epsilon);
  }

  const int timeSeriesDim = M_i.dims(1);
  double unlaggedDist = 0.0;
  {
    const int numUnlaggedExtras = M.E_extras - M.E_lagged_extras;

    array eitherMissing = (x == M.missing || y == M.missing);

    array cond = metricOpts(seq(timeSeriesDim, numUnlaggedExtras + timeSeriesDim - 1));
    array ulDists = (cond * af::abs(x - y) + (1 - cond) * (x != y).as(cType));
    array dists = (eitherMissing * opts.missingdistance + (1 - eitherMissing) * ulDists);
    unlaggedDist = sum<double>(dists);
  }

  if (opts.panelMode && opts.idw > 0) {
    unlaggedDist += opts.idw * arePanelIdsSame;
  }

  const seq timeSeries(timeSeriesDim);
  array costMatrix = af::constant(unlaggedDist, len_j, len_i, cType);

  array cmIsMetricDiff = metricOpts(timeSeries);
  cmIsMetricDiff = tile(moddims(cmIsMetricDiff, 1, 1, timeSeriesDim), len_j, len_i);

  if (skipMissing) {
    // In this case: Unless both entries are available, no need to process anything else
    array idxM_i_missing = where(!M_i_missing);
    array idxMp_j_missing = where(!Mp_j_missing);

    array Mp_j_k = moddims(Mp_j(idxMp_j_missing, timeSeries), len_j, 1, timeSeriesDim);
    array M_i_k = moddims(M_i(idxM_i_missing, timeSeries), 1, len_i, timeSeriesDim);
    array cmMp_j = tile(Mp_j_k, 1, len_i);
    array cmM_i = tile(M_i_k, len_j);
    array cmDiff = select(cmIsMetricDiff, af::abs(cmMp_j - cmM_i), (cmM_i != cmMp_j).as(cType));

    if (M.E_dt > 0) {
      // For time series k = 1, scale by gamma
      cmDiff(span, span, 1) *= gamma;
    }
    costMatrix += sum(cmDiff, 2); // Add results to unlaggedDist
  } else {
    array Mp_j_missingT = tile(Mp_j_missing, 1, M_i_missing.dims(0)); // cost matrix shape
    array M_i_missingT = tile(M_i_missing.T(), Mp_j_missing.dims(0)); // cost matrix shape
    array eitherMissing = (M_i_missingT || Mp_j_missingT);            // one of the entries missing
    array eitherMissingT = tile(eitherMissing, 1, 1, timeSeriesDim);  // [len_j len_i timeSeriesDim 1]

    array Mp_j_k = moddims(Mp_j(span, timeSeries), len_j, 1, timeSeriesDim);
    array M_i_k = moddims(M_i(span, timeSeries), 1, len_i, timeSeriesDim);
    array cmMp_j = tile(Mp_j_k, 1, len_i);
    array cmM_i = tile(M_i_k, len_j);
    array cmDiff = select(cmIsMetricDiff, af::abs(cmMp_j - cmM_i), (cmM_i != cmMp_j).as(cType));
    array cmDist = select(eitherMissingT, opts.missingdistance, cmDiff);

    if (M.E_dt > 0) {
      // For time series k = 1, scale by gamma
      cmDist(span, span, 1) *= gamma;
    }
    costMatrix += sum(cmDist, 2); // Add results to unlaggedDist
  }

  return costMatrix;
}

DistanceIndexPairs afWassersteinDistances(int Mp_i, const Options& opts, const Manifold& hostM, const Manifold& hostMp,
                                          const ManifoldOnGPU& M, const ManifoldOnGPU& Mp,
                                          const std::vector<int>& inpInds, const af::array& metricOpts)
{
  using af::anyTrue;
  using af::seq;

  const bool skipMissing = (opts.missingdistance == 0);

  // Precompute values that are iteration invariant
  const seq mpjRange0(Mp.E_x);
  const seq mpjRange1(Mp_i, Mp_i + 1 + (Mp.E_dt > 0) + Mp.E_lagged_extras / Mp.E_x);

  const array Mp_j = Mp.mdata(mpjRange0, mpjRange1);

  const array Mp_j_missing = anyTrue(Mp_j == Mp.missing, 1);

  const int numUnlaggedExtrasEnd = M.E_extras - M.E_lagged_extras;

  const seq xseq(M.E_x + M.E_dt + M.E_lagged_extras, M.E_x + M.E_dt + M.E_lagged_extras + numUnlaggedExtrasEnd - 1);
  const seq yseq(Mp.E_x + Mp.E_dt + Mp.E_lagged_extras,
                 Mp.E_x + Mp.E_dt + Mp.E_lagged_extras + numUnlaggedExtrasEnd - 1);

  const array y = Mp.mdata(yseq, Mp_i);
  const int len_j = (skipMissing ? Mp.E_x - af::sum<int>(Mp_j_missing) : Mp.E_x);

  // Since both len_i and len_j should be greater than zero
  // to collect valid indices and respective distances, just return empty handed
  if (len_j <= 0) {
    return {};
  }

  // Return Items
  std::vector<int> inds;
  std::vector<double> dists;

  // Compare every observation in the M manifold to the Mp_i'th observation in the Mp manifold.
  for (int i : inpInds) {
    const seq miRange0(M.E_x);
    const seq miRange1(i, i + 1 + (M.E_dt > 0) + M.E_lagged_extras / M.E_x);

    array M_i = M.mdata(miRange0, miRange1);

    array M_i_missing = anyTrue(M_i == M.missing, 1);

    array x = M.mdata(xseq, i);

    const int len_i = (skipMissing ? M.E_x - af::sum<int>(M_i_missing) : M.E_x);

    if (len_i > 0) { // Short-check for len_j already passed
      // TODO I think afWassersteinCostMatrix can be further vectorized to run for all i's
      array cm = afWassersteinCostMatrix(skipMissing, opts, metricOpts, M, M_i, M_i_missing, x, len_i, Mp_j,
                                         Mp_j_missing, y, len_j, (hostM.panel(i) != hostMp.panel(Mp_i)));
      std::vector<double> C(cm.elements());
      cm.host(C.data());
      double dist_i = wasserstein(C.data(), len_i, len_j);

      // Alternative: approximate version based on Sinkhorn's algorithm
      // double dist_i = approx_wasserstein(C.get(), len_i, len_j, 0.1, 0.1)
      // In that case, the "std::isnormal" is really needed on the next line, as some
      // instability gives us some 'nan' distances using that method.
      if (dist_i != 0 && std::isnormal(dist_i)) {
        dists.push_back(dist_i);
        inds.push_back(i);
      }
    }
  }

  return { inds, dists };
}
#endif
