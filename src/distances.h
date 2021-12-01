#pragma once

#include "common.h"

std::vector<Metric> expand_metrics(const ManifoldGenerator& generator, int E, Distance distance,
                                   const std::vector<Metric>& metrics);

DistanceIndexPairs eager_lp_distances(int Mp_i, const Options& opts, const Manifold& M, const Manifold& Mp);

DistanceIndexPairs lazy_lp_distances(int Mp_i, const Options& opts, const Manifold& M, const Manifold& Mp);

DistanceIndexPairs wasserstein_distances(int Mp_i, const Options& opts, const Manifold& M, const Manifold& Mp);

#if defined(WITH_ARRAYFIRE)
DistanceIndexPairsOnGPU afLPDistances(const int numPredictions, const Options& opts, const ManifoldOnGPU& M,
                                      const ManifoldOnGPU& Mp, const af::array& metricOpts);
DistanceIndexPairs afWassersteinDistances(int Mp_i, const Options& opts, const Manifold& hostM, const Manifold& hostMp,
                                          const ManifoldOnGPU& M, const ManifoldOnGPU& Mp,
                                          const std::vector<int>& inpInds, const af::array& metricOpts);
#endif
