#pragma once

#include "common.h"

DistanceIndexPairs lp_distances(int Mp_i, const Options& opts, const Manifold& M, const Manifold& Mp,
                                std::vector<int> inds);
DistanceIndexPairs wasserstein_distances(int Mp_i, const Options& opts, const Manifold& M, const Manifold& Mp,
                                         std::vector<int> inds);

#if defined(WITH_ARRAYFIRE)
DistanceIndexPairsOnGPU afLPDistances(const int numPredictions, const Options& opts, const ManifoldOnGPU& M,
                                      const ManifoldOnGPU& Mp, const af::array& metricOpts);
DistanceIndexPairs afWassersteinDistances(int Mp_i, const Options& opts, const Manifold& hostM, const Manifold& hostMp,
                                          const ManifoldOnGPU& M, const ManifoldOnGPU& Mp,
                                          const std::vector<int>& inpInds, const af::array& metricOpts);
#endif
