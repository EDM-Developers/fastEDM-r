#pragma once

#include "common.h"

std::vector<std::future<Prediction>> launch_task_group(const ManifoldGenerator& generator, Options opts,
                                                       const std::vector<int>& Es, const std::vector<int>& libraries,
                                                       int k, int numReps, int crossfold, bool explore, bool full,
                                                       bool saveFinalPredictions, bool saveFinalCoPredictions,
                                                       bool saveSMAPCoeffs, bool copredictMode,
                                                       const std::vector<bool>& usable, const std::string& rngState,
                                                       IO* io, bool keep_going(), void all_tasks_finished());

// Below are the 'private' members of edm.cpp; they are added here just so they can be accessed for testing.

using MatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using MatrixXi = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

std::future<Prediction> launch_edm_task(const ManifoldGenerator& generator, Options opts, int E,
                                        const std::vector<bool>& trainingRows, const std::vector<bool>& predictionRows,
                                        IO* io, bool keep_going(), void all_tasks_finished());

Prediction edm_task(const Options opts, const Manifold M, const Manifold Mp, const std::vector<bool> predictionRows,
                    IO* io, bool keep_going(), void all_tasks_finished());

void make_prediction(int Mp_i, const Options& opts, const Manifold& M, const Manifold& Mp, Eigen::Map<MatrixXd> ystar,
                     Eigen::Map<MatrixXi> rc, Eigen::Map<MatrixXd> coeffs, int* kUsed, bool keep_going());

std::vector<int> potential_neighbour_indices(int Mp_i, const Options& opts, const Manifold& M, const Manifold& Mp);

DistanceIndexPairs kNearestNeighbours(const DistanceIndexPairs& potentialNeighbours, int k);

void simplex_prediction(int Mp_i, int t, const Options& opts, const Manifold& M, const std::vector<double>& dists,
                        const std::vector<int>& kNNInds, Eigen::Map<MatrixXd> ystar, Eigen::Map<MatrixXi> rc,
                        int* kUsed);

void smap_prediction(int Mp_i, int t, const Options& opts, const Manifold& M, const Manifold& Mp,
                     const std::vector<double>& dists, const std::vector<int>& kNNInds, Eigen::Map<MatrixXd> ystar,
                     Eigen::Map<MatrixXd> coeffs, Eigen::Map<MatrixXi> rc, int* kUsed);

#if defined(WITH_ARRAYFIRE)
void af_make_prediction(const int numPredictions, const Options& opts, const Manifold& hostM, const Manifold& hostMp,
                        const ManifoldOnGPU& M, const ManifoldOnGPU& Mp, const af::array& metricOpts,
                        Eigen::Map<MatrixXd> ystar, Eigen::Map<MatrixXi> rc, Eigen::Map<MatrixXd> coeffs,
                        std::vector<int>& kUseds, bool keep_going());
#endif
