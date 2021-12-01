#pragma once

#include "common.h"

std::vector<std::future<PredictionResult>> launch_task_group(
  const std::shared_ptr<ManifoldGenerator> generator, Options opts, const std::vector<int>& Es,
  const std::vector<int>& libraries, int k, int numReps, int crossfold, bool explore, bool full, bool shuffle,
  bool saveFinalPredictions, bool saveFinalCoPredictions, bool saveSMAPCoeffs, bool copredictMode,
  const std::vector<bool>& usable, const std::string& rngState, IO* io, bool keep_going(), void all_tasks_finished());
