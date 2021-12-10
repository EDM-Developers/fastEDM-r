#include "common.h"

#ifdef JSON

void to_json(json& j, const Algorithm& a)
{
  if (a == Algorithm::Simplex) {
    j = "Simplex";
  } else {
    j = "SMap";
  }
}
void from_json(const json& j, Algorithm& a)
{
  if (j == "Simplex") {
    a = Algorithm::Simplex;
  } else {
    a = Algorithm::SMap;
  }
}

void to_json(json& j, const Distance& d)
{
  if (d == Distance::MeanAbsoluteError) {
    j = "MeanAbsoluteError";
  } else if (d == Distance::Euclidean) {
    j = "Euclidean";
  } else {
    j = "Wasserstein";
  }
}

void from_json(const json& j, Distance& d)
{
  if (j == "MeanAbsoluteError") {
    d = Distance::MeanAbsoluteError;
  } else if (j == "Euclidean") {
    d = Distance::Euclidean;
  } else {
    d = Distance::Wasserstein;
  }
}

void to_json(json& j, const Metric& m)
{
  if (m == Metric::Diff) {
    j = "Diff";
  } else {
    j = "CheckSame";
  }
}

void from_json(const json& j, Metric& m)
{
  if (j == "Diff") {
    m = Metric::Diff;
  } else {
    m = Metric::CheckSame;
  }
}

void to_json(json& j, const Options& o)
{
  j = json{ { "copredict", o.copredict },
            { "forceCompute", o.forceCompute },
            { "savePrediction", o.savePrediction },
            { "saveSMAPCoeffs", o.saveSMAPCoeffs },
            { "k", o.k },
            { "nthreads", o.nthreads },
            { "missingdistance", o.missingdistance },
            { "panelMode", o.panelMode },
            { "idw", o.idw },
            { "thetas", o.thetas },
            { "algorithm", o.algorithm },
            { "taskNum", o.taskNum },
            { "numTasks", o.numTasks },
            { "configNum", o.configNum },
            { "calcRhoMAE", o.calcRhoMAE },
            { "aspectRatio", o.aspectRatio },
            { "distance", o.distance },
            { "metrics", o.metrics },
            { "cmdLine", o.cmdLine } };
}

void from_json(const json& j, Options& o)
{
  j.at("copredict").get_to(o.copredict);
  j.at("forceCompute").get_to(o.forceCompute);
  j.at("savePrediction").get_to(o.savePrediction);
  j.at("saveSMAPCoeffs").get_to(o.saveSMAPCoeffs);
  j.at("k").get_to(o.k);
  j.at("nthreads").get_to(o.nthreads);
  j.at("missingdistance").get_to(o.missingdistance);
  j.at("panelMode").get_to(o.panelMode);
  j.at("idw").get_to(o.idw);
  j.at("thetas").get_to(o.thetas);
  j.at("algorithm").get_to(o.algorithm);
  j.at("taskNum").get_to(o.taskNum);
  j.at("numTasks").get_to(o.numTasks);
  j.at("configNum").get_to(o.configNum);
  j.at("calcRhoMAE").get_to(o.calcRhoMAE);
  j.at("aspectRatio").get_to(o.aspectRatio);
  j.at("distance").get_to(o.distance);
  j.at("metrics").get_to(o.metrics);
  j.at("cmdLine").get_to(o.cmdLine);
}

void to_json(json& j, const PredictionStats& s)
{
  j = json{ { "mae", s.mae }, { "rho", s.rho } };
}

void from_json(const json& j, PredictionStats& s)
{
  j.at("mae").get_to(s.mae);
  j.at("rho").get_to(s.rho);
}

void to_json(json& j, const PredictionResult& p)
{
  std::vector<double> predictionsVec, coeffsVec;
  if (p.predictions != nullptr) {
    predictionsVec = std::vector<double>(p.predictions.get(), p.predictions.get() + p.numThetas * p.numPredictions);
  }
  if (p.coeffs != nullptr) {
    coeffsVec = std::vector<double>(p.coeffs.get(), p.coeffs.get() + p.numPredictions * p.numCoeffCols);
  }

  j = json{ { "rc", p.rc },
            { "numThetas", p.numThetas },
            { "numPredictions", p.numPredictions },
            { "numCoeffCols", p.numCoeffCols },
            { "predictions", predictionsVec },
            { "coeffs", coeffsVec },
            { "stats", p.stats },
            { "predictionRows", p.predictionRows },
            { "kMin", p.kMin },
            { "kMax", p.kMax },
            { "cmdLine", p.cmdLine },
            { "configNum", p.configNum } };
}

void from_json(const json& j, PredictionResult& p)
{
  j.at("rc").get_to(p.rc);
  j.at("numThetas").get_to(p.numThetas);
  j.at("numPredictions").get_to(p.numPredictions);
  j.at("numCoeffCols").get_to(p.numCoeffCols);
  j.at("predictionRows").get_to(p.predictionRows);
  j.at("stats").get_to(p.stats);
  j.at("kMin").get_to(p.kMin);
  j.at("kMax").get_to(p.kMax);
  j.at("cmdLine").get_to(p.cmdLine);
  j.at("configNum").get_to(p.configNum);

  // TODO: Test this coeffs/predictions loading works as expected
  std::vector<double> predictions = j.at("predictions");
  if (predictions.size()) {
    p.predictions = std::make_unique<double[]>(predictions.size());
    for (int i = 0; i < predictions.size(); i++) {
      p.predictions[i] = predictions[i];
    }
  } else {
    p.predictions = nullptr;
  }

  std::vector<double> coeffs = j.at("coeffs");
  if (coeffs.size()) {
    p.coeffs = std::make_unique<double[]>(coeffs.size());

    for (int i = 0; i < coeffs.size(); i++) {
      p.coeffs[i] = coeffs[i];
    }
  } else {
    p.coeffs = nullptr;
  }
}

#endif