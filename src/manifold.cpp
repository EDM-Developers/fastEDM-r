#pragma warning(disable : 4018)

#include "manifold.h"

// Recursive function to return gcd of a and b
// Lifted from https://www.geeksforgeeks.org/program-find-gcd-floating-point-numbers/
double gcd(double a, double b)
{
  if (a < b)
    return gcd(b, a);

  // base case
  if (fabs(b) < 0.001)
    return a;

  else
    return (gcd(b, a - floor(a / b) * b));
}

double ManifoldGenerator::calculate_time_increment() const
{
  // Find the units which time is measured in.
  // E.g. if time variables are 1, 2, 3, ... then the 'unit' is 1
  // Whereas if time is like 1000, 2000, 4000, 20000 ... then the 'unit' is perhaps 1000.
  double unit = -1;

  // Go through the supplied time index and find the greatest common divisor of the differences between consecutive time
  // points.
  for (int i = 1; i < _t.size(); i++) {

    double timeDiff = _t[i] - _t[i - 1];

    // In the panel data case, we may get consecutive times which are negative at the boundary of panels.
    if (timeDiff <= 0 || _t[i] == MISSING_D || _t[i - 1] == MISSING_D) {
      continue;
    }

    // For the first time, just replace sentinel value with the time difference.
    if (unit < 0) {
      unit = timeDiff;
      continue;
    }

    unit = gcd(timeDiff, unit);
  }

  return unit;
}

void ManifoldGenerator::setup_observation_numbers()
{
  _observation_number.clear();

  if (!(_dt || _reldt)) {
    // In normal situations (non-dt)
    double unit = calculate_time_increment();
    double minT = *std::min_element(_t.begin(), _t.end());
    // Note, missing values are huge (1e100 or so) so minT won't be affected by missing times.

    // Create a time index which is a discrete count of the number of 'unit' time units.
    for (int i = 0; i < _t.size(); i++) {
      if (_t[i] != MISSING_D) {
        _observation_number.push_back(std::round((_t[i] - minT) / unit));
      } else {
        _observation_number.push_back(MISSING_I);
      }
    }
  } else {
    // In 'dt' mode
    int countUp = 0;
    for (int i = 0; i < _t.size(); i++) {
      if (_t[i] != MISSING_D && (_allow_missing || (_x[i] != MISSING_D))) { // TODO: What about co_x missing here?
        _observation_number.push_back(countUp);
        countUp += 1;
      } else {
        _observation_number.push_back(MISSING_I);
      }
    }
  }
}

bool ManifoldGenerator::find_observation_num(int target, int& k, int direction, int panel) const
{
  // Loop either forward or back until we find the right index or give up.
  while (k >= 0 && k < _observation_number.size()) {

    if (_panel_mode) {
      // Skip over garbage rows which don't have a panel id recorded.
      if (_panelIDs[k] == MISSING_I) {
        k += direction;
        continue;
      }

      // If in panel mode, make sure we don't wander over a panel boundary.
      if (panel != _panelIDs[k]) {
        return false;
      }
    }

    // Skip over garbage rows which don't have a time recorded.
    if (_observation_number[k] == MISSING_I) {
      k += direction;
      continue;
    }

    // If we found the desired row at index k then stop here and report the success.
    if (_observation_number[k] == target) {
      return true;
    }

    // If we've gone past it & therefore this target doesn't exist, give up.
    if (direction > 0 && _observation_number[k] > target) {
      return false;
    }
    if (direction < 0 && _observation_number[k] < target) {
      return false;
    }

    k += direction;
  }

  return false;
}

#if defined(WITH_ARRAYFIRE)
ManifoldOnGPU Manifold::toGPU(const bool useFloat) const
{
  using af::array;

  if (useFloat) {
    return ManifoldOnGPU{ array(_E_actual, _numPoints, _flat.get()).as(f32),
                          (_targets.size() > 0 ? array(_numPoints, _targets.data()) : array()).as(f32),
                          (_panelIDs.size() > 0 ? array(_numPoints, _panelIDs.data()) : array()),
                          _numPoints,
                          _E_x,
                          _E_dt,
                          _E_extras,
                          _E_lagged_extras,
                          _E_actual,
                          MISSING_F };
  } else {
    return ManifoldOnGPU{ array(_E_actual, _numPoints, _flat.get()),
                          (_targets.size() > 0 ? array(_numPoints, _targets.data()) : array()),
                          (_panelIDs.size() > 0 ? array(_numPoints, _panelIDs.data()) : array()),
                          _numPoints,
                          _E_x,
                          _E_dt,
                          _E_extras,
                          _E_lagged_extras,
                          _E_actual,
                          MISSING_D };
  }
}
#endif

void ManifoldGenerator::fill_in_point(int i, int E, bool copredictionMode, bool predictionSet, double dtWeight,
                                      double* point) const
{
  int panel = _panel_mode ? _panelIDs[i] : -1;
  bool use_co_x = copredictionMode && predictionSet;

  double tPred;
  if (_dt || _reldt) {
    int targetIndex = i;
    get_target(i, copredictionMode, predictionSet, targetIndex);
    tPred = (targetIndex >= 0) ? _t[targetIndex] : MISSING_D;
  }

  // For obs i, which indices correspond to looking back 0, tau, ..., (E-1)*tau observations.
  int laggedIndex = i;
  int pointStartObsNum = _observation_number[i];

  // Start by going back one index
  int k = i - 1;

  int prevLaggedIndex = laggedIndex;

  for (int j = 0; j < E; j++) {

    bool foundLaggedObs = (j == 0);

    if (!foundLaggedObs) {
      // Find the discrete time we're searching for.
      int targetObsNum = pointStartObsNum - j * _tau;

      if (find_observation_num(targetObsNum, k, -1, panel)) {
        foundLaggedObs = true;
        laggedIndex = k;
      }
    }

    // Fill in the lagged embedding of x (or co_x) in the first columns
    if (!foundLaggedObs) {
      point[j] = MISSING_D;
    } else if (use_co_x) {
      point[j] = _co_x[laggedIndex];
    } else {
      point[j] = _x[laggedIndex];
    }

    // Put the lagged embedding of dt in the next columns
    if (_dt || _reldt) {
      if (!foundLaggedObs) {
        point[j + E] = MISSING_D;
      } else {
        double tNow = _t[laggedIndex];
        if (j == 0 || _reldt) {
          if (tNow != MISSING_D && tPred != MISSING_D) {
            point[E + j] = dtWeight * (tPred - tNow);
          } else {
            point[E + j] = MISSING_D;
          }
        } else {
          double tNext = _t[prevLaggedIndex];
          if (tNext != MISSING_D && tNow != MISSING_D) {
            point[E + j] = dtWeight * (tNext - tNow);
          } else {
            point[E + j] = MISSING_D;
          }
        }
      }
    }

    // Finally put the extras in the last columns
    for (int k = 0; k < _num_extras; k++) {
      double extra_k = foundLaggedObs ? _extras[k][laggedIndex] : MISSING_D;
      if (k < _num_extras_lagged) {
        // Adding the lagged extras
        point[E + E_dt(E) + k * E + j] = extra_k;
      } else if (j == 0) {
        // Add in any unlagged extras as the end
        point[E + E_dt(E) + _num_extras_lagged * E + (k - _num_extras_lagged)] = extra_k;
      }
    }

    prevLaggedIndex = laggedIndex;
  }
}

double ManifoldGenerator::get_dt(int i, bool copredictionMode, bool predictionSet, double dtWeight) const
{
  int panel = _panel_mode ? _panelIDs[i] : -1;

  // For obs i, which indices correspond to looking back 0, tau, ..., (E-1)*tau observations.
  int laggedIndex, prevLaggedIndex = i;
  int pointStartObsNum = _observation_number[i];

  // Start by going back one index
  int k = i - 1;

  // Find the discrete time we're searching for.
  int targetObsNum = pointStartObsNum - _tau;

  if (find_observation_num(targetObsNum, k, -1, panel)) {
    laggedIndex = k;
  } else {
    return MISSING_D;
  }

  double tNow = _t[laggedIndex];
  if (_reldt) {
    int targetIndex = i;
    get_target(i, copredictionMode, predictionSet, targetIndex);
    double tPred = (targetIndex >= 0) ? _t[targetIndex] : MISSING_D;

    if (tNow != MISSING_D && tPred != MISSING_D) {
      return dtWeight * (tPred - tNow);
    } else {
      return MISSING_D;
    }
  } else {
    double tNext = _t[prevLaggedIndex];
    if (tNext != MISSING_D && tNow != MISSING_D) {
      return dtWeight * (tNext - tNow);
    } else {
      return MISSING_D;
    }
  }
}

double ManifoldGenerator::get_target(int i, bool copredictionMode, bool predictionSet, int& targetIndex) const
{
  int panel = _panel_mode ? _panelIDs[i] : -1;
  bool use_co_x = copredictionMode && predictionSet;

  // What is the target of this point in the manifold?

  if (_p != 0) {
    // At what time does the prediction occur?
    int targetObsNum = _observation_number[targetIndex] + _p;
    int direction = _p > 0 ? 1 : -1;
    if (!find_observation_num(targetObsNum, targetIndex, direction, panel)) {
      targetIndex = -1;
    }
  }

  double target;

  if (targetIndex >= 0) {
    if (use_co_x) {
      target = _co_x[targetIndex];
    } else if (_xmap_mode) {
      target = _xmap[targetIndex];
    } else {
      target = _x[targetIndex];
    }
  } else {
    target = MISSING_D;
  }

  return target;
}

bool is_usable(double* point, int E_actual, bool allowMissing)
{
  if (allowMissing) {
    // If we are allowed to have missing values in the points, just
    // need to ensure that we don't allow a 100% missing point.
    for (int j = 0; j < E_actual; j++) {
      if (point[j] != MISSING_D) {
        return true;
      }
    }
    return false;

  } else {
    // Check that there are no missing values in the points.
    for (int j = 0; j < E_actual; j++) {
      if (point[j] == MISSING_D) {
        return false;
      }
    }
    return true;
  }
}

std::vector<bool> ManifoldGenerator::generate_usable(int maxE, bool copredictionMode) const
{
  const double USABLE_DTWEIGHT = 1.0;

  std::vector<bool> usable(_t.size());

  int sizeOfPoint = E_actual(maxE);
  auto point = std::make_unique<double[]>(sizeOfPoint);

  for (int i = 0; i < _t.size(); i++) {
    if (_t[i] == MISSING_D) {
      usable[i] = false;
      continue;
    }
    if (_panel_mode && _panelIDs[i] == MISSING_I) {
      usable[i] = false;
      continue;
    }
    fill_in_point(i, maxE, copredictionMode, copredictionMode, USABLE_DTWEIGHT, point.get());
    usable[i] = is_usable(point.get(), sizeOfPoint, _allow_missing);
  }

  return usable;
}

#ifdef JSON

void to_json(json& j, const ManifoldGenerator& g)
{
  j = json{ { "_dt", g._dt },
            { "_reldt", g._reldt },
            { "_dtWeight", g._dtWeight },
            { "_panel_mode", g._panel_mode },
            { "_xmap_mode", g._xmap_mode },
            { "_tau", g._tau },
            { "_p", g._p },
            { "_num_extras", g._num_extras },
            { "_num_extras_lagged", g._num_extras_lagged },
            { "_x", g._x },
            { "_xmap", g._xmap },
            { "_co_x", g._co_x },
            { "_t", g._t },
            { "_observation_number", g._observation_number },
            { "_extras", g._extras },
            { "_panelIDs", g._panelIDs } };
}

void from_json(const json& j, ManifoldGenerator& g)
{
  j.at("_dt").get_to(g._dt);
  j.at("_reldt").get_to(g._reldt);
  j.at("_dtWeight").get_to(g._dtWeight);
  j.at("_panel_mode").get_to(g._panel_mode);
  j.at("_xmap_mode").get_to(g._xmap_mode);
  j.at("_tau").get_to(g._tau);
  j.at("_p").get_to(g._p);
  j.at("_num_extras").get_to(g._num_extras);
  j.at("_num_extras_lagged").get_to(g._num_extras_lagged);
  j.at("_x").get_to(g._x);
  j.at("_xmap").get_to(g._xmap);
  j.at("_co_x").get_to(g._co_x);
  j.at("_t").get_to(g._t);
  j.at("_observation_number").get_to(g._observation_number);
  j.at("_extras").get_to(g._extras);
  j.at("_panelIDs").get_to(g._panelIDs);
}

#endif