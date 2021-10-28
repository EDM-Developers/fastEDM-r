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
  if (!_dt) {
    // In normal situations (non-dt)
    double unit = calculate_time_increment();
    double minT = *std::min_element(_t.begin(), _t.end());

    // Create a time index which is a discrete count of the number of 'unit' time units.
    for (int i = 0; i < _t.size(); i++) {
      if (_t[i] != MISSING_D) {
        _observation_number.push_back(std::round((_t[i] - minT) / unit));
      } else {
        _observation_number.push_back(-1);
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
        _observation_number.push_back(-1);
      }
    }
  }
}

bool ManifoldGenerator::find_observation_num(int target, int& k, int direction, int panel) const
{
  // Loop either forward or back until we find the right index or give up.
  while (k >= 0 && k < _observation_number.size()) {
    // If in panel mode, make sure we don't wander over a panel boundary.
    if (_panel_mode) {
      if (panel != _panel_ids[k]) {
        return false;
      }
    }

    // Skip over garbage rows which don't have a time recorded.
    if (_observation_number[k] < 0) {
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

std::vector<int> ManifoldGenerator::get_lagged_indices(int startIndex, int E, int panel) const
{

  std::vector<int> laggedIndices(E);
  std::fill_n(laggedIndices.begin(), E, -1);

  // For obs i, which indices correspond to looking back 0, tau, ..., (E-1)*tau observations.
  laggedIndices[0] = startIndex;
  int pointStartObsNum = _observation_number[startIndex];

  // Start by going back one index
  int k = startIndex - 1;

  for (int j = 1; j < E; j++) {
    // Find the discrete time we're searching for.
    int targetObsNum = pointStartObsNum - j * _tau;

    if (find_observation_num(targetObsNum, k, -1, panel)) {
      laggedIndices[j] = k;
    }
  }

  return laggedIndices;
}

#if defined(WITH_ARRAYFIRE)
ManifoldOnGPU Manifold::toGPU(const bool useFloat) const
{
  using af::array;

  if (useFloat) {
    return ManifoldOnGPU{ array(_E_actual, _nobs, _flat.get()).as(f32),
                          (_y.size() > 0 ? array(_nobs, _y.data()) : array()).as(f32),
                          (_panel_ids.size() > 0 ? array(_nobs, _panel_ids.data()) : array()),
                          _nobs,
                          _E_x,
                          _E_dt,
                          _E_extras,
                          _E_lagged_extras,
                          _E_actual,
                          MISSING_F };
  } else {
    return ManifoldOnGPU{ array(_E_actual, _nobs, _flat.get()),
                          (_y.size() > 0 ? array(_nobs, _y.data()) : array()),
                          (_panel_ids.size() > 0 ? array(_nobs, _panel_ids.data()) : array()),
                          _nobs,
                          _E_x,
                          _E_dt,
                          _E_extras,
                          _E_lagged_extras,
                          _E_actual,
                          MISSING_D };
  }
}
#endif

Manifold ManifoldGenerator::create_manifold(int E, const std::vector<bool>& filter, bool copredict, bool prediction,
                                            double dtWeight, bool skipMissing) const
{
  bool takeEveryPoint = filter.size() == 0;

  int nobs = 0;
  std::vector<int> pointNumToStartIndex;
  for (int i = 0; i < _t.size(); i++) {
    if (takeEveryPoint || filter[i]) {
      pointNumToStartIndex.push_back(i);
      nobs += 1;
    }
  }

  auto flat = std::make_unique<double[]>(nobs * E_actual(E));

  std::vector<double> y;
  std::vector<int> panelIDs;

  // Fill in the manifold row-by-row (point-by-point)
  int M_i = 0;
  double target;

  for (int i = 0; i < nobs; i++) {
    double* point = &(flat[M_i * E_actual(E)]);
    fill_in_point(pointNumToStartIndex[i], E, copredict, prediction, dtWeight, point, target);

    // Erase this point if we don't want missing values in the resulting manifold
    if (skipMissing) {
      bool foundMissing = false;
      for (int j = 0; j < E_actual(E); j++) {
        if (point[j] == MISSING_D) {
          foundMissing = true;
          break;
        }
      }

      if (foundMissing) {
        continue;
      }
    }

    y.push_back(target);
    if (_panel_mode) {
      panelIDs.push_back(_panel_ids[pointNumToStartIndex[i]]);
    }

    M_i += 1;
  }

  nobs = M_i;

  return { flat, y, panelIDs, nobs, E, E_dt(E), E_extras(E), E * numExtrasLagged(), E_actual(E) };
}

void ManifoldGenerator::fill_in_point(int i, int E, bool copredict, bool prediction, double dtWeight, double* point,
                                      double& target) const
{
  int panel = _panel_mode ? _panel_ids[i] : -1;
  bool use_co_x = copredict && prediction;

  std::vector<int> laggedIndices = get_lagged_indices(i, E, panel);

  auto lookup_vec = [&laggedIndices](const std::vector<double>& vec, int j) {
    if (laggedIndices[j] < 0) {
      return MISSING_D;
    } else {
      return vec[laggedIndices[j]];
    }
  };

  // What is the target of this point in the manifold?
  int targetIndex = i;

  if (_p != 0) {
    // At what time does the prediction occur?
    int targetObsNum = _observation_number[targetIndex] + _p;
    int direction = _p > 0 ? 1 : -1;
    if (!find_observation_num(targetObsNum, targetIndex, direction, panel)) {
      targetIndex = -1;
    }
  }

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

  // Fill in the lagged embedding of x (or co_x) in the first columns
  for (int j = 0; j < E; j++) {
    if (use_co_x) {
      point[j] = lookup_vec(_co_x, j);
    } else {
      point[j] = lookup_vec(_x, j);
    }
  }

  // Put the lagged embedding of dt in the next columns
  if (_dt) {
    double tPred = (targetIndex >= 0) ? _t[targetIndex] : MISSING_D;

    for (int j = 0; j < E_dt(E); j++) {
      double tNow = lookup_vec(_t, j);
      if (j == 0 || _reldt) {
        if (tNow != MISSING_D && tPred != MISSING_D) {
          point[E + j] = dtWeight * (tPred - tNow);
        } else {
          point[E + j] = MISSING_D;
        }
      } else {
        double tNext = lookup_vec(_t, j - 1);
        if (tNext != MISSING_D && tNow != MISSING_D) {
          point[E + j] = dtWeight * (tNext - tNow);
        } else {
          point[E + j] = MISSING_D;
        }
      }
    }
  }

  // Finally put the extras in the last columns
  int offset = 0;
  for (int k = 0; k < _num_extras; k++) {
    int numLags = (k < _num_extras_lagged) ? E : 1;
    for (int j = 0; j < numLags; j++) {
      point[E + E_dt(E) + offset + j] = lookup_vec(_extras[k], j);
    }
    offset += numLags;
  }
}

bool is_usable(double* point, double target, int E_actual, bool allowMissing, bool targetRequired)
{
  if (targetRequired && target == MISSING_D) {
    return false;
  }

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

std::vector<bool> ManifoldGenerator::generate_usable(int maxE, bool coprediction) const
{
  const double USABLE_DTWEIGHT = 1.0;

  bool targetRequired = true;

  std::vector<bool> usable(_t.size());

  int E = E_actual(maxE);
  auto point = std::make_unique<double[]>(E);
  double target;

  for (int i = 0; i < _t.size(); i++) {
    fill_in_point(i, maxE, coprediction, coprediction, USABLE_DTWEIGHT, point.get(), target);
    usable[i] = is_usable(point.get(), target, E, _allow_missing, targetRequired);
  }

  return usable;
}

void to_json(json& j, const ManifoldGenerator& g)
{
  j = json{ { "_dt", g._dt },
            { "_dt0", g._dt0 },
            { "_reldt", g._reldt },
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
            { "_panel_ids", g._panel_ids } };
}

void from_json(const json& j, ManifoldGenerator& g)
{
  j.at("_dt").get_to(g._dt);
  j.at("_dt0").get_to(g._dt0);
  j.at("_reldt").get_to(g._reldt);
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
  j.at("_panel_ids").get_to(g._panel_ids);
}
