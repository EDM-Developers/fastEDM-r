#pragma once

#include <vector>

#include "mersennetwister.h"
#include "stats.h"

class TrainPredictSplitter
{
private:
  bool _explore, _full;
  int _crossfold, _numObsUsable;
  std::vector<bool> _usable;
  std::vector<bool> _trainingRows, _predictionRows;
  std::vector<int> _crossfoldURank;
  MtRng64 _rng;

public:
  TrainPredictSplitter(bool explore, bool full, int crossfold, std::vector<bool> usable, const std::string& rngState)
    : _explore(explore)
    , _full(full)
    , _crossfold(crossfold)
    , _usable(usable)
  {
    if (!rngState.empty()) {
      // Sync the local random number generator with Stata's
      set_rng_state(rngState);
    } else {
      _rng.init((unsigned long long)0);
    }

    _numObsUsable = std::accumulate(usable.begin(), usable.end(), 0);

    if (crossfold > 0) {
      std::vector<double> u;

      for (int i = 0; i < _numObsUsable; i++) {
        u.push_back(_rng.getReal2());
      }

      _crossfoldURank = rank(u);
    }
  }

  static bool requiresRandomNumbers(int crossfold, bool full) { return crossfold > 0 || !full; }

  void set_rng_state(const std::string& rngState)
  {
    unsigned long long state[312];

    // Set up the rng at the beginning on this batch (given by the 'state' array)
    for (int i = 0; i < 312; i++) {
      state[i] = std::stoull(rngState.substr(3 + i * 16, 16), nullptr, 16);
      _rng.state_[i] = state[i];
    }

    _rng.left_ = 312;
    _rng.next_ = _rng.state_;

    // Burn all the rv's which are already used
    std::string countStr = rngState.substr(3 + 312 * 16 + 4, 8);
    long long numUsed = std::stoull(countStr, nullptr, 16);

    for (int i = 0; i < numUsed; i++) {
      _rng.getReal2();
    }
  }

  // Assuming this is called in explore mode
  int next_training_size(int crossfoldIter) const
  {
    int trainSize = 0;
    if (_crossfold > 0) {
      for (int obsNum = 0; obsNum < _numObsUsable; obsNum++) {
        if ((obsNum + 1) % _crossfold != (crossfoldIter - 1)) {
          trainSize += 1;
        }
      }
      return trainSize;
    } else if (_full) {
      return _numObsUsable;
    } else {
      return _numObsUsable / 2;
    }
  }

  void update_train_predict_split(int library, int crossfoldIter)
  {
    if (_explore && _full) {
      _trainingRows = _usable;
      _predictionRows = _usable;
      return;
    }

    _trainingRows = std::vector<bool>(_usable.size());
    _predictionRows = std::vector<bool>(_usable.size());

    if (_explore && _crossfold > 0) {
      int obsNum = 0;
      for (int i = 0; i < _trainingRows.size(); i++) {
        if (_usable[i]) {
          if (_crossfoldURank[obsNum] % _crossfold == (crossfoldIter - 1)) {
            _trainingRows[i] = false;
            _predictionRows[i] = true;
          } else {
            _trainingRows[i] = true;
            _predictionRows[i] = false;
          }
          obsNum += 1;
        } else {
          _trainingRows[i] = false;
          _predictionRows[i] = false;
        }
      }

      return;
    }

    std::vector<double> u;

    for (int i = 0; i < _numObsUsable; i++) {
      u.push_back(_rng.getReal2());
    }

    if (_explore) {
      double med = median(u);

      int obsNum = 0;
      for (int i = 0; i < _trainingRows.size(); i++) {
        if (_usable[i]) {
          if (u[obsNum] < med) {
            _trainingRows[i] = true;
            _predictionRows[i] = false;
          } else {
            _trainingRows[i] = false;
            _predictionRows[i] = true;
          }
          obsNum += 1;
        } else {
          _trainingRows[i] = false;
          _predictionRows[i] = false;
        }
      }
    } else {
      double uCutoff = 1.0;
      if (library < u.size()) {
        std::vector<double> uCopy(u);
        const auto uCutoffIt = uCopy.begin() + library;
        std::nth_element(uCopy.begin(), uCutoffIt, uCopy.end());
        uCutoff = *uCutoffIt;
      }

      int obsNum = 0;
      for (int i = 0; i < _trainingRows.size(); i++) {
        if (_usable[i]) {
          _predictionRows[i] = true;
          if (u[obsNum] < uCutoff) {
            _trainingRows[i] = true;
          } else {
            _trainingRows[i] = false;
          }
          obsNum += 1;
        } else {
          _trainingRows[i] = false;
          _predictionRows[i] = false;
        }
      }
    }
  }

  std::vector<bool> trainingRows() const { return _trainingRows; }
  std::vector<bool> predictionRows() const { return _predictionRows; }
};
