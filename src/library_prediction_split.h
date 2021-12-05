#pragma once

#include <numeric>
#include <vector>

#include "mersennetwister.h"
#include "stats.h"

enum class Set
{
  Neither,
  Library,
  Prediction,
  Both
};

class LibraryPredictionSetSplitter
{
private:
  bool _explore, _full, _shuffle;
  int _crossfold, _numObsUsable;
  std::vector<bool> _usable;
  std::vector<bool> _libraryRows, _predictionRows;
  std::vector<int> _crossfoldGroup;
  MtRng64 _rng;

public:
  LibraryPredictionSetSplitter(bool explore, bool full, bool shuffle, int crossfold, std::vector<bool> usable,
                               const std::string& rngState = "")
    : _explore(explore)
    , _full(full)
    , _shuffle(shuffle)
    , _crossfold(crossfold)
    , _usable(usable)
  {
    if (!rngState.empty() && shuffle) {
      // Sync the local random number generator with Stata's
      set_rng_state(rngState);
    } else {
      _rng.init((unsigned long long)0);
    }

    _numObsUsable = std::accumulate(usable.begin(), usable.end(), 0);

    if (crossfold > 0) {
      setup_crossfold_groups();
    }
  }

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

  void setup_crossfold_groups()
  {
    std::vector<int> uRank;

    if (_shuffle) {
      std::vector<double> u;
      for (int i = 0; i < _numObsUsable; i++) {
        u.push_back(_rng.getReal2());
      }
      uRank = rank(u);
    }

    int sizeOfEachFold = std::round(((float)_numObsUsable) / _crossfold);
    int obsNum = 0;

    _crossfoldGroup = std::vector<int>(_usable.size());

    for (int i = 0; i < _usable.size(); i++) {
      if (_usable[i]) {
        if (_shuffle) {
          _crossfoldGroup[i] = uRank[obsNum] % _crossfold;
        } else {
          _crossfoldGroup[i] = obsNum / sizeOfEachFold;
        }
        obsNum += 1;
      } else {
        _crossfoldGroup[i] = -1;
      }
    }
  }

  // Assuming this is called in explore mode
  int next_library_size(int crossfoldIter) const
  {
    int librarySize = 0;
    if (_crossfold > 0) {
      for (int obsNum = 0; obsNum < _numObsUsable; obsNum++) {
        if ((obsNum + 1) % _crossfold != (crossfoldIter - 1)) {
          librarySize += 1;
        }
      }
      return librarySize;
    } else if (_full) {
      return _numObsUsable;
    } else {
      return _numObsUsable / 2;
    }
  }

  std::vector<bool> libraryRows() const { return _libraryRows; }
  std::vector<bool> predictionRows() const { return _predictionRows; }
  std::vector<int> crossfold_groups() const { return _crossfoldGroup; }

  std::vector<Set> setMemberships() const
  {

    std::vector<Set> m;

    for (int i = 0; i < _libraryRows.size(); i++) {
      if (_libraryRows[i] && _predictionRows[i]) {
        m.push_back(Set::Both);
      } else if (_libraryRows[i]) {
        m.push_back(Set::Library);
      } else if (_predictionRows[i]) {
        m.push_back(Set::Prediction);
      } else {
        m.push_back(Set::Neither);
      }
    }

    return m;
  }

  void update_library_prediction_split(int library = -1, int crossfoldIter = -1)
  {
    if (_explore && _full) {
      _libraryRows = _usable;
      _predictionRows = _usable;
    } else if (_explore && _crossfold > 0) {
      crossfold_split(crossfoldIter);
    } else if (_explore) {
      half_library_prediction_split();
    } else {
      fixed_size_library(library);
    }
  }

  void crossfold_split(int crossfoldIter)
  {
    _libraryRows = std::vector<bool>(_usable.size());
    _predictionRows = std::vector<bool>(_usable.size());

    for (int i = 0; i < _libraryRows.size(); i++) {
      if (_usable[i]) {
        if (_crossfoldGroup[i] == (crossfoldIter - 1)) {
          _libraryRows[i] = false;
          _predictionRows[i] = true;
        } else {
          _libraryRows[i] = true;
          _predictionRows[i] = false;
        }
      } else {
        _libraryRows[i] = false;
        _predictionRows[i] = false;
      }
    }
  }

  void half_library_prediction_split()
  {
    _libraryRows = std::vector<bool>(_usable.size());
    _predictionRows = std::vector<bool>(_usable.size());

    if (_shuffle) {
      std::vector<double> u;
      for (int i = 0; i < _numObsUsable; i++) {
        u.push_back(_rng.getReal2());
      }

      double med = median(u);

      int obsNum = 0;
      for (int i = 0; i < _libraryRows.size(); i++) {
        if (_usable[i]) {
          if (u[obsNum] < med) {
            _libraryRows[i] = true;
            _predictionRows[i] = false;
          } else {
            _libraryRows[i] = false;
            _predictionRows[i] = true;
          }
          obsNum += 1;
        } else {
          _libraryRows[i] = false;
          _predictionRows[i] = false;
        }
      }
    } else {
      int librarySize = _numObsUsable / 2;

      int obsNum = 0;
      for (int i = 0; i < _libraryRows.size(); i++) {
        if (_usable[i]) {
          if (obsNum < librarySize) {
            _libraryRows[i] = true;
            _predictionRows[i] = false;
          } else {
            _libraryRows[i] = false;
            _predictionRows[i] = true;
          }
          obsNum += 1;
        } else {
          _libraryRows[i] = false;
          _predictionRows[i] = false;
        }
      }
    }
  }

  void fixed_size_library(int library)
  {
    _libraryRows = std::vector<bool>(_usable.size());
    _predictionRows = _usable;

    if (_shuffle) {
      std::vector<double> u;
      for (int i = 0; i < _numObsUsable; i++) {
        u.push_back(_rng.getReal2());
      }

      double uCutoff = 1.0;
      if (library < u.size()) {
        std::vector<double> uCopy(u);
        const auto uCutoffIt = uCopy.begin() + library;
        std::nth_element(uCopy.begin(), uCutoffIt, uCopy.end());
        uCutoff = *uCutoffIt;
      }

      int obsNum = 0;
      for (int i = 0; i < _libraryRows.size(); i++) {
        if (_usable[i]) {
          _predictionRows[i] = true;
          if (u[obsNum] < uCutoff) {
            _libraryRows[i] = true;
          } else {
            _libraryRows[i] = false;
          }
          obsNum += 1;
        } else {
          _libraryRows[i] = false;
        }
      }
    } else {
      int obsNum = 0;
      for (int i = 0; i < _libraryRows.size(); i++) {
        if (_usable[i]) {
          if (obsNum < library) {
            _libraryRows[i] = true;
          } else {
            _libraryRows[i] = false;
          }
          obsNum += 1;
        } else {
          _libraryRows[i] = false;
        }
      }
    }

    int numInLibrary = 0;
    for (int i = 0; i < _libraryRows.size(); i++) {
      if (_libraryRows[i]) {
        numInLibrary += 1;
      }
    }
  }
};
