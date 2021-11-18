#pragma once

template<typename T>
void cuLPDistances(char * const valids, T * const distances,
                   const int numPredictions, const bool isDistanceMeanAbsoluteError, const bool isPanelMode,
                   const double idw, const double missingDistance,
                   const int eacts, const int numLibraryPoints,
                   const T* mData, const int* mPanelIds,
                   const T* mpData, const int* mpPanelIds,
                   const char* metricOptions, const cudaStream_t);
