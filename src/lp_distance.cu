#include "lp_distance.cuh"

#define divup(a, b) (((a) + (b)-1) / (b))

template<typename T>
__device__ constexpr T getMissingConstant()
{
  return 1.0e+100; // double
}

template<>
__device__ float getMissingConstant()
{
  return 1.0e+30;
}

template<typename T>
__device__ inline T getValue(T value)
{
  return value;
}

template<>
__device__ inline float getValue(float value)
{
  return isinf(value) ? getMissingConstant<float>() : value;
}

template<typename T, int BLOCK_DIM_X, int BLOCK_DIM_Y>
__device__ T reduceEacts(int tid, T* smem)
{
  constexpr unsigned int BLK_SIZE = BLOCK_DIM_X * BLOCK_DIM_Y;

  if (BLK_SIZE > 128) {
    if (tid < 128) { smem[tid] += smem[tid + 128]; }
    __syncthreads();
  }
  if (BLK_SIZE > 64) {
    if (tid < 64) { smem[tid] += smem[tid + 64]; }
    __syncthreads();
  }
  if (BLK_SIZE > 32) {
    if (tid < 32) { smem[tid] += smem[tid + 32]; }
    __syncthreads();
  }

  T retVal = smem[ tid ];
  __syncthreads();

  return retVal;
}

template<typename T, bool isDMAE, int BLOCK_DIM_X, int BLOCK_DIM_Y>
__global__
void lpDistances(char * const valids, T * const distances,
                 const int numPredictions, const bool isPanelMode,
                 const double idw, const double missingDistance,
                 const int eacts, const int numLibraryPoints,
                 const T* mData, const int* mPanelIds,
                 const T* mpData, const int* mpPanelIds,
                 const char* mopts)
{
  const T MISSING = getMissingConstant<T>();

  const int p = blockIdx.y; //nth prediction

  if (p < numPredictions)
  {
    __shared__ T dists[BLOCK_DIM_X * BLOCK_DIM_Y];
    __shared__ bool markers[BLOCK_DIM_X * BLOCK_DIM_Y];

    const bool isZero = (missingDistance == 0);
    const T* predsMp  = mpData + p * eacts;

    const int tid = BLOCK_DIM_X * threadIdx.y + threadIdx.x;
    const int nob = BLOCK_DIM_X * blockIdx.x + threadIdx.x;

    if (nob < numLibraryPoints)
    {
      const T* predsM   = mData + nob * eacts;
      bool anyEAmissing = false;
      T dist_i          = T(0);

      if ( threadIdx.y == 0 && isPanelMode && idw > 0 ) {
        dist_i += (idw * (mPanelIds[nob] != mpPanelIds[p]));
      }
      for (int e = threadIdx.y; e < eacts; e += BLOCK_DIM_Y)
      {
        T M_ij    = getValue(predsM[e]);
        T Mp_ij   = getValue(predsMp[e]);
        bool mopt = mopts[e];
        bool msng = (M_ij == MISSING || Mp_ij == MISSING);
        T diffM   = M_ij - Mp_ij;
        T compM   = M_ij != Mp_ij;
        T distM   = mopt * diffM + (1 - mopt) * compM;
        T dist_ij = msng * (1 - isZero) * missingDistance + (1 - msng) * distM;

        if (isDMAE) {
          dist_i += abs(dist_ij) / eacts;
        } else {
          dist_i += dist_ij * dist_ij;
        }
        anyEAmissing = (anyEAmissing || msng);
      }
      __syncthreads();

      dists[ tid ] = dist_i;
      __syncthreads();
      dist_i = reduceEacts<T, BLOCK_DIM_X, BLOCK_DIM_Y>(tid, dists);

      markers[ tid ] = anyEAmissing;
      __syncthreads();
      anyEAmissing = reduceEacts<bool, BLOCK_DIM_X, BLOCK_DIM_Y>(tid, markers) > 0;

      // Only first warp writes numLibraryPoints results
      if (tid < BLOCK_DIM_X) {
        anyEAmissing = anyEAmissing && isZero;

        dist_i = anyEAmissing * MISSING + (1 - anyEAmissing) * dist_i;

        bool isValid = dist_i != 0 && dist_i != MISSING;

        dist_i = (isDMAE * dist_i + (1 - isDMAE) * sqrt(dist_i));

        valids[nob + p * numLibraryPoints] = (char)isValid;
        distances[nob + p * numLibraryPoints] = dist_i;
      }
    }
  }
}

template<typename T, int BLOCK_DIM_Y>
void lpDistances(char * const valids, T * const distances,
                 const int numPredictions, const bool isDMAE, const bool isPanelMode,
                 const double idw, const double missingDistance,
                 const int eacts, const int numLibraryPoints,
                 const T* mData, const int* mPanelIds,
                 const T* mpData, const int* mpPanelIds,
                 const char* mopts, const cudaStream_t stream,
                 const dim3 blocks, const dim3 threads)
{
  if (isDMAE) {
    lpDistances<T, true, 32, BLOCK_DIM_Y> <<<blocks, threads, 0, stream>>>(
            valids, distances, numPredictions, isPanelMode, idw, missingDistance, eacts, numLibraryPoints,
            mData, mPanelIds, mpData, mpPanelIds, mopts);
  } else {
    lpDistances<T, false, 32, BLOCK_DIM_Y> <<<blocks, threads, 0, stream>>>(
            valids, distances, numPredictions, isPanelMode, idw, missingDistance, eacts, numLibraryPoints,
            mData, mPanelIds, mpData, mpPanelIds, mopts);
  }
}

inline unsigned int powerOf2LE(unsigned int value)
{
    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;

    return value ^ (value >> 1);
}

template<typename T>
void cuLPDistances(char * const valids, T * const distances,
                   const int numPredictions, const bool isDMAE, const bool isPanelMode,
                   const double idw, const double missingDistance,
                   const int eacts, const int numLibraryPoints,
                   const T* mData, const int* mPanelIds,
                   const T* mpData, const int* mpPanelIds,
                   const char* mopts, const cudaStream_t stream)
{
  dim3 threads(32, powerOf2LE(eacts));

  threads.y = (threads.y > 8 ? 8 : threads.y);

  dim3 blocks(divup(numLibraryPoints, threads.x), numPredictions);

  switch(threads.y) {
    case 8:
      lpDistances<T, 8>(valids, distances, numPredictions, isDMAE, isPanelMode, idw, missingDistance,
              eacts, numLibraryPoints, mData, mPanelIds, mpData, mpPanelIds, mopts, stream, blocks, threads);
      break;
    case 4:
      lpDistances<T, 4>(valids, distances, numPredictions, isDMAE, isPanelMode, idw, missingDistance,
              eacts, numLibraryPoints, mData, mPanelIds, mpData, mpPanelIds, mopts, stream, blocks, threads);
      break;
    default:
      lpDistances<T, 2>(valids, distances, numPredictions, isDMAE, isPanelMode, idw, missingDistance,
              eacts, numLibraryPoints, mData, mPanelIds, mpData, mpPanelIds, mopts, stream, blocks, threads);
      break;
  }
}

#define INSTANTIATE(T)                                                                      \
template void cuLPDistances(char * const, T * const, const int, const bool, const bool,     \
        const double, const double, const int, const int, const T*, const int*, const T*,   \
        const int*, const char*, const cudaStream_t);

INSTANTIATE(float)
INSTANTIATE(double)
