#include "cpu.h"

#ifdef _MSC_VER
#include <windows.h>
#elif defined __APPLE__
#include <sys/sysctl.h>
#endif

size_t num_logical_cores()
{
#ifdef _MSC_VER
  // Inspired by https://github.com/ninja-build/ninja/blob/master/src/util.cc
  return GetActiveProcessorCount(ALL_PROCESSOR_GROUPS);
#else
  return std::thread::hardware_concurrency();
#endif
}

// Adapted from https://github.com/giampaolo/psutil/blob/master/psutil/arch/windows/cpu.c
// and https://github.com/giampaolo/psutil/blob/master/psutil/_psutil_osx.c

/*
 * Return the number of physical CPU cores (hyper-thread CPUs count
 * is excluded).
 */
size_t num_physical_cores()
{
#ifdef _MSC_VER
  DWORD rc;
  PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX buffer = NULL;
  PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX ptr = NULL;
  DWORD length = 0;
  DWORD offset = 0;
  DWORD ncpus = 0;
  DWORD prev_processor_info_size = 0;

  // GetLogicalProcessorInformationEx() is available from Windows 7
  // onward. Differently from GetLogicalProcessorInformation()
  // it supports process groups, meaning this is able to report more
  // than 64 CPUs. See:
  // https://bugs.python.org/issue33166
  if (GetLogicalProcessorInformationEx == NULL) {
    return 0;
  }

  while (1) {
    rc = GetLogicalProcessorInformationEx(RelationAll, buffer, &length);
    if (rc == FALSE) {
      if (GetLastError() == ERROR_INSUFFICIENT_BUFFER) {
        if (buffer) {
          free(buffer);
        }
        buffer = (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX)malloc(length);
        if (NULL == buffer) {
          return 0;
        }
      } else {
        if (buffer) {
          free(buffer);
        }
        return 0;
      }
    } else {
      break;
    }
  }

  ptr = buffer;
  while (offset < length) {
    // Advance ptr by the size of the previous
    // SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX struct.
    ptr = (SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*)(((char*)ptr) + prev_processor_info_size);

    if (ptr->Relationship == RelationProcessorCore) {
      ncpus += 1;
    }

    // When offset == length, we've reached the last processor
    // info struct in the buffer.
    offset += ptr->Size;
    prev_processor_info_size = ptr->Size;
  }

  if (buffer != NULL) {
    free(buffer);
  }

  if (ncpus != 0) {
    return ncpus;
  } else {
    return 0;
  }

#elif defined __APPLE__

  int num;
  size_t size = sizeof(int);

  if (sysctlbyname("hw.physicalcpu", &num, &size, NULL, 0))
    return 0;
  return num;

#else
  // TODO: Is there a proper 'linux' way to get # physical cores?
  return std::thread::hardware_concurrency();
#endif
}

// Adapted from https://chrisgreendevelopmentblog.wordpress.com/2017/08/29/thread-pools-and-windows-processor-groups/
void distribute_threads(std::vector<std::thread>& threads)
{
#ifdef _MSC_VER

  DWORD htPerCore = (DWORD)(num_logical_cores() / num_physical_cores());

  int nNumGroups = GetActiveProcessorGroupCount();
  if (nNumGroups > 1) {
    WORD nCurGroup = 0;
    DWORD nNumRemaining = GetMaximumProcessorCount(nCurGroup);
    for (int i = 0; i < threads.size(); i++) {
      auto hndl = threads[i].native_handle();
      GROUP_AFFINITY oldaffinity;
      if (GetThreadGroupAffinity(hndl, &oldaffinity)) {
        GROUP_AFFINITY affinity;
        affinity = oldaffinity;
        if (affinity.Group != nCurGroup) {
          affinity.Group = nCurGroup;
          SetThreadGroupAffinity(hndl, &affinity, nullptr);
          nNumRemaining -= htPerCore;
          if (nNumRemaining <= 0) {
            nCurGroup = (nCurGroup + 1) % nNumGroups;
            nNumRemaining = GetMaximumProcessorCount(nCurGroup);
          }
        }
      }
    }
  }

#else
// TODO: https://eli.thegreenplace.net/2016/c11-threads-affinity-and-hyperthreading/
#endif
}