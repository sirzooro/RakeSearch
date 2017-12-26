#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include "boinc_api.h"
#include "MovePairSearch.h"
#if defined(__i386__) || defined (__x86_64__)
#include <cpuid.h>
#include <stdio.h>
#endif
#ifdef __arm__
#include <sys/auxv.h>
#include <asm/hwcap.h>
#include <stdio.h>
#endif

using namespace std;

static const bool isDebug = false;

// Checking the file for existence
inline bool file_exists (const std::string& name)
{
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

// Computing
int Compute(string wu_filename, string result_filename)
{
  string localWorkunit;
  string localResult;
  string localCheckpoint; 
  string localTmpCheckpoint;

  string initStartFileName;
  string initResultFileName;
  string initCheckpointFileName;
  string initTmpCheckpointFileName;

  MovePairSearch search;

    // Checking the workunit file, the checkpoint file, the result file for existence
    localWorkunit = wu_filename; 
    localResult = result_filename; 
    localCheckpoint = "checkpoint.txt";
    localTmpCheckpoint = "tmp_checkpoint.txt";

    // Starting computing from a checkpoint
    if (file_exists(localCheckpoint))
    {
      // Checking the workunit file for existence
      if (file_exists(localWorkunit))
      {
        // Starting computing
        initStartFileName  = localWorkunit;
        initResultFileName = localResult;
        initCheckpointFileName = localCheckpoint;       
        initTmpCheckpointFileName = localTmpCheckpoint;
 
        if(isDebug) cout << "Start from checkpoint of workunit " << localWorkunit << endl;

        search.InitializeMoveSearch(initStartFileName, initResultFileName, 
                      initCheckpointFileName, initTmpCheckpointFileName);
        search.StartMoveSearch();
      }
      else
      {
        cerr << "Error: detected a checkpoint file " << localCheckpoint;
        cerr << " without workunit file!" << endl;
        return -1;
      }
    }

    // Starting computing from a workunit file
    if (!file_exists(localCheckpoint) && file_exists(localWorkunit))
    {
      // Starting computing from a workunit file, 
      // present without the files
      // of a checkpoint and a result
      initStartFileName  = localWorkunit;
      initResultFileName = localResult;
      initCheckpointFileName = localCheckpoint;      
      initTmpCheckpointFileName = localTmpCheckpoint;

      if(isDebug) cout << "Start from workunit file " << localWorkunit << endl;
      search.InitializeMoveSearch(initStartFileName, initResultFileName, 
                      initCheckpointFileName, initTmpCheckpointFileName);
      search.StartMoveSearch();
    }

  return 0;
}


int main(int argumentsCount, char* argumentsValues[])
{
  string wu_filename = "workunit.txt";
  string result_filename = "result.txt";
  string resolved_in_name;  // Variables for working with logical
  string resolved_out_name; // and physical filenames in BOINC

  int retval;

  boinc_init(); // Initialize BOINC API for single-threaded application 
  // Set minimal number of seconds between writing checkpoints
  boinc_set_min_checkpoint_period(60); 
  
  // Convert logical filename into physical filename.
  // We do it at the high level, passing already resolved filenames further.
  retval = boinc_resolve_filename_s(wu_filename.c_str(), resolved_in_name);
  if (retval) 
  { 
    if(isDebug) cerr << "can't resolve IN filename!" << endl; 
    boinc_finish(retval); return 0;
  }
  retval = boinc_resolve_filename_s(result_filename.c_str(), resolved_out_name);
  if (retval) 
  { 
    if(isDebug) cerr << "can't resolve OUT filename" << endl; 
    boinc_finish(retval); return 0;
  }

  boinc_fraction_done(0.0); // Tell the BOINC client the fraction of workunit completion
    retval = Compute(resolved_in_name, resolved_out_name); // Launch computing
  boinc_fraction_done(1.0); // Tell the BOINC client the fraction of workunit completion

  // Tell the BOINC client the status of workunit completion (does not return)
  boinc_finish(retval); 
                   // In case we need to show a message to the user, use function
                   // boinc_finish_message(int status, const char* msg, bool is_notice); 
                   // If is_notice is true, the message will be shown as a notice 
                   // in the GUI (works with 7.5+ clients; for others, no message 
                   // will be shown). 
  return 0;
}

// Do not check x86 CPU features if SSE2 is not required. Very old CPUs may not support CPUID.
// Code below uses C stdio because C++ iostreams are not initialized yet and app would crash.
// This function is called before main() starts.
#if (defined(__i386__) || defined (__x86_64__)) && defined(__SSE2__)
__attribute__((constructor(101), target("no-sse")))
void CheckCpuid()
{
  unsigned int a, b, c, d;

  if (!__get_cpuid(1, &a, &b, &c, &d))
  {
    fprintf(stderr, "CPUID instruction is not supported by your CPU!\n");
    exit(1);
  }

  if (0 == (d & bit_SSE2))
  {
    fprintf(stderr, "SSE2 instructions are not supported by your CPU!\n");
    exit(1);
  }

#ifdef __AVX__
  if (0 == (c & bit_AVX))
  {
    fprintf(stderr, "AVX instructions are not supported by your CPU!\n");
    exit(1);
  }

  // AVX also needs OS support, check for it
  if (0 == (c & bit_OSXSAVE))
  {
    fprintf(stderr, "OSXSAVE instructions are not supported by your CPU!\n");
    exit(1);
  }

  unsigned int eax, edx;
  unsigned int ecx = 0; // _XCR_XFEATURE_ENABLED_MASK
  __asm__ ("xgetbv" : "=a" (eax), "=d" (edx) : "c" (ecx));
  if (0x6 != (eax & 0x6)) // XSTATE_SSE | XSTATE_YMM
  {
    fprintf(stderr, "AVX instructions are not supported by your OS!\n");
    exit(1);
  }
#endif

#ifdef __AVX2__
  //if (!__get_cpuid(7, &a, &b, &c, &d)) {
  if (__get_cpuid_max(0, 0) < 7)
  {
    fprintf(stderr, "CPUID level 7 is not supported by your CPU!\n");
    exit(1);
  }

  __cpuid_count(7, 0, a, b, c, d);

  if (0 == (b & bit_AVX2))
  {
    fprintf(stderr, "AVX2 instructions are not supported by your CPU!\n");
    exit(1);
  }
#endif

#ifdef __BMI2__
#ifndef __AVX2__
#error AVX2 is not enabled!
#endif
  const unsigned int bmibits = bit_BMI | bit_BMI2;
  if (bmibits != (b & bmibits))
  {
    fprintf(stderr, "BMI1&2 instructions are not supported by your CPU!\n");
    exit(1);
  }
#endif

#ifdef __AVX512F__
  // During comppilation following AVX512 instruction sets are enabled: F, CD, VL, BW, DQ
  // Code directly uses F and VL. Compiler may emit instruction from other ones, so check them all

  // We need data loaded during AVX and AVX2 checks, make sure we have it
#if !defined(__AVX__) || !defined(__AVX2__)
#error AVX or AVX2 is not enabled!
#endif

  const unsigned int avx512bits = bit_AVX512F | bit_AVX512DQ | bit_AVX512CD | bit_AVX512BW | bit_AVX512VL;

  if (avx512bits != (b & avx512bits))
  {
    fprintf(stderr, "AVX512 instructions are not supported by your CPU!\n");
    exit(1);
  }

  // AVX512 also needs OS support, check for it
  if (0xe6 != (eax & 0xe6)) // XSTATE_SSE | XSTATE_YMM | XSTATE_OPMASK | XSTATE_ZMM | XSTATE_HI_ZMM
  {
    fprintf(stderr, "AVX512 instructions are not supported by your OS!\n");
    exit(1);
  }
#endif
}
#endif

// Check if ARM CPU supports NEON instructions if required. For AARCH64 they are always present,
// so no need to check them here.
#if defined(__arm__) && defined(__ARM_NEON)
__attribute__((constructor(101)))
void CheckArmHwcap()
{
  if (HWCAP_NEON != (getauxval(AT_HWCAP) & HWCAP_NEON))
  {
    fprintf(stderr, "NEON instructions are not supported by your CPU!\n");
    exit(1);
  }
}
#endif
