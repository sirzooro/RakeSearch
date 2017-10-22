# include <iostream>
# include <fstream>
# include <string>
# include <sys/stat.h>
# include "boinc_api.h"

# include "MovePairSearch.h"

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
