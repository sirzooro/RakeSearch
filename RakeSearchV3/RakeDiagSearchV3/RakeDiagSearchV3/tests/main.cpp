#include "TestRakeSearch.h"
#include <iostream>

int main(int argc, char* argv[])
{
    TestRakeSearch search;
    
    TestNum testNum = TestNum::Test1;
    if (argc > 1)
    {
        testNum = (TestNum)atoi(argv[1]);
    }
    search.testNum = testNum;
    
    const char* wu_filename = "workunit.txt";
    const char* result_filename = "result.txt";
    const char* localCheckpoint = "checkpoint.txt";
    const char* localTmpCheckpoint = "tmp_checkpoint.txt";

    search.Initialize(wu_filename, result_filename, localCheckpoint, localTmpCheckpoint);
    
    try
    {
        search.Start();
    }
    catch(const EndTest&)
    {}
    catch(const std::exception& e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
    }

    return 0;
}
