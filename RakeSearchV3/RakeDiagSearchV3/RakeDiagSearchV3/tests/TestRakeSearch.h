#ifndef TEST_RAKE_SEARCH_H
#define TEST_RAKE_SEARCH_H

#include "../RakeSearch.h"
#include <exception>

class EndTest : public std::exception
{};

enum class TestNum
{
    Test1 = 1,
    Test2 = 2,
    Test3 = 3,
};

class TestRakeSearch : public RakeSearch
{
public:
    void PermuteRows() override;
    void ProcessSquare() override;
    void ProcessOrthoSquare() override;
    
    void CallBasePermuteRows();
    void CallBaseProcessSquare();
    void CallBaseProcessOrthoSquare();
    
    int counter = 0;
    TestNum testNum = TestNum::Test1;
};

#endif // TEST_RAKE_SEARCH_H
