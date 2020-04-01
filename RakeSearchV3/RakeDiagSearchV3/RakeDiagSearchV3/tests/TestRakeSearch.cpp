#include "TestRakeSearch.h"

// Call order:
// Start() - generate squares
//   ProcessSquare() - copy square, calculate masks
//     PermuteRows() - generate orthogonal squares
//       ProcessOrthoSquare() - check OrthoDegree, write result
//     CheckMutualOrthogonality() - ???


void TestRakeSearch::ProcessSquare()
{
    std::cout << "{GenSquare\n";
    for (int n = 0; n < Rank; ++n)
    {
        for (int k = 0; k < Rank; ++k)
        {
            std:: cout << squareA[n][k] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "}\n";

    if (TestNum::Test1 == testNum)
    {
        if (100 == ++counter)
            throw EndTest();
    }
    else
        RakeSearch::ProcessSquare();
}

void TestRakeSearch::PermuteRows()
{
    if (TestNum::Test2 == testNum)
    {
        // TODO: print calculated masks
        if (100 == ++counter)
            throw EndTest();
    }
    else
        RakeSearch::PermuteRows();
}

void TestRakeSearch::ProcessOrthoSquare()
{
    if (TestNum::Test3 == testNum)
    {
        Square a(squareA);
        Square b(squareB);
        int orthoDegree = Square::OrthoDegree(a, b);
        
        std::cout << "{PermSquare Degree " << orthoDegree << "\n";
        for (int n = 0; n < Rank; ++n)
        {
            for (int k = 0; k < Rank; ++k)
            {
                std:: cout << squareB[n][k] << " ";
            }
            std::cout << "\n";
        }

        std::cout << "}\n";

        if (100 == ++counter)
            throw EndTest();
    }
    else
        RakeSearch::ProcessOrthoSquare();
}

//---------------------------------------------------------

void TestRakeSearch::CallBasePermuteRows()
{
    RakeSearch::PermuteRows();
}

void TestRakeSearch::CallBaseProcessSquare()
{
    RakeSearch::ProcessSquare();
}

void TestRakeSearch::CallBaseProcessOrthoSquare()
{
    RakeSearch::ProcessOrthoSquare();
}

#include "../RakeSearch.cpp"
