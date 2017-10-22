#include "Square.h"

using namespace std;

// Default constructor. "Zeroing" all values
Square::Square()
{
	Reset();
}


// Constructing a square by the given matrix
Square::Square(int source[Rank][Rank])
{
	Initialize(source);
}


// Copy constructor
Square::Square(Square& source)
{
	Initialize(source.Matrix);
}


// Initialization of internal structures
void Square::Initialize(int source[Rank][Rank])
{
	for (int rowId = 0; rowId < Rank; rowId++)
	{
		for (int columnId = 0; columnId < Rank; columnId++)
		{
			Matrix[rowId][columnId] = source[rowId][columnId];
		}
	}
}


// Resetting the values of internal variables
void Square::Reset()
{
	for (int rowId = 0; rowId < Rank; rowId++)
	{
		for (int columnId = 0; columnId < Rank; columnId++)
		{
			Matrix[rowId][columnId] = Empty;
		}
	}
}


// Comparison operator
int Square::operator == (Square& value)
{
	int isEqual = 1;

	for (int rowId = 0; rowId < Rank && isEqual; rowId++)
	{
		for (int columnId = 0; columnId < Rank && isEqual; columnId++)
		{
			if (Matrix[rowId][columnId] != value.Matrix[rowId][columnId])
			{
				isEqual = 0;
			}
		}
	}

return isEqual;
}


// Assignment operator
Square& Square::operator = (Square& value)
{
	Initialize(value.Matrix);

return *this;
}


// Operator of square data output
std::ostream& operator << (std::ostream& os, Square& value)
{
	value.Write(os);

return os;
}


// Operator of square data input
std::istream& operator >> (std::istream& is, Square& value)
{
	value.Read(is);

return is;
}


// Reading the square from the stream
void Square::Read(std::istream& is)
{
	char readedChar;	// "Buffer" used for char-by-char reading

	// Reading the start token of data block
	do
	{
		is >> readedChar;
	}
	while (readedChar != HeadToken);

	// Reading the matrix of the square
	for (int rowId = 0; rowId < Rank; rowId++)
	{
		for (int columnId = 0; columnId < Rank; columnId++)
		{
			is >> Matrix[rowId][columnId];
		}
	}

	// Reading the end token of data block
	do
	{
		is >> readedChar;
	}
	while (readedChar != TailToken);
}


// Writing the square into the stream
void Square::Write(std::ostream& os)
{
	// Writing the start token of data block
	os << HeadToken << endl;

	// Writing the matrix of the square
	for (int rowId = 0; rowId < Rank; rowId++)
	{
		for (int columnId = 0; columnId < Rank; columnId++)
		{
			os << Matrix[rowId][columnId] << " ";
		}
		os << endl;
	}

	// Writing the start token of data block
	os << TailToken << endl;
}

// Checking the square for being diagonal Latin
int Square::IsDiagonal()
{
	int isDiagonal = 1;

	// Checking the main diagonal: [0;0] - [rank;rank]
	for (int itemId = 0; itemId < Rank && isDiagonal; itemId++)
	{
        // Checking the element [itemId; itemId] to match any other element of the main diagonal
		for (int comparedId = itemId + 1; comparedId < Rank && isDiagonal; comparedId++)
		{
			if (Matrix[itemId][itemId] == Matrix[comparedId][comparedId])
			{
				isDiagonal = 0;
			}
		}
	}

	// Checking the secondary diagonal: [rank - itemId - 1; itemId]
	for (int itemId = 0; itemId < Rank && isDiagonal; itemId ++)
	{
		// Checking the element [rank - itemId - 1; itemId] to match any other element of the secondary diagonal 
		for (int comparedId = itemId + 1; comparedId < Rank && isDiagonal; comparedId++)
		{
			if (Matrix[(Rank - itemId - 1)][itemId] == Matrix[(Rank - comparedId - 1)][comparedId])
			{
				isDiagonal = 0;
			}
		}
	}

	return isDiagonal;
}


// Checking the square for being Latin
int Square::IsLatin()
{
	int isLatin = 1;

	// Checking the columns
	for (int columnId = 0; columnId < Rank && isLatin; columnId++)
	{
		// Checking correctness of the column columnId
		for (int rowId = 0; rowId < Rank && isLatin; rowId++)
		{
			// Checking the element [rowId; columnId] to match any other element of the column columndId
			for (int comparedRowId = rowId + 1; comparedRowId < Rank && isLatin; comparedRowId++)
			{
				if (Matrix[comparedRowId][columnId] == Matrix[rowId][columnId])
				{
					isLatin = 0;
				}
			}
		}
	}

	// Checking the rows
	for (int rowId = 0; rowId < Rank && isLatin; rowId++)
	{
		// Checking correctness of the row rowId
		for (int columnId = 0; columnId < Rank && isLatin; columnId++)
		{
			// Checking the element [rowId; columnId] to match any other element of the row rowId
			for (int comparedColumnId = columnId + 1; comparedColumnId < Rank && isLatin; comparedColumnId++)
			{
				if (Matrix[rowId][columnId] == Matrix[rowId][comparedColumnId])
				{
					isLatin = 0;
				}
			}
		}
	}

	return isLatin;
}


// Checking the orthogonality of the squares a and b
int Square::OrthoDegree(Square a, Square b)
{
	int degree = 0;				// Orthogonality degree
	int freePair[Rank][Rank];	// Array of usage of the values pairs in the generated Greco-Latin square

	// Initialize all pairs as non-used, "free"
	for (int i = 0; i < Rank; i++)
	{
		for (int j = 0; j < Rank; j++)
		{
			freePair[i][j] = 1;
		}
	}

	// Mark the pairs used in the Greco-Latin square
	for (int rowId = 0; rowId < Rank; rowId++)
	{
		for (int columnId = 0; columnId < Rank; columnId++)
		{
			if (freePair[a.Matrix[rowId][columnId]][b.Matrix[rowId][columnId]])
			{
				freePair[a.Matrix[rowId][columnId]][b.Matrix[rowId][columnId]] = 0;
				degree++;
			}
		}
	}

return degree;
}
