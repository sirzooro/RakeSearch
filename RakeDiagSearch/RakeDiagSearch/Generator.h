// DLS generator

# if !defined Generator_h
# define Generator_h

# include <iostream>
# include <fstream>
# include <string>

# include "Square.h"

using namespace std;

class MovePairSearch;

class Generator
{
public:
	Generator();					// Default constructor
	Generator(const Generator& source);	// Copy constructor
	void Start();					// Start the squares generation
	void Reset();					// Reset all values of internal structures
	void SetFileNames(const string& start, const string& result, const string& checkpoint, const string& temp);	// Set names for the files of parameters and checkpoints
	void Initialize(const string& start, const string& result, const string& checkpoint, const string& temp);	// Initialize the generator
	void Subscribe(MovePairSearch *search);		// Subscribe to the event of square generation
	void Unsubscribe();				            // Unsubscribe from the event of square generation

 	Generator& operator = (const Generator&  value);									// Assignment operator
	friend std::ostream& operator << (std::ostream& os, Generator& value);		// Operator of writing the generator state
	friend std::istream& operator >> (std::istream& is, Generator& value);		// Operator of reading the generator state

//protected:
private:
	static const int Rank = Square::Rank;	// Rank (for convenience)
	static const int Free = 1;				// Flag of the value free for usage
	static const int Used = 0;				// Flag of the value being used in a diagonal | row | column 
	static const int Yes = 1;				// Flag "Yes"
	static const int No = 0;				// Flag "No"
	static const int MaxCellsInPath = Rank*Rank;	// Maximal number of processed cells
	int cellsInPath;				       // Number of processed cells

	void CopyState(const Generator& source);		// Copy the state
	void Read(std::istream& is);			// Read the generator state from stream
	void Write(std::ostream& os);			// Write the generator state into stream

	Square newSquare;				        // Generated square

	int path[MaxCellsInPath][2];			// Path to fill the square matrix: path[i][0] is the row on step i, path[i][1] is the column
	int keyRowId;					        // Row ID of the key cell, the value of which stops computation  
	int keyColumnId;						// Column ID of the key cell
	int keyValue;							// Value of the key cell which will stop computation

	int primary;							// Contents of the main diagonal
	int secondary;							// Contents of the secondary diagonal
	int columns[Rank];						// Matrix of the values used in columns: columns[column][value] = 0 | 1; 0 means the value has been used, 1 - it is free.
											// Note: this matrix is transposed copy of columns input data.
	int rows[Rank];							// Matrix of the values used in rows: rows[row][value] = 0 | 1
	int cellsHistory[Rank][Rank];			// Cube of the values used for generating the part of the square - cellsHistory[row][column][value]

	string startParametersFileName;			// Name of the file with start parameters for computing 
	string resultFileName;					// Name of the file with results
	string checkpointFileName;				// Name of the file with checkpoint
	string tempCheckpointFileName;			// Temporary name of the file with a new checkpoint

	int isInitialized;						// Flag of the search initialization
	int squaresCount;						// The number of discovered DLS
	int rowId;								// Row ID of the processed cell
	int columnId;							// Column ID of the processed cell
	int cellId;								// Cell ID in the path of the square bypassing

private:
	void CreateCheckpoint();				// Create a checkpoint
	void ProcessSquare();					// Process the found square

	string generatorStateHeader;			// Header preceding the data about search state in start parameters file or a checkpoint file
	MovePairSearch* subscriber;			    // Object of DLS search by rows permutation
};

# endif
