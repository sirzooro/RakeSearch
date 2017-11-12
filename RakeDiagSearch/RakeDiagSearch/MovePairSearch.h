// Search for pairs of diagonal Latin squares by the method of rows permutation

# if !defined MovePairSearch_h
# define MovePairSearch_h

# include <iostream>
# include <string>

# include "Generator.h"

using namespace std;

class MovePairSearch
{
public:
  static const int Rank = Square::Rank;

  MovePairSearch();              // Default constructor
  void Reset();                  // Reset search settings
  void ClearBeforeNextSearch();  // Reset the variables before the next search
  void InitializeMoveSearch(string start, string result, string checkpoint, string temp);  // Search initialization
  void StartMoveSearch();        // Start the search for orthogonal squares by the method of rows permutation
  void OnSquareGenerated(Square newSquare);  // Event processor of DLS generation, will start the search for its pair

private:                              
  static const int CheckpointInterval = 1000000;  // Interval for checkpoint creation
  static const int OrhoSquaresCacheSize = 32;     // Cache size to store the squares orthogonal to the processed one

  void MoveRows();                  // Permute the rows of the given DLS, trying to find ODLS for it
  void ProcessOrthoSquare();        // Process the found orthogonal square
  void CheckMutualOrthogonality();  // Check the mutual orthogonality of a set of squares found in the current search
  void CreateCheckpoint();          // Create a checkpoint
  void Read(std::istream& is);      // Read the search state from stream
  void Write(std::ostream& os);     // Write the search state into stream
  void ShowSearchTotals();          // Display the total results of the search
  
  static void CopyRow(int* __restrict dst, int* __restrict src);
  static void SetRow(int* dst, int val);

  Generator squareAGenerator;       // DLS generator
  int squareA[Rank][Rank];          // Initial DLS, whose rows will be permuted
  int squareB[Rank][Rank];          // Generated DLS, the rows inside which will be permuted 
  int rowsHistory[Rank];      // Array of the history of rows usage; rowsHistory[number of the row][value] = 0 | 1, where 0 means the row with the number "value" has been used for the row "number of the row" of the generated square; 1 - the row can be used.
  int currentSquareRows[Rank];      // Array listing the current rows used in the square. The number of the used row is at the i-th position
  int pairsCount;                   // The number of discovered diagonal squares in rows permutations of the found squareA
  int totalPairsCount;              // The total number of discovered diagonal squares, within the whole search
  int totalSquaresWithPairs;        // The total number of squares having at least one orthogonal found for them
  int totalProcessedSquaresLarge;   // The number of processed DLS received from the generator, billions
  int totalProcessedSquaresSmall;   // The number of processed DLS received from the generator, less than a billion

  Square orthoSquares[OrhoSquaresCacheSize];    // Cache to store the squares orthogonal to the initial one

  string startParametersFileName;   // Name of the file with start parameters for computing  
  string resultFileName;            // Name of the file with results
  string checkpointFileName;        // Name of the file with checkpoint
  string tempCheckpointFileName;    // Temporary name of the file with a new checkpoint

  int isInitialized;                // Flag of the search initialization
  int isStartFromCheckpoint;        // Flag of starting from the checkpoint

  string moveSearchGlobalHeader;    // Header preceding the data about search state
  string moveSearchComponentHeader; // Header preceding the data about the state of the component of rows permutation
  static const bool isDebug = false; // Flag of displaying debug information
};

# endif
