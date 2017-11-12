// Search for pairs of diagonal Latin squares by the method of rows permutation

#include "MovePairSearch.h"
#include "boinc_api.h"

#ifdef __SSE2__
#include "immintrin.h"
#endif

#ifndef _MSC_VER
#define ffs __builtin_ffs
#else
#include <intrin.h>
inline int ffs(int i)
{
  unsigned long index;
  if (_BitScanForward(&index, i))
    return index + 1;
  else
    return 0;
}
#endif

#define SetBit(bitfield, bitno) ((bitfield) |= (1u << (bitno)))
#define ClearBit(bitfield, bitno) ((bitfield) &= ~(1u << (bitno)))

#define AllBitsMask(numbits) ((1u << (numbits)) - 1)

inline void MovePairSearch::CopyRow(int* __restrict dst, int* __restrict src)
{
  int n = 0;
#ifdef __AVX__
  for (; n < Rank - 7; n += 8)
  {
    __m256i v = _mm256_loadu_si256((const __m256i*)&src[n]);
    _mm256_storeu_si256((__m256i*)&dst[n], v);
  }
#endif
#ifdef __SSE2__
  for (; n < Rank - 3; n += 4)
  {
    __m128i v = _mm_loadu_si128((const __m128i*)&src[n]);
    _mm_storeu_si128((__m128i*)&dst[n], v);
  }
#endif
  for (; n < Rank; n++)
  {
    dst[n] = src[n];
  }
}

inline void MovePairSearch::SetRow(int* dst, int val)
{
  int n = 0;
#ifdef __AVX__
  __m256i v_avx = _mm256_set1_epi32(val);
  for (; n < Rank - 7; n += 8)
  {
    _mm256_storeu_si256((__m256i*)&dst[n], v_avx);
  }
#endif
#ifdef __SSE2__
#ifdef __AVX__
  __m128i v_sse = _mm256_castsi256_si128(v_avx);
#else
  __m128i v_sse = _mm_set1_epi32(val);
#endif
  for (; n < Rank - 3; n += 4)
  {
    _mm_storeu_si128((__m128i*)&dst[n], v_sse);
  }
#endif
  for (; n < Rank; n++)
  {
    dst[n] = val;
  }
}

// Default constructor
MovePairSearch::MovePairSearch()
{
  Reset();
}

// Reset search settings
void MovePairSearch::Reset()
{
  // Reset the DLS generator settings
  squareAGenerator.Reset();

  // Reset the next search settings
  ClearBeforeNextSearch();

  // Reset the values of global counters
  totalPairsCount = 0;
  totalSquaresWithPairs = 0;
  totalProcessedSquaresSmall = 0;
  totalProcessedSquaresLarge = 0;

  startParametersFileName = "start_parameters.txt";
  resultFileName = "result.txt";
  checkpointFileName = "checkpoint.txt";
  tempCheckpointFileName = "tmp_checkpoint.txt";

  // Set the header constant for the file of parameters or checkpoint
  moveSearchGlobalHeader = "# Move search of pairs OLDS status";
  moveSearchComponentHeader = "# Move search component status";

  // Reset the initialization flag
  isInitialized = 0;
}


// Reset the variables before the next search
void MovePairSearch::ClearBeforeNextSearch()
{
  // Reset the values of matrices of squares A and B, 
  // also the matrix of rows usage when forming square B
  for (int i = 0; i < Rank; i++)
  {
    for (int j = 0; j < Rank; j++)
    {
      squareA[i][j] = -1;
      squareB[i][j] = -1;
    }
    rowsHistory[i] = AllBitsMask(Rank);
  }

  // Reset the values in the vectors of rows usage in the next permutation
  // and the rows numbers used for the current square
  for (int i = 0; i < Rank; i++)
  {
    currentSquareRows[i] = -1;
  }

  // Reset the counter of pairs found for the given DLS
  pairsCount = 0;
}


// Search initialization
void MovePairSearch::InitializeMoveSearch(string start, string result, 
                                       string checkpoint, string temp)
{
  fstream startFile;
  fstream checkpointFile;

  // Read the filenames
  startParametersFileName = start;
  resultFileName = result;
  checkpointFileName = checkpoint;
  tempCheckpointFileName = temp;

  // Read the generator state and the search state from the checkpoint file or initial values 
    // Open files with initial parameters and the checkpoint file
    startFile.open(startParametersFileName.c_str(), std::ios_base::in);
    checkpointFile.open(checkpointFileName.c_str(), std::ios_base::in);

    // Read the state
	if (checkpointFile.is_open())
	{
		// Read the state from the checkpoint file
		try
		{
			Read(checkpointFile);
			isStartFromCheckpoint = 1;
		}
		catch (...)
		{
			cerr << "Error opening checkpoint file! Starting with workunit start parameters." << endl;
			isStartFromCheckpoint = 0;
		}
    }

	if (isStartFromCheckpoint != 1)
    {
		// Read the state from the initial parameters file
		Read(startFile);
		isStartFromCheckpoint = 0;
    }

    // Close the files
    startFile.close();
    checkpointFile.close();
}


// Read the search state from stream
void MovePairSearch::Read(istream& is)
{
  string marker;

  // Reset the initialization flag
  isInitialized = 0;

  // Read the search state
    // Find the start marker of the state
    do
	{
		std::getline(is, marker);

		if (is.eof())
		{
			throw ("Expected start marker, but EOF found.");
		}
    }
    while (marker != moveSearchGlobalHeader);
    
    // Read the state of the DLS generator
    is >> squareAGenerator;

    // Find the marker of permutation component
	do
	{
		std::getline(is, marker);

		if (is.eof())
		{
			throw ("Expected start marker, but EOF found.");
		}
	}
    while (marker != moveSearchComponentHeader);

    // Read the search variables by permutation (actually, the variables containing statistics)
    is >> pairsCount;
    is >> totalPairsCount;
    is >> totalSquaresWithPairs;
    is >> totalProcessedSquaresLarge;
    is >> totalProcessedSquaresSmall;

  // Set the initialization flag
  isInitialized = 1;
}


// Write the search state into stream
void MovePairSearch::Write(ostream& os)
{
  // Write the search state
    // Writing the header
    os << moveSearchGlobalHeader << endl;
    os << endl;

    // Write the state of the DLS generator
    os << squareAGenerator;
    os << endl;

    // Write the header of permutation block
    os << moveSearchComponentHeader << endl;
    os << endl;

    // Writing statistical values
    os << pairsCount << " " << totalPairsCount << " " << totalSquaresWithPairs << endl;
    os << totalProcessedSquaresLarge << " " << totalProcessedSquaresSmall << endl;
    os << endl;
}


// Create a checkpoint
void MovePairSearch::CreateCheckpoint()
{
  ofstream checkpointFile;

  checkpointFile.open(tempCheckpointFileName.c_str(), std::ios_base::out);
  if (checkpointFile.is_open())
  {
    Write(checkpointFile);
    checkpointFile.close();
    remove(checkpointFileName.c_str());
    rename(tempCheckpointFileName.c_str(), checkpointFileName.c_str());
  }
  else
    cerr << "Error opening checkpoint file!" << endl;
}


// Start the search for orthogonal squares by the method of rows permutation
void MovePairSearch::StartMoveSearch()
{
  // Subscribe to the event of finding a DLS
  squareAGenerator.Subscribe(this); 

  // Start DLS generation
  squareAGenerator.Start();

  // Unsubscribe from the event of finding a DLS
  squareAGenerator.Unsubscribe();

  // Output the search results
  ShowSearchTotals();
}


// Event processor of DLS generation, will start the search for its pair 
void MovePairSearch::OnSquareGenerated(Square newSquare)
{
  // Reset before the search for orthogonal squares
  ClearBeforeNextSearch();

  // Write the found square
  for (int i = 0; i < Rank; i++)
  {
    for (int j = 0; j < Rank; j++)
    {
      squareA[i][j] = newSquare.Matrix[i][j];
    }
  }

  // Start the rows permutation
  MoveRows();

  // Check the mutual orthogonality of the squares
  if (pairsCount > 0)
  {
    CheckMutualOrthogonality();
  }

  // Gather statistics on the processed squares
  totalProcessedSquaresSmall++;

  // Fix the information about state of processing 
  if (totalProcessedSquaresSmall % CheckpointInterval == 0)
  {
    // Update the completion progress for the BOINC client
    double fraction_done;
    if(Rank == 8)
      fraction_done = (double)(totalProcessedSquaresSmall)/5000000.0;
    else 
    {
      if(Rank == 9)
        fraction_done = (double)(totalProcessedSquaresSmall)/275000000.0;
      else
        fraction_done = 0.999999999;
    }

    if(fraction_done >=1) fraction_done = (double)(totalProcessedSquaresSmall)/300000000.0;
    if(fraction_done >=1) fraction_done = 0.999999999;

    boinc_fraction_done(fraction_done); // Tell the BOINC client the completion fraction

    // Check if the BOINC client is able to create a checkpoint,
    // and if yes, start the function of its creation
    if (boinc_time_to_checkpoint()) {
      CreateCheckpoint();
      boinc_checkpoint_completed(); // BOINC client knows the checkpoint has been created
    }  

    if(isDebug)
    {
      cout << "# ------------------------" << endl;
      cout << "# Processed " << totalProcessedSquaresLarge << " milliards and " << totalProcessedSquaresSmall << " squares." << endl;
      cout << "# Last processed square:" << endl;
      cout << endl;
      cout << newSquare;
      cout << "# ------------------------" << endl;
    }
  }
}

// Permute the rows of the given DLS, trying to find ODLS for it
void MovePairSearch::MoveRows()
{
  int currentRowId = 1;
  int isRowGet = 0;
  int gettingRowId = -1;
  int oldRowId = -1;

  int diagonalValues, diagonalValues2;
  int duplicationDetected = 0;
  
  int diagonalValuesHistory[Rank][2];
  
  int rowsUsage; // Flags of the rows usage at the current moment; rowsUsage[number of the row] = 0 | 1, where 0 means the row is already used, 1 - not.

  // Write the 1st row of square A into square B for the search of normalized squares
  CopyRow(&squareB[0][0], &squareA[0][0]);

  // Mark the usage of the 1st row, because it is fixed
  rowsUsage = AllBitsMask(Rank) & ~1u;
  ClearBit(rowsHistory[0], 0);
  currentSquareRows[0] = 0;
  
  // Set bits for diagonal values in 1st row
  diagonalValuesHistory[0][0] = 1u << squareB[0][0];
  diagonalValuesHistory[0][1] = 1u << squareB[0][Rank - 1];
  
  diagonalValues = diagonalValuesHistory[0][0];
  diagonalValues2 = diagonalValuesHistory[0][1];

  while (1)
  {
    // Select a row from the initial square for the position currentRowId of the generated square
    // Process the search result
    gettingRowId = rowsUsage & rowsHistory[currentRowId];
    if (gettingRowId)
    {
      gettingRowId = ffs(gettingRowId) - 1;
      // Process the new found row

      // Mark the row in the history of the used rows
      ClearBit(rowsHistory[currentRowId], gettingRowId);

      // Check diagonality of the generated part of the square
      // Check the main diagonal and secondary diagonal
      // Calculate bits for current row
      int bit1 = 1u << squareA[gettingRowId][currentRowId];
      int bit2 = 1u << squareA[gettingRowId][Rank - 1 - currentRowId];

      // Duplicate check
      duplicationDetected = ((0 != (diagonalValues & bit1)) || (0 != (diagonalValues2 & bit2)));

        // Process the results of checking the square for diagonality
        if (!duplicationDetected)
        {
          // Write the new row into the square
          CopyRow(&squareB[currentRowId][0], &squareA[gettingRowId][0]);
          
          // Write the row into the array of the current rows
          currentSquareRows[currentRowId] = gettingRowId;

          // Step forward depending on the current position
          if (currentRowId == Rank - 1)
          {
            // Process the found square
            ProcessOrthoSquare();
          }
          else
          {
            // Save new bitmasks for diagonal values for further use
            diagonalValues |= bit1;
            diagonalValues2 |= bit2;
            diagonalValuesHistory[currentRowId][0] = diagonalValues;
            diagonalValuesHistory[currentRowId][1] = diagonalValues2;

            // Mark the row in the array of the used rows
            ClearBit(rowsUsage, gettingRowId);

            // Step forward
            currentRowId++;
          }
        }
    }
    else
    {
      // Process not-founding of the new row: step backward, clear the flags of usage,
      // the history of usage, the list of current rows and clear the square itself

        // Clear the current row
        SetRow(&squareB[currentRowId][0], -1);
        // Clear the current square
        currentSquareRows[currentRowId] = -1;
        // Clear the history of work with this cell
        rowsHistory[currentRowId] = AllBitsMask(Rank);
        // Step backward
        currentRowId--;
        // Check if we are done
        if (0 == currentRowId)
          break;
        // Get saved values for previous row
        diagonalValues = diagonalValuesHistory[currentRowId-1][0];
        diagonalValues2 = diagonalValuesHistory[currentRowId-1][1];
        // Clear the flag of row usage
        SetBit(rowsUsage, currentSquareRows[currentRowId]);
    }
  }
}


// Process the found orthogonal square
void MovePairSearch::ProcessOrthoSquare()
{
  int isDifferent = 0;      // The number of differences in the rows with the initial square (to avoid generating its copy)
  ofstream resultFile;      // The stream for output into the results file

  Square a(squareA);        // Square A as an object
  Square b(squareB);        // Square B as an object

  int orthoMetric = Rank*Rank;  // The value of orthogonality metric saying that the squares are fully orthogonal

  // Check the square for being a copy of the initial one
  isDifferent = 0;
  
  for (int i = 0; i < Rank; i++)
  {
    if (currentSquareRows[i] != i)
    {
      isDifferent = 1;
      break;
    }
  }

  // Process the found square
  if (isDifferent && Square::OrthoDegree(a, b) == orthoMetric 
      && b.IsDiagonal() && b.IsLatin() && a.IsDiagonal() && b.IsLatin())
  {
    // Write the information about the found square
      // Increase the squares counter
      pairsCount++;
      totalPairsCount++;

      // Remember the basic square
      if (pairsCount == 1)
      {
        orthoSquares[pairsCount - 1] = a;
        totalSquaresWithPairs++;
      }

      // Remember the pair square
      if (pairsCount < OrhoSquaresCacheSize)
      {
        orthoSquares[pairsCount] = b;
      }

      // Output the header
      if (pairsCount == 1)
      {
        if(isDebug)
        {
          // Output the information about the 1st square in a pair as a header 
          cout << "{" << endl;
          cout << "# ------------------------" << endl;
          cout << "# Detected pair for the square: " << endl;
          cout << "# ------------------------" << endl;
          cout << a;
          cout << "# ------------------------" << endl;
        }
        // Write the information into file
        resultFile.open(resultFileName.c_str(), std::ios_base::binary | std::ios_base::app);
        if (resultFile.is_open())
        {
          resultFile << "{" << endl;
          resultFile << "# ------------------------" << endl;
          resultFile << "# Detected pair for the square: " << endl;
          resultFile << "# ------------------------" << endl;
          resultFile << a;
          resultFile << "# ------------------------" << endl;
          resultFile.close();
        }
        else
        {
          std::cerr << "Error opening file!";
        }
      }

      // Output the information about the found pair
        if(isDebug)
        {
          // Output the information into the console
          cout << b << endl;
        }

        // Output the information into the file
        resultFile.open(resultFileName.c_str(), std::ios_base::binary | std::ios_base::app);
        if (resultFile.is_open())
        { 
          resultFile << b << endl;
          resultFile.close();
        }
        else
        {
          std::cerr << "Error opening file!";
        }

  }
}


// Check the mutual orthogonality of a set of squares found in the current search
void MovePairSearch::CheckMutualOrthogonality()
{
  int orthoMetric = Rank*Rank;
  int maxSquareId;
  ofstream resultFile;

  // Detect the upper bound of the squares to process
  if (pairsCount < OrhoSquaresCacheSize)
  {
    maxSquareId = pairsCount;
  }
  else
  {
    maxSquareId = OrhoSquaresCacheSize - 1;
  }

  // Open the file with results
  resultFile.open(resultFileName.c_str(), std::ios_base::binary | std::ios_base::app);
  if (!resultFile.is_open()) { cout << "Error opening file!"; return; }

  // Check the mutual orthogonality of a set of squares 
  for (int i = 0; i <= maxSquareId; i++)
  {
    for (int j = i + 1; j <= maxSquareId; j++)
    {
      if (Square::OrthoDegree(orthoSquares[i], orthoSquares[j]) == orthoMetric)
      {
        if(isDebug) cout << "# Square " << i << " # " << j << endl;
        resultFile << "# Square " << i << " # " << j << endl;
      }
    }
  }
  if(isDebug) cout << endl;
  resultFile << endl;

  // Write the total number of the found DLS pairs
  if(isDebug) cout << "# Pairs found: " << pairsCount << endl;
  resultFile << "# Pairs found: " << pairsCount << endl;

  // Set the end marker of the results section
  if(isDebug) cout << "}" << endl;
  resultFile << "}" << endl;

  // Close the file with results
  resultFile.close();
}

// Display the total results of the search
void MovePairSearch::ShowSearchTotals()
{
  ofstream resultFile;

  if(isDebug)
  {
    // Write the results into the console
    cout << "# ------------------------" << endl;
    cout << "# Total pairs found: " << totalPairsCount << endl;
    cout << "# Total squares with pairs: " << totalSquaresWithPairs << endl;
    cout << "# ------------------------" << endl;
  }

  // Write the results into the file
  resultFile.open(resultFileName.c_str(), std::ios_base::binary | std::ios_base::app);
  if (resultFile.is_open())
  {
    resultFile << "# ------------------------" << endl;
    resultFile << "# Total pairs found: " << totalPairsCount << endl;
    resultFile << "# Total squares with pairs: " << totalSquaresWithPairs << endl;
    resultFile << "# Processes " << totalProcessedSquaresLarge << " milliards " << totalProcessedSquaresSmall << " squares" << endl;
    resultFile << "# ------------------------" << endl;
    resultFile.close();
  }
  else cerr << "Error opening file!" << endl;
}
