// Search for pairs of diagonal Latin squares by the method of rows permutation

# include "MovePairSearch.h"
# include "boinc_api.h"

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
      rowsHistory[i][j] = 1;
    }
  }

  // Reset the values in the vectors of rows usage in the next permutation
  // and the rows numbers used for the current square
  for (int i = 0; i < Rank; i++)
  {
    rowsUsage[i] = 1;
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

  // Read the generator status and the search status from the checkpoint file or initial values 
    // Open files with initial parameters and the checkpoint file
    startFile.open(startParametersFileName.c_str(), std::ios_base::in);
    checkpointFile.open(checkpointFileName.c_str(), std::ios_base::in);

    // Read the status
	if (checkpointFile.is_open())
	{
		// Read the status from the checkpoint file
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
		// Read the status from the initial parameters file
		Read(startFile);
		isStartFromCheckpoint = 0;
    }

    // Close the files
    startFile.close();
    checkpointFile.close();
}


// Read the search status from stream
void MovePairSearch::Read(istream& is)
{
  string marker;

  // Reset the initialization flag
  isInitialized = 0;

  // Read the search status
    // Find the start marker of the status
    do
	{
		std::getline(is, marker);

		if (is.eof())
		{
			throw ("Expected start marker, but EOF found.");
		}
    }
    while (marker != moveSearchGlobalHeader);
    
    // Read the status of the DLS generator
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


// Write the search status into stream
void MovePairSearch::Write(ostream& os)
{
  // Write the search status
    // Writing the header
    os << moveSearchGlobalHeader << endl;
    os << endl;

    // Write the status of the DLS generator
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

  // Fix the information about status of processing 
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

  int diagonalValues[Rank];
  int duplicationDetected = 0;

  // Write the 1st row of square A into square B for the search of normalized squares
  for (int j = 0; j < Rank; j++)
  {
    squareB[0][j] = squareA[0][j];
  }

  // Mark the usage of the 1st row, because it is fixed
  rowsUsage[0] = 0;
  rowsHistory[0][0] = 0;
  currentSquareRows[0] = 0;

  while (currentRowId > 0)
  {
    // Select a row from the initial square for the position currentRowId of the generated square
    isRowGet = 0;
    gettingRowId = -1;
    for (int i = 0; i < Rank; i++)
    {
      // Check the i-th row of the initial square
      if (rowsUsage[i] && rowsHistory[currentRowId][i])
      {
        isRowGet = 1;
        gettingRowId = i;

        break;
      }
    }

    // Process the search result
    if (isRowGet)
    {
      // Process the new found row
        // Write the new row into the square
          // Read the number of the row which is now in the square
          oldRowId = currentSquareRows[currentRowId];
          // Write the new row into the square, the flags array of the used rows 
          // into the history of the used rows, and the array of the current rows
            // Write the new row into the square
            for (int j = 0; j < Rank; j++)
            {
              squareB[currentRowId][j] = squareA[gettingRowId][j];
            }
            // Mark the row in the array of the used rows
            rowsUsage[gettingRowId] = 0;
            // Mark the row in the history of the used rows
            rowsHistory[currentRowId][gettingRowId] = 0;
            // Write the row into the array of the current rows
            currentSquareRows[currentRowId] = gettingRowId;

        // Clear usage flags of the previous row
          // Clear the mark in the array of the used rows
          if (oldRowId != -1)
          {
            rowsUsage[oldRowId] = 1;
          }

        // Check diagonality of the generated part of the square 
          // Clear the flag signalizing about duplicates on diagonals
          duplicationDetected = 0;
          // Check the main diagonal
            // Clear the flags of the used values
            for (int i = 0; i < Rank; i++)
            {
              diagonalValues[i] = 1;
            }
            // Check the values on the main diagonal
            for (int i = 0; i <= currentRowId; i++)
            {
              // Check the i-th element on the main diagonal - cell (i, i)
              if (diagonalValues[squareB[i][i]])
              {
                diagonalValues[squareB[i][i]] = 0;
              }
              else
              {
                duplicationDetected = 1;
                break;
              }
            }
          // Check the secondary diagonal if needed
          if (!duplicationDetected)
          {
            // Check the secondary diagonal
              // Clear the flags of the used values
              for (int i = 0; i < Rank; i++)
              {
                diagonalValues[i] = 1;
              }
              // Check the values on the secondary diagonal starting from its "tail"
              for (int i = 0; i <= currentRowId; i++)
              {

                // Check the i-th value on the secondary diagonal - cell (i, rank - 1 - i)
                if (diagonalValues[squareB[i][Rank - 1 - i]]) 
                {
                  diagonalValues[squareB[i][Rank - 1 - i]] = 0;
                }
                else
                {
                  duplicationDetected = 1;
                  break;
                }
              }
          }

        // Process the results of checking the square for diagonality
        if (!duplicationDetected)
        {
          // Step forward depending on the current position
          if (currentRowId == Rank - 1)
          {
            // Process the found square
            ProcessOrthoSquare();
          }
          else
          {
            // Step forward
            currentRowId++;
          }
        }
    }
    else
    {
      // Process not-founding of the new row: step backward, clear the flags of usage,
      // the history of usage, the list of current rows and clear the square itself
        // Read the number of the current row
        oldRowId = currentSquareRows[currentRowId];
        // Clear the current row
        for (int j = 0; j < Rank; j++)
        {
          squareB[currentRowId][j] = -1;
        }
        // Clear the current square
        currentSquareRows[currentRowId] = -1;
        // Clear the flag of possible usage
        rowsUsage[oldRowId] = 1;
        // Clear the history of work with this cell
        for (int i = 0; i < Rank; i++)
        {
          rowsHistory[currentRowId][i] = 1;
        }
        // Step backward
        currentRowId--;
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
