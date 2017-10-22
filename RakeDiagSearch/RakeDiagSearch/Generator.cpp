// DLS generator

# include "Generator.h"
# include "MovePairSearch.h"

using namespace std;

// Default constructor
Generator::Generator()
{
  // Reset the settings
  Reset();

  // Set the text constants
  generatorStateHeader = "# Generation of DLS status";
}


// Copy constructor
Generator::Generator(Generator& source)
{
  CopyState(source);
}


// Reset all values of internal structures
void Generator::Reset()
{
  // Reset internal values of the square
  newSquare.Reset();

  // Reset the values of the square generating structures
    // Reset the values corresponding to the key cell
    keyRowId = Square::Empty;
    keyColumnId = Square::Empty;
    keyValue = Square::Empty;

    // Reset the values connected with the path of filling the cells 
    for (int i = 0; i < MaxCellsInPath; i++)
    {
      path[i][0] = Square::Empty;
      path[i][1] = Square::Empty;
    }

    // Reset the values in the vectors of diagonal elements usage
    for (int i = 0; i < Rank; i++)
    {
      primary[i] = Free;
      secondary[i] = Free;
    }

    // Reset the values in the matrices of rows/columns elements usage
    for (int i = 0; i < Rank; i++)
    {
      for (int j = 0; j < Rank; j++)
      {
        columns[i][j] = Free;
        rows[i][j] = Free;
      }
    }

    // Reset the values in the cube of the history of cell values usage
    for (int i = 0; i < Rank; i++)
    {
      for (int j = 0; j < Rank; j++)
      {
        for (int h = 0; h < Rank; h++)
        {
          cellsHistory[i][j][h] = Free;
        }
      }
    }

    // Reset coordinates of the processed cell
    rowId = Square::Empty;
    columnId = Square::Empty;

    // Reset filenames
    checkpointFileName.clear();
    tempCheckpointFileName.clear();
    resultFileName.clear();

    // Reset the number of generated squares
    squaresCount = 0;

    // Reset the initialization flag
    isInitialized = No;

    // Reset the pointer to a subscriber
    subscriber = 0;
}


// Initialize the generator
void Generator::Initialize(string start, string result, string checkpoint, string temp)
{
  fstream startFile;
  fstream checkpointFile;

  // Reset the values of internal structures 
  Reset();

  // Remember filenames of the configuration, checkpoints and results
  startParametersFileName = start;
  checkpointFileName = checkpoint;
  resultFileName = result;
  tempCheckpointFileName = temp;

  // Read the settings
  startFile.open(startParametersFileName.c_str(), std::ios_base::in);
  checkpointFile.open(checkpointFileName.c_str(), std::ios_base::in);

  if (checkpointFile.is_open())
  {
    // Read the data from checkpoint file
    Read(checkpointFile);
  }
  else
  {
    // Read the data from start parameters file
    if (startFile.is_open())
    {
      Read(startFile);
    }
  }

  startFile.close();
  checkpointFile.close();
}

// Operator of writing the generator state
std::ostream& operator << (std::ostream& os, Generator& value)
{
  value.Write(os);

return os;
}


// Operator of reading the generator state
std::istream& operator >> (std::istream& is, Generator& value)
{
  value.Read(is);

return is;
}


// Read the generator state from stream
void Generator::Read(std::istream& is)
{
  int rankToVerify;
  int result = Yes;
  string marker;

  // Reset initialization flag
  isInitialized = No;

  // Search for header
  do
  {
    std::getline(is, marker);
  }
  while (marker != generatorStateHeader);

  // Read the rank from stream
  is >> rankToVerify;

  // Read the data for search of desired rank
  if (rankToVerify == Square::Rank)
  {
    // Read the square from stream
    is >> newSquare;

    // Read the number of cells in the path 
    is >> cellsInPath;

    // Read the path of cells bypassing
    for (int i = 0; i < cellsInPath; i++)
    {
      is >> path[i][0];
      is >> path[i][1];
    }

    // Read the information about the key cell
    is >> keyRowId;
    is >> keyColumnId;
    is >> keyValue;

    // Read the information about the processed cell
    is >> rowId;
    is >> columnId;
    is >> cellId;

    // Read the information about used values and the history of values
      // Read the information about the main diagonal values
      for (int i = 0; i < Rank; i++)
      {
        is >> primary[i];
      }

      // Read the information about the secondary diagonal values
      for (int i = 0; i < Rank; i++)
      {
        is >> secondary[i];
      }

      // Read the information about values in rows
      for (int i = 0; i < Rank; i++)
      {
        for (int j = 0; j < Rank; j++)
        {
          is >> rows[i][j];
        }
      }

      // Read the information about values in columns
      for (int i = 0; i < Rank; i++)
      {
        for (int j = 0; j < Rank; j++)
        {
          is >> columns[i][j];
        }
      }

      // Read the information about the history of values in square cells
      for (int h = 0; h < Rank; h++)
      {
        for (int i = 0; i < Rank; i++)
        {
          for (int j = 0; j < Rank; j++)
          {
            is >> cellsHistory[i][j][h];
          }
        }  
      }  

    // Read the number of generated squares
    is >> squaresCount;

    // Set initialization flag
    isInitialized = Yes;
  }
}


// Write the generator state into stream
void Generator::Write(std::ostream& os)
{
  // Write the header
  os << generatorStateHeader << endl << endl;

  // Write the rank
  os << Square::Rank << endl;

  // Write the square
  os << newSquare;

  // Write the number of cells in path
  os << cellsInPath << endl;
  os << endl;
  
  // Write the path of cells bypassing
  for (int i = 0; i < cellsInPath; i++)
  {
    os << path[i][0] << " ";
    os << path[i][1] << " ";
    os << endl;
  }
  os << endl;

  // Write the information about the key cell
  os << keyRowId << " " << keyColumnId << " " << keyValue << endl;

  // Write the information about the current cell
  os << rowId << " " << columnId << " " << cellId  << endl;

  // Write an empty line for convenience
  os << endl;

  // Write the information about used values and the history of values
    // Write the information about main diagonal values
    for (int i = 0; i < Rank; i++)
    {
      os << primary[i] << " ";
    }
    os << endl;

    // Write the information about secondary diagonal values
    for (int i = 0; i < Rank; i++)
    {
      os << secondary[i] << " ";
    }
    os << endl;

    // Additional empty line
    os << endl;

    // Write the information about rows values
    for (int i = 0; i < Rank; i++)
    {
      for (int j = 0; j < Rank; j++)
      {
        os << rows[i][j] << " ";
      }
      os << endl;
    }
    os << endl;

    // Write the information about columns values
    for (int i = 0; i < Rank; i++)
    {
      for (int j = 0; j < Rank; j++)
      {
        os << columns[i][j] << " ";
      }
      os << endl;
    }
    os << endl;

    // Write the information about the history of values in cells
    for (int h = 0; h < Rank; h++)
    {
      for (int i = 0; i < Rank; i++)
      {
        for (int j = 0; j < Rank; j++)
        {
          os << cellsHistory[i][j][h] << " ";
        }
        os << endl;
      }
      os << endl;
    }
    os << endl;

  // Write the number of generated squares
  os << squaresCount << endl;
}


// Assignment operator
Generator& Generator::operator = (Generator& value)
{
  CopyState(value);

return *this;
}


// Copy the state from the given object
void Generator::CopyState(Generator& source)
{
  // Copy variables connected with the path of cells bypassing
  for (int i = 0; i < cellsInPath; i++)
  {
    path[i][0] = source.path[i][0];
    path[i][1] = source.path[i][1];
  }

  keyRowId = source.keyRowId;
  keyColumnId = source.keyColumnId;
  keyValue = source.keyValue;

  // Copy the flag arrays of used values
  for (int i = 0; i < Rank; i++)
  {
    primary[i] = source.primary[i];
    secondary[i] = source.secondary[i];
  }

  for (int i = 0; i < Rank; i++)
  {
    for (int j = 0; j < Rank; j++)
    {
      columns[i][j] = source.columns[i][j];
      rows[i][j] = source.rows[i][j];
    }
  }

  for (int i = 0; i < Rank; i++)
  {
    for (int j = 0; j < Rank; j++)
    {
      for (int h = 0; h < Rank; h++)
      {
        cellsHistory[i][j][h] = source.cellsHistory[i][j][h];
      }
    }
  }

  // Copy the filenames
  startParametersFileName = source.startParametersFileName;
  resultFileName = source.resultFileName;
  checkpointFileName = source.checkpointFileName;
  tempCheckpointFileName = source.tempCheckpointFileName;

  // Copy the variables of the current state
  isInitialized = source.isInitialized;
  squaresCount = source.squaresCount;
  rowId = source.rowId;
  columnId = source.columnId;
  cellId = source.cellId;

  // Copy the addresses of text constants
  generatorStateHeader = source.generatorStateHeader;

  // Copy the link to the subscriber
  subscriber = source.subscriber;
}


// Set names for the files of parameters and checkpoints
void Generator::SetFileNames(string start, string result, string checkpoint, string temp)
{
  startParametersFileName = start;
  resultFileName = result;
  checkpointFileName = checkpoint;
  tempCheckpointFileName = temp;
}


// Create a checkpoint
void Generator::CreateCheckpoint()
{
  fstream newCheckpointFile;

  // Write settings into a new file of checkpoint
  newCheckpointFile.open(tempCheckpointFileName.c_str(), std::ios_base::out);

  if (newCheckpointFile.is_open())
  {
    Write(newCheckpointFile);
    newCheckpointFile.close();
    remove(checkpointFileName.c_str());
    rename(tempCheckpointFileName.c_str(), checkpointFileName.c_str());
  }
}


// Start the squares generation
void Generator::Start()
{
  int isGet;        // Flag of getting new value for the cell
  int cellValue;    // New value for the cell
  int oldCellValue; // Old value from the cell

  int stop = 0;     // Flag of finishing the computing

  if (isInitialized == Yes)
  {
    // Selection of the cells values
    do
    {
      // Selection of the value for the next cell
        // Read coordinates of the cell
        rowId = path[cellId][0];
        columnId = path[cellId][1];

        // Generate new value for the cell (rowId, columnId)
          // Reset variables values
          isGet = 0;
          cellValue = Square::Empty;

          // Select the value for the cell
          for (int i = 0; i < Rank && !isGet; i++)
          {
            // Check the i value for possibility to be written into the cell (rowId, columnId)
            if (columns[i][columnId] && rows[rowId][i] && cellsHistory[rowId][columnId][i])
            {
              // The value is not used in rows and columns, the diagonals are to check
                // Set the flag which may be unset by diagonal checking
                isGet = 1;
                // Test the value: has it been used in diagonals
                  // Test the main diagonal
                  if(columnId == rowId)
                  {
                    if (!primary[i])
                    {
                      isGet = 0;
                    }
                  }

                  // Test the secondary diagonal
                  if (rowId == Rank - 1 - columnId)
                  {
                    if (!secondary[i])
                    {
                      isGet = 0;
                    }
                  }
            }

            // Remember the value found in the cycle
            if (isGet)
            {
              cellValue = i;
            }
          }

        // Process the search result
        if (isGet)
        {
          // Process the new found value
            // Read the current value
            oldCellValue = newSquare.Matrix[rowId][columnId];
            // Write the new value
              // Write the value into the square
              newSquare.Matrix[rowId][columnId] = cellValue;
              // Mark the value in columns
              columns[cellValue][columnId] = Used;
              // Mark the value in rows
              rows[rowId][cellValue] = Used;
              // Mark the value in diagonals
              if (rowId == columnId)
              {
                primary[cellValue] = Used;
              }
              if (rowId == Rank - 1 - columnId)
              {
                secondary[cellValue] = Used;
              }
              // Mark the value in the history of cell values
              cellsHistory[rowId][columnId][cellValue] = Used;

            // Restore the previous value without clearing the history (because we are working with this cell)
            if (oldCellValue != Square::Empty)
            {
              // Restore the value into columns
              columns[oldCellValue][columnId] = Free;
              // Restore the value into rows
              rows[rowId][oldCellValue] = Free;
              // Restore the value into diagonals
              if (rowId == columnId)
              {
                primary[oldCellValue] = Free;
              }
              if (rowId == Rank - 1 - columnId)
              {
                secondary[oldCellValue] = Free;
              }
            }

            // Process the finish of the square generation
            if (cellId == cellsInPath - 1)
            {
              // Process the found square
              ProcessSquare();
            }
            else
            {
              // Step forward
              cellId++;
            }
        }
        else
        {
          // Process the fact of not-founding a new value in the cell (rowId; columnId)
            // Restore the previous value from the square into arrays 
              // Read the current value
              cellValue = newSquare.Matrix[rowId][columnId];
              // Restore the value into auxilary arrays
              if (cellValue != Square::Empty)
              {
                // Restore the value into columns
                columns[cellValue][columnId] = Free;
                // Restore the value into rows
                rows[rowId][cellValue] = Free;
                // Restore the value into diagonals
                if (rowId == columnId)
                {
                  primary[cellValue] = Free;
                }
                if (rowId == Rank - 1 - columnId)
                {
                  secondary[cellValue] = Free;
                }
                // Reset the cell of the square
                newSquare.Matrix[rowId][columnId] = Square::Empty;
                // Clear the history of the cell (rowId, columnId)
                for (int i = 0; i < Rank; i++)
                {
                  cellsHistory[rowId][columnId][i] = 1;
                }
              }

            // Step backward
            cellId--;
        }

        // Check the finish condition of search 
        if (keyValue == Square::Empty)
        {
          // Set the flag if the terminal value is "-1" which means we must leave the cell
          if (newSquare.Matrix[keyRowId][keyColumnId] == keyValue && cellId < 0)
          {
            stop = Yes;
          }
        }
        else
        {
          // Set the flag if the terminal value is other
          if (newSquare.Matrix[keyRowId][keyColumnId] == keyValue)
          {
            stop = Yes;
          }
        }
    }
    while (!stop);
  }
}


// Subscribe to the event of square generation
void Generator::Subscribe(MovePairSearch *search)
{
  subscriber = search;
}


// Unsubscribe from the event of square generation
void Generator::Unsubscribe()
{
  subscriber = 0;
}


// Process the square
void Generator::ProcessSquare()
{
  // Increase the counter of found squares
  squaresCount++;

  // Signal about square generation
  if (subscriber != 0)
  {
    subscriber->OnSquareGenerated(newSquare);
  }
}

