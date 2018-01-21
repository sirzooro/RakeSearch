// DLS generator

#include "Generator.h"
#include "MovePairSearch.h"
#include <string.h>
#include <type_traits>

#define GetBit(bitfield, bitno) ((bitfield) & (1u << (bitno)))
#define SetBit(bitfield, bitno) ((bitfield) |= (1u << (bitno)))
#define ClearBit(bitfield, bitno) ((bitfield) &= ~(1u << (bitno)))

#define AllBitsMask(numbits) ((1u << (numbits)) - 1)

// Used = 0, Free = 1, Code now uses bits, so dedicated macros would be helpful.
#define SetUsed(bitfield, bitno) ClearBit(bitfield, bitno)
#define SetFree(bitfield, bitno) SetBit(bitfield, bitno)

#define IsUsed(bitfield, bitno) (0 == GetBit(bitfield, bitno))
#define IsFree(bitfield, bitno) (0 != GetBit(bitfield, bitno))

#define AllFree AllBitsMask(Rank)

#define GetBit01(bitfield, bitno) (GetBit(bitfield, bitno) ? 1 : 0)

#define ffs __builtin_ffs


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
Generator::Generator(const Generator& source)
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
    primary = AllFree;
    secondary = AllFree;

    // Reset the values in the matrices of rows/columns elements usage
    for (int i = 0; i < Rank; i++)
    {
      columns[i] = AllFree;
      rows[i] = AllFree;
    }

    // Reset the values in the cube of the history of cell values usage
    for (int i = 0; i < Rank; i++)
    {
      for (int j = 0; j < Rank; j++)
      {
        cellsHistory[i][j] = AllFree;
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
void Generator::Initialize(const string& start, const string& result, const string& checkpoint, const string& temp)
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
  string marker;
  int val;

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
      primary = 0;
      for (int i = 0; i < Rank; i++)
      {
        is >> val;
        if (val)
          SetBit(primary, i);
      }

      // Read the information about the secondary diagonal values
      secondary = 0;
      for (int i = 0; i < Rank; i++)
      {
        is >> val;
        if (val)
          SetBit(secondary, i);
      }

      // Read the information about values in rows
      memset(rows, 0, sizeof(rows));
      for (int i = 0; i < Rank; i++)
      {
        for (int j = 0; j < Rank; j++)
        {
          is >> val;
          if (val)
            SetBit(rows[i], j);
        }
      }

      // Read the information about values in columns
      memset(columns, 0, sizeof(columns));
      for (int i = 0; i < Rank; i++)
      {
        for (int j = 0; j < Rank; j++)
        {
          is >> val;
          if (val)
            SetBit(columns[j], i);
        }
      }

      // Read the information about the history of values in square cells
      memset(cellsHistory, 0, sizeof(cellsHistory));
      for (int h = 0; h < Rank; h++)
      {
        for (int i = 0; i < Rank; i++)
        {
          for (int j = 0; j < Rank; j++)
          {
            is >> val;
            if (val)
              SetBit(cellsHistory[i][j], h);
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
      os << GetBit01(primary, i) << " ";
    }
    os << endl;

    // Write the information about secondary diagonal values
    for (int i = 0; i < Rank; i++)
    {
      os << GetBit01(secondary, i) << " ";
    }
    os << endl;

    // Additional empty line
    os << endl;

    // Write the information about rows values
    for (int i = 0; i < Rank; i++)
    {
      for (int j = 0; j < Rank; j++)
      {
        os << GetBit01(rows[i], j) << " ";
      }
      os << endl;
    }
    os << endl;

    // Write the information about columns values
    for (int i = 0; i < Rank; i++)
    {
      for (int j = 0; j < Rank; j++)
      {
        os << GetBit01(columns[j], i) << " ";
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
          os << GetBit01(cellsHistory[i][j], h) << " ";
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
Generator& Generator::operator = (const Generator& value)
{
  CopyState(value);

return *this;
}


// Copy the state from the given object
void Generator::CopyState(const Generator& source)
{
  // Copy variables connected with the path of cells bypassing
  memcpy(path, source.path, sizeof(path[0][0]) * 2 * cellsInPath);

  keyRowId = source.keyRowId;
  keyColumnId = source.keyColumnId;
  keyValue = source.keyValue;

  // Copy the flag arrays of used values
  primary = source.primary;
  secondary = source.secondary;

  memcpy(columns, source.columns, sizeof(columns));
  memcpy(rows, source.rows, sizeof(rows));

  memcpy(cellsHistory, source.cellsHistory, sizeof(cellsHistory));

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
void Generator::SetFileNames(const string& start, const string& result, const string& checkpoint, const string& temp)
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
  // Check value of keyValue and pass result as a type to StartImpl
  if (keyValue == Square::Empty)
    StartImpl<true_type>();
  else
    StartImpl<false_type>();
}

// Actual implementation of the squares generation
// Note: values on diagonal are preset in WU, so corresponding parts of code are commented out.
// It turned out that it was quite costly to have instructions which were doing nothing.
template<typename IsKeyValueEmpty>
inline void Generator::StartImpl()
{
  int cellValue;    // New value for the cell
  int oldCellValue; // Old value from the cell

  // Create constant copies of used fields to speedup calculations
  const int_fast32_t cellsInPath = this->cellsInPath;
  const int keyValue = this->keyValue;
  const int_fast32_t keyRowId = this->keyRowId;
  const int_fast32_t keyColumnId = this->keyColumnId;

  // Use registers for local variables instead of memory
  int_fast32_t rowId, columnId;
  int_fast32_t cellId = this->cellId;

  // Checkpoint may be written after new ODLS is created only.
  // Class members moved to registers above are constant in checkpoint
  // file, so they can be set to proper values here.
  this->rowId = path[cellsInPath - 1][0];
  this->columnId = path[cellsInPath - 1][1];
  this->cellId = cellsInPath - 1;

  if (isInitialized == Yes)
  {
    // Selection of the cells values
    while(1)
    {
      // Selection of the value for the next cell
        // Read coordinates of the cell
        rowId = path[cellId][0];
        columnId = path[cellId][1];

        // Generate new value for the cell (rowId, columnId)
          // Select the value for the cell
          // Check the i value for possibility to be written into the cell (rowId, columnId)
          cellValue = columns[columnId] & rows[rowId] & cellsHistory[rowId][columnId];

          // Test the value: has it been used in diagonals
          // Test the main diagonal
          /*if(columnId == rowId)
          {
            cellValue &= primary;
          }

          // Test the secondary diagonal
          if (rowId == Rank - 1 - columnId)
          {
            cellValue &= secondary;
          }*/

        // Process the search result
        if (cellValue)
        {
          // Extract lowest bit set
          int bits = (-cellValue) & cellValue;
          // Mark the value in the history of cell values
          cellsHistory[rowId][columnId] ^= bits;

          // Read the current value
          oldCellValue = newSquare.Matrix[rowId][columnId];
          // If it is there, set corresponding bit it bits to restore the previous value too.
          // Both bits will be swapped using XOR instructions.
          // History is not cleared, because we are working with this cell.
          if (oldCellValue != Square::Empty)
            bits |= 1 << oldCellValue;

          // Mark/restore the value in columns
          columns[columnId] ^= bits;
          // Mark/restore the value in rows
          rows[rowId] ^= bits;
          // Mark/restore the value in diagonals
          /*if (rowId == columnId)
          {
            primary ^= bits;
          }
          if (rowId == Rank - 1 - columnId)
          {
            secondary ^= bits;
          }*/
          // Write the value into the square
          newSquare.Matrix[rowId][columnId] = __builtin_ctz(cellValue);

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
                SetFree(columns[columnId], cellValue);
                // Restore the value into rows
                SetFree(rows[rowId], cellValue);
                // Restore the value into diagonals
                /*if (rowId == columnId)
                {
                  SetFree(primary, cellValue);
                }
                if (rowId == Rank - 1 - columnId)
                {
                  SetFree(secondary, cellValue);
                }*/
                // Reset the cell of the square
                newSquare.Matrix[rowId][columnId] = Square::Empty;
                // Clear the history of the cell (rowId, columnId)
                cellsHistory[rowId][columnId] = AllBitsMask(Rank);
              }

            // Step backward
            cellId--;

            // Check the finish condition of search
            if (IsKeyValueEmpty::value)
            {
              // Set the flag if the terminal value is "-1" which means we must leave the cell
              if (cellId < 0 && newSquare.Matrix[keyRowId][keyColumnId] == Square::Empty)
              {
                break;
              }
            }
        }

        // Check the finish condition of search
        if (!IsKeyValueEmpty::value)
        {
          // Set the flag if the terminal value is other
          if (newSquare.Matrix[keyRowId][keyColumnId] == keyValue)
          {
            break;
          }
        }
    }
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
