# include <iostream>
# include <fstream>
# include <string>
#include <sys/stat.h>

# include "MovePairSearch.h"
# include "PairSearch.h"

using namespace std;

// Проверка существования файла
inline bool file_exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

// Выполнение вычислений
int Compute(string wu_filename, string result_filename)
{
  string localWorkunit;
  string localResult;
  string localCheckpoint;
  string pathLocalWorkunit;
  string pathLocalResult;
  string pathLocalCheckpoint;

  string initStartFileName;
  string initResultFileName;
  string initCheckpointFileName;
  string initTempCheckpointFileName;

  MovePairSearch search;

  {
    // Проверка наличия файла задания, контрольной точки, результата
    localWorkunit = wu_filename; 
    localResult = result_filename; 
    localCheckpoint = "checkpoint_" + wu_filename;

    pathLocalWorkunit = localWorkunit;
    pathLocalResult = localResult;
    pathLocalCheckpoint = localCheckpoint;

    // Результат существует, удалить результат и контрольные точки
    if (file_exists(localResult))
    {
      // Отправка результата и удаление файла с заданием и контрольной точкой
        // Удаление файла с заданием
        std::cout << "Remove a workunit file: " << localWorkunit << endl;
        remove(pathLocalWorkunit.c_str());

        // Удаление файла с контрольной точкой
        std::cout << "Remove a checkpoint file: " << localCheckpoint << endl;
        remove(pathLocalCheckpoint.c_str());

        std::cout << "Result removed" << endl;
        return 0;
    }

    // Запуск вычислений с контрольной точки
    if (file_exists(localCheckpoint) && !file_exists(localResult))
    {
      // Проверка наличия файла с заданием
      if (file_exists(localWorkunit))
      {
        // Запускаем расчёт
        initStartFileName  = localWorkunit;
        initResultFileName = localResult;
        
        std::cout << "Start from checkpoint of workunit " << localWorkunit << endl;

        search.InitializeMoveSearch(initStartFileName, initResultFileName, 
                      initCheckpointFileName, initTempCheckpointFileName);
        search.StartMoveSearch();
      }
      else
      {
        std::cout << "Error: delected a checkpoint file " << localCheckpoint << " without workunit file!" << endl;
        return -1;
      }
    }

    // Запуск вычислений с файла задания
    if (!file_exists(localCheckpoint) && !file_exists(localResult) 
                                    && file_exists(localWorkunit))
    {
      // Запуск вычислений с файла задания, присутствующего без файлов 
      // контрольной точки и результата
      initStartFileName  = localWorkunit;
      initResultFileName = localResult;
      initCheckpointFileName = "checkpoint_" + localWorkunit;      
      initTempCheckpointFileName = "checkpoint_new_" + localWorkunit;

      std::cout << "Start from workunit file " << localWorkunit << endl;
      search.InitializeMoveSearch(initStartFileName, initResultFileName, 
                      initCheckpointFileName, initTempCheckpointFileName);
      search.StartMoveSearch();
    }
  }

  return 0;
}


int main(int argumentsCount, char* argumentsValues[])
{
  string wu_filename;
  string result_filename;

  if (argumentsCount == 3)
  {
    wu_filename = argumentsValues[1];
    result_filename = argumentsValues[2];
    Compute(wu_filename, result_filename);
  }
  else
  {
    std::cout << "Please, specify WU filename and result filename." << endl << endl;
    return -1;
  }

  //cout << "Press any key to exit ... " << endl;
  //cin.get();

  return 0;
}

/*
Sample of Pair search

  int aMatrix[9][9] =
  {
    {0, 1, 2, 3, 4, 5, 6, 7, 8},
    {4, 2, 7, 6, 8, 1, 3, 5, 0},
    {3, 5, 1, 0, 7, 8, 4, 6, 2},
    {6, 3, 8, 4, 1, 7, 0, 2, 5},
    {5, 6, 0, 7, 3, 2, 8, 4, 1},
    {7, 8, 4, 1, 5, 6, 2, 0, 3},
    {8, 7, 6, 2, 0, 3, 5, 1, 4},
    {1, 0, 3, 5, 2, 4, 7, 8, 6},
    {2, 4, 5, 8, 6, 0, 1, 3, 7}
  };

  Square a(aMatrix);
  PairSearch search;

  search.Initialize("start_parameters.txt", "result.txt", "checkpoint.txt", "checkpoint_new.txt");
  search.OnSquareGenerated(a);

*/
