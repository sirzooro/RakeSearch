# include <iostream>
# include <fstream>
# include <string>
# include <sys/stat.h>
# include "boinc_api.h"

// g++ -c Square.cpp -o Square.o
// g++ -c Generator.cpp -o Generator.o
// g++ -c PairSearch.cpp -o PairSearch.o
// g++ -c MovePairSearch.cpp -o MovePairSearch.o
// g++ -c main.cpp -o main.o -I/usr/include/boinc
// g++ -o all *.o -L/usr/lib64 -lboinc_api

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
  string wu_filename = "workunit.txt";
  string result_filename = "result.txt";

  boinc_init();
  
  Compute(wu_filename, result_filename);

  boinc_finish(0);

  return 0;
}

