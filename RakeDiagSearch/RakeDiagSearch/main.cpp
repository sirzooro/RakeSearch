# include <iostream>
# include <fstream>
# include <string>
# include <sys/stat.h>
# include "boinc_api.h"

# include "MovePairSearch.h"

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
  string resolved_in_name;  // Переменные для работы с логическими
  string resolved_out_name; // и физическими именами файлов в BOINC

  int retval;

  boinc_init(); // Инициализация BOINC API для однопоточного приложения

  boinc_set_min_checkpoint_period(600); // Минимальное число секунд между записью контрольных точек
  
  // Преобразовать логическое имя файла в физическое.
  // Мы делаем это на верхнем уровне, передавая дальше уже преобразованные имена.
  retval = boinc_resolve_filename_s(wu_filename.c_str(), resolved_in_name);
  if (retval) { cout << "can't resolve IN filename!" << endl; boinc_finish(-1); return 0;}

  retval = boinc_resolve_filename_s(result_filename.c_str(), resolved_out_name);
  if (retval) { cout << "can't resolve OUT filename" << endl; boinc_finish(-1); return 0;}

  boinc_fraction_done(0.0); // Сообщить клиенту BOINC о доле выполнения задания

  // Запустить расчет
  Compute(resolved_in_name, resolved_out_name);

  boinc_fraction_done(1.0); // Сообщить клиенту BOINC о доле выполнения задания

  boinc_finish(0); // Сообщить клиенту BOINC о статусе завершения расчета
                   // (не делает return)
                   // Если нужно вывести сообщение пользователю, используем функцию
                   // boinc_finish_message(int status, const char* msg, bool is_notice); 
                   // If is_notice is true, the message will be shown as a notice 
                   // in the GUI (works with 7.5+ clients; for others, no message 
                   //  will be shown). 
  return 0;
}

