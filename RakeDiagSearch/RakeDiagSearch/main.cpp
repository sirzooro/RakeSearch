# include <iostream>
# include <fstream>
# include <string>
# include <sys/stat.h>
# include "boinc_api.h"

# include "MovePairSearch.h"

using namespace std;

static const bool isDebug = false;

// Проверка существования файла
inline bool file_exists (const std::string& name)
{
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

// Выполнение вычислений
int Compute(string wu_filename, string result_filename)
{
  string localWorkunit;
  string localResult;
  string localCheckpoint; 
  string localTmpCheckpoint;

  string initStartFileName;
  string initResultFileName;
  string initCheckpointFileName;
  string initTmpCheckpointFileName;

  MovePairSearch search;

    // Проверка наличия файла задания, контрольной точки, результата
    localWorkunit = wu_filename; 
    localResult = result_filename; 
    localCheckpoint = "checkpoint.txt";
    localTmpCheckpoint = "tmp_checkpoint.txt";

    // Запуск вычислений с контрольной точки
    if (file_exists(localCheckpoint))
    {
      // Проверка наличия файла с заданием
      if (file_exists(localWorkunit))
      {
        // Запускаем расчёт
        initStartFileName  = localWorkunit;
        initResultFileName = localResult;
        initCheckpointFileName = localCheckpoint;       
        initTmpCheckpointFileName = localTmpCheckpoint;
 
        if(isDebug) cout << "Start from checkpoint of workunit " << localWorkunit << endl;

        search.InitializeMoveSearch(initStartFileName, initResultFileName, 
                      initCheckpointFileName, initTmpCheckpointFileName);
        search.StartMoveSearch();
      }
      else
      {
        cerr << "Error: detected a checkpoint file " << localCheckpoint;
        cerr << " without workunit file!" << endl;
        return -1;
      }
    }

    // Запуск вычислений с файла задания
    if (!file_exists(localCheckpoint) && file_exists(localWorkunit))
    {
      // Запуск вычислений с файла задания, 
      // присутствующего без файлов 
      // контрольной точки и результата
      initStartFileName  = localWorkunit;
      initResultFileName = localResult;
      initCheckpointFileName = localCheckpoint;      
      initTmpCheckpointFileName = localTmpCheckpoint;

      if(isDebug) cout << "Start from workunit file " << localWorkunit << endl;
      search.InitializeMoveSearch(initStartFileName, initResultFileName, 
                      initCheckpointFileName, initTmpCheckpointFileName);
      search.StartMoveSearch();
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

  boinc_init(); // Инициализировать BOINC API для однопоточного приложения
  // Установить минимальное число секунд между записью контрольных точек 
  boinc_set_min_checkpoint_period(60); 
  
  // Преобразовать логическое имя файла в физическое.
  // Мы делаем это на верхнем уровне, передавая дальше уже преобразованные имена.
  retval = boinc_resolve_filename_s(wu_filename.c_str(), resolved_in_name);
  if (retval) 
  { 
    if(isDebug) cerr << "can't resolve IN filename!" << endl; 
    boinc_finish(retval); return 0;
  }
  retval = boinc_resolve_filename_s(result_filename.c_str(), resolved_out_name);
  if (retval) 
  { 
    if(isDebug) cerr << "can't resolve OUT filename" << endl; 
    boinc_finish(retval); return 0;
  }

  boinc_fraction_done(0.0); // Сообщить клиенту BOINC о доле выполнения задания
  // Запустить расчет
  retval = Compute(resolved_in_name, resolved_out_name);
  boinc_fraction_done(1.0); // Сообщить клиенту BOINC о доле выполнения задания

  // Сообщить клиенту BOINC о статусе завершения расчета (не делает return)
  boinc_finish(retval); 
                   // Если нужно вывести сообщение пользователю, используем функцию
                   // boinc_finish_message(int status, const char* msg, bool is_notice); 
                   // If is_notice is true, the message will be shown as a notice 
                   // in the GUI (works with 7.5+ clients; for others, no message 
                   //  will be shown). 
  return 0;
}