// Генератор Диагональных латинских квадратов

# include <iostream>
# include <fstream>
# include <iomanip>
# include <string>
# include <sstream>

# include "Generator.h"

using namespace std;

// Конструктор по умолчанию
Generator::Generator()
{
  // Сброс настроек
  Reset();

  // Задание текстовых констант
  generatorStateHeader = "# Workunits generator for Rake-search status";
}


// Конструктор копирования
Generator::Generator(Generator& source)
{
  CopyState(source);
}


// Генерация пакета заданий
void Generator::GenerateWorkunits(string start, string result, string checkpoint, string temp, int workunits)
{
  // Сброс состояния генератора
  Reset();

  // Инициализация генератора
  Initialize(start, result, checkpoint, temp, workunits);

  // Генерация пакета заданий
  Start();

  // Создание контрольной точки после завершения генерации
  CreateCheckpoint();
}


// Сброс всех значений внутренних структур
void Generator::Reset()
{
  // Сброс внутренних значений квадрата
  newSquare.Reset();

  // Сброс значений структур генерации квадратов
    // Сброс значений, соответствующих ключевой клетке
    keyRowId = Square::Empty;
    keyColumnId = Square::Empty;
    keyValue = Square::Empty;

    // Сброс значений, связанных с путём заполнения клеток
    for (int i = 0; i < CellsInPath; i++)
    {
      path[i][0] = Square::Empty;
      path[i][1] = Square::Empty;
    }

    // Сброс значений в векторах использования элементов на диагонали
    for (int i = 0; i < Rank; i++)
    {
      primary[i] = Free;
      secondary[i] = Free;
    }

    // Сброс значений в матрицах использования элементов в столбцах и строках
    for (int i = 0; i < Rank; i++)
    {
      for (int j = 0; j < Rank; j++)
      {
        columns[i][j] = Free;
        rows[i][j] = Free;
      }
    }

    // Сброс значений в кубе истории использования значений в клетках
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

    // Сброс координат обрабатываемой клетки
    rowId = Square::Empty;
    columnId = Square::Empty;

    // Сброс названий файлов
    checkpointFileName.clear();
    tempCheckpointFileName.clear();
    resultFileName.clear();

    // Сброс номера сгенерированного снимка и счётчиков заданий
    snapshotNumber = 0;
    workunitsCount = 0;
    workunitsToGenerate = 0;

    // Сброс адреса ключевой клетки
    workunitKeyRowId = 0;
    workunitKeyColumnId = 0;

    // Сброс флага инициализированности
    isInitialized = No;
}


// Инициализация генератора
void Generator::Initialize(string start, string result, string checkpoint, string temp, int workunits)
{
  fstream startFile;
  fstream checkpointFile;

  // Сбрасываем значения внутренних структур
  Reset();

  // Запоминаем имена файлов конфигурации, контрольных точек и результатов
  startParametersFileName = start;
  checkpointFileName = checkpoint;
  resultFileName = result;
  tempCheckpointFileName = temp;
  workunitsToGenerate = workunits;

  // Считываем настройки
  startFile.open(startParametersFileName.c_str(), std::ios_base::in);
  checkpointFile.open(checkpointFileName.c_str(), std::ios_base::in);

  if (checkpointFile.is_open())
  {
    // Считываем данные из файла контрольной точки
    Read(checkpointFile);
  }
  else
  {
    // Считываем данные из файла параметров запуска
    if (startFile.is_open())
    {
      Read(startFile);
    }
  }

  startFile.close();
  checkpointFile.close();
}

// Оператор записи состояния генератора
std::ostream& operator << (std::ostream& os, Generator& value)
{
  value.Write(os);

return os;
}


// Оператор считывания состояния генератора
std::istream& operator >> (std::istream& is, Generator& value)
{
  value.Read(is);

return is;
}


// Считывание состояние генератора из потока
void Generator::Read(std::istream& is)
{
  int rankToVerify;
  int result = Yes;
  string marker;

  // Сбрасываем флаг инициализированности
  isInitialized = No;

  // Поиск заголовка
  do
  {
    std::getline(is, marker);
  }
  while (marker != generatorStateHeader);

  // Считывание из потока ранга квадрата
  is >> rankToVerify;

  // Считывание данных поиска нужного нам ранга
  if (rankToVerify == Square::Rank)
  {
    // Считывание из потока квадрата
    is >> newSquare;

    // Считывание из потока пути обхода клеток
    for (int i = 0; i < CellsInPath; i++)
    {
      is >> path[i][0];
      is >> path[i][1];
    }

    // Считывание из потока информации о ключевой клетке
    is >> keyRowId;
    is >> keyColumnId;
    is >> keyValue;

    // Считывание информации об обрабатываемой клетке
    is >> rowId;
    is >> columnId;
    is >> cellId;

    // Считывание из потока информации о задействованных значениях и истории значений
      // Считывание информации о значениях на главной диагонали
      for (int i = 0; i < Rank; i++)
      {
        is >> primary[i];
      }

      // Считывание информации о значениях на побочной диагонали
      for (int i = 0; i < Rank; i++)
      {
        is >> secondary[i];
      }

      // Считывание информации о значениях в строках
      for (int i = 0; i < Rank; i++)
      {
        for (int j = 0; j < Rank; j++)
        {
          is >> rows[i][j];
        }
      }

      // Считывание информации о значениях в столбцах
      for (int i = 0; i < Rank; i++)
      {
        for (int j = 0; j < Rank; j++)
        {
          is >> columns[i][j];
        }
      }

      // Считывание информации об истории значений в клетках квадрата
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

    // Считываем номер последнего сгенерированного снимка
    is >> snapshotNumber;

    // Считывание строки и столбца ключевой клетки
    is >> workunitKeyRowId;
    is >> workunitKeyColumnId;

    // Выставляем флаг инициализированности
    isInitialized = Yes;
  }
}


// Запись состояния генератора в поток
void Generator::Write(std::ostream& os)
{
  // Запись заголовка
  os << generatorStateHeader << endl << endl;

  // Запись в поток ранга квадрата
  os << Square::Rank << endl;

  // Запись в поток квадрата
  os << newSquare;
  
  // Запись в поток пути обхода клеток
  for (int i = 0; i < CellsInPath; i++)
  {
    os << path[i][0] << " ";
    os << path[i][1] << " ";
    os << endl;
  }
  os << endl;

  // Запись в поток информации о ключевой клетке
  os << keyRowId << " " << keyColumnId << " " << keyValue << endl;

  // Запись информации о текущей клетке
  os << rowId << " " << columnId << " " << cellId  << endl;

  // Записываем пустую строку для удобства
  os << endl;

  // Запись информации о задействованных значениях и истории значений
    // Запись информации о значениях на главной диагонали
    for (int i = 0; i < Rank; i++)
    {
      os << primary[i] << " ";
    }
    os << endl;

    // Запись информации о значениях на побочной диагонали
    for (int i = 0; i < Rank; i++)
    {
      os << secondary[i] << " ";
    }
    os << endl;

    // Дополнительная пустая строка
    os << endl;

    // Запись информации о значениях в строках
    for (int i = 0; i < Rank; i++)
    {
      for (int j = 0; j < Rank; j++)
      {
        os << rows[i][j] << " ";
      }
      os << endl;
    }
    os << endl;

    // Запись информации о значениях в столбцах
    for (int i = 0; i < Rank; i++)
    {
      for (int j = 0; j < Rank; j++)
      {
        os << columns[i][j] << " ";
      }
      os << endl;
    }
    os << endl;

    // Запись информации об истории значений в клетках квадрата
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

  // Запись в поток информации о номере последнего сгенерированного снимка
  os << snapshotNumber << endl;
  os << endl;

  // Запис строки и столбца ключевой клетки для задания
  os << workunitKeyRowId << " " << workunitKeyColumnId << endl;
  os << endl;
}


// Оператор присваивания
Generator& Generator::operator = (Generator& value)
{
  CopyState(value);

return *this;
}


// Копирование состояния с заданного объекта
void Generator::CopyState(Generator& source)
{
  // Копируем переменные связанные с путём заполнения ячеек
  for (int i = 0; i < CellsInPath; i++)
  {
    path[i][0] = source.path[i][0];
    path[i][1] = source.path[i][1];
  }

  keyRowId = source.keyRowId;
  keyColumnId = source.keyColumnId;
  keyValue = source.keyValue;

  // Копируем массивы флагов использования значений
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

  // Копирование названий файлов
  startParametersFileName = source.startParametersFileName;
  resultFileName = source.resultFileName;
  checkpointFileName = source.checkpointFileName;
  tempCheckpointFileName = source.tempCheckpointFileName;

  // Копирование переменных текущего состояния
  isInitialized = source.isInitialized;
  snapshotNumber = source.snapshotNumber;
  rowId = source.rowId;
  columnId = source.columnId;
  cellId = source.cellId;

  // Копирование адресов текстовых констант
  generatorStateHeader = source.generatorStateHeader;
}


// Создание контрольной точки
void Generator::CreateCheckpoint()
{
  fstream newCheckpointFile;

  // Записываем настройки в новый файл контрольной точки
  newCheckpointFile.open(tempCheckpointFileName.c_str(), std::ios_base::out);

  if (newCheckpointFile.is_open())
  {
    Write(newCheckpointFile);
    newCheckpointFile.close();
    remove(checkpointFileName.c_str());
    rename(tempCheckpointFileName.c_str(), checkpointFileName.c_str());
  }
}


// Запуск генерации заданий
void Generator::Start()
{
  int isGet;         // Флаг получения нового значения для клетки
  int cellValue;     // Новое значение для клетки
  int oldCellValue;  // Старое значение, стоявшее в клетке

  int stop = 0;      // Флаг достижения окончания расчёта

  if (isInitialized == Yes)
  {
    // Подбор значений клеток квадрата
    do
    {
      // Подбор значения для очередной клетки квадрата
        // Считываем координаты клетки
        rowId = path[cellId][0];
        columnId = path[cellId][1];

        // Генерируем новое значение для клетки (rowId, columnId)
          // Сбрасываем значения переменных
          isGet = 0;
          cellValue = Square::Empty;

          // Подбираем значение для клетки
          for (int i = 0; i < Rank && !isGet; i++)
          {
            // Проверяем значение i на возможность записи в клетку (rowId, columnId)
            if (columns[i][columnId] && rows[rowId][i] && cellsHistory[rowId][columnId][i])
            {
              // Значение не занято в столбцах и строках, но надо ещё проверить диагонали
                // Выставляем флаг, который, возможно, будет сброшен диагональной проверкой
                isGet = 1;
                // Проверяем значение - не попадалось ли оно на диагоналях
                  // Проверка первой диагонали
                  if(columnId == rowId)
                  {
                    if (!primary[i])
                    {
                      isGet = 0;
                    }
                  }

                  // Проверка второй диагонали
                  if (rowId == Rank - 1 - columnId)
                  {
                    if (!secondary[i])
                    {
                      isGet = 0;
                    }
                  }
            }

            // Запоминание значения, найденного в цикле
            if (isGet)
            {
              cellValue = i;
            }
          }

        // Обработка результата поиска
        if (isGet)
        {
          // Обработка найденного нового значения

            // Считывание текущего значения
              oldCellValue = newSquare.Matrix[rowId][columnId];
            
            // Запись нового значения
              // Записываем значение в квадрат
              newSquare.Matrix[rowId][columnId] = cellValue;
              // Отмечаем значение в столбцах
              columns[cellValue][columnId] = Used;
              // Отмечаем значение в строках
              rows[rowId][cellValue] = Used;
              // Отмечаем значение в диагоналях
              if (rowId == columnId)
              {
                primary[cellValue] = Used;
              }
              if (rowId == Rank - 1 - columnId)
              {
                secondary[cellValue] = Used;
              }
              // Отмечаем значение в истории значений клетки
              cellsHistory[rowId][columnId][cellValue] = Used;

            // Возвращение предыдущего значения без зачистки истории 
            // (так как мы работаем с этой клеткой)
            if (oldCellValue != Square::Empty)
            {
              // Возвращаем значение в столбцы
              columns[oldCellValue][columnId] = Free;
              // Возвращаем значение в строки
              rows[rowId][oldCellValue] = Free;
              // Возвращаем значение в диагонали
              if (rowId == columnId)
              {
                primary[oldCellValue] = Free;
              }
              if (rowId == Rank - 1 - columnId)
              {
                secondary[oldCellValue] = Free;
              }
            }

            // Обработка окончания формирования квадрата
            if (cellId == CellsInPath - 1)
            {
              // Генерируем задание
              ProcessWorkunit();
            }
            else
            {
              // Делаем шаг вперёд
              cellId++;
            }
        }
        else
        {
          // Обработка факта ненахождения нового значения в клетке (rowId; columnId)
            // Возвращаем текущее значение из квадрата в массивы
              // Считываем текущее значение
              cellValue = newSquare.Matrix[rowId][columnId];
              // Возвращаем значение в служебные массивы
              if (cellValue != Square::Empty)
              {
                // Возвращаем значение в столбцы
                columns[cellValue][columnId] = Free;
                // Возвращаем значение в строки
                rows[rowId][cellValue] = Free;
                // Возвращаем значение в диагонали
                if (rowId == columnId)
                {
                  primary[cellValue] = Free;
                }
                if (rowId == Rank - 1 - columnId)
                {
                  secondary[cellValue] = Free;
                }
                // Сбрасываем клетку квадрата
                newSquare.Matrix[rowId][columnId] = Square::Empty;
                // Зачищаем историю клетки (rowId, columnId)
                for (int i = 0; i < Rank; i++)
                {
                  cellsHistory[rowId][columnId][i] = 1;
                }
              }

            // Делаем шаг назад
            cellId--;
        }

        // Проверяем условия окончания поиска
        if (keyValue == Square::Empty)
        {
          // Выставление флага при завершаеющем значении "-1" при котором производится уход из клетки
          if (newSquare.Matrix[keyRowId][keyColumnId] == keyValue && cellId < 0)
          {
            stop = Yes;
          }
        }
        else
        {
          // Выставление флага при обычном завершаеющем значении
          if (newSquare.Matrix[keyRowId][keyColumnId] == keyValue)
          {
            stop = Yes;
          }
        }

        // Проверяем условие окончания генерации заданий
        if (workunitsCount == workunitsToGenerate)
        {
          stop = Yes;
        }
    }
    while (!stop);
  }
}


// Обработка квадрата
void Generator::ProcessWorkunit()
{
  switch (Rank)
  {
    case 8:
      ProcessWorkunitR8();
    break;

    case 9:
      ProcessWorkunitR9();
    break;
  }
}


// Формирование очередного задания для квадратов 8-го ранга
void Generator::ProcessWorkunitR8()
{
  fstream workunitFile;
  string workunitFileName;
  std::stringstream wuNameBuilder;
  int workunitNumber;
  int workunitKeyCellValue;

  // Увеличиваем счётчик номера снимка
  snapshotNumber++;

  // Формируем задание
    // Определение номера задания и стопорного значения ключевой клетки
    workunitNumber = snapshotNumber;
    workunitKeyCellValue = newSquare.Matrix[workunitKeyRowId][workunitKeyColumnId];

    // Формируем название файла
    wuNameBuilder << "wu_" << std::setw(6) << std::setfill('0') << workunitNumber << ".txt";
    workunitFileName = wuNameBuilder.str();

    // Формируем файл с заданием
    cout << "Workunit " << workunitFileName << ":" << endl;
    cout << newSquare << endl;

    // Открываем файл
    workunitFile.open(workunitFileName.c_str(), std::ios_base::out);
  
    // Записываем заголовки и квадрат, с которого надо начинать генерацию
    workunitFile << "# Move search of pairs OLDS status" << endl
          << endl
          << "# Generation of DLS status" << endl
          << endl
          << Rank << endl
          << endl
          << newSquare
          << endl;
  
    // Записываем число клеток в пути
    workunitFile << (Rank - 1)*Rank - CellsInPath << endl
           << endl;

    // Записываем путь переходов по клеткам
    workunitFile << "7 0" << endl << "6 1" << endl << "5 2" << endl << "4 3" << endl << "3 4" << endl << "2 5" << endl << "1 6" << endl
           << "1 0" << endl << "1 2" << endl << "1 3" << endl << "1 4" << endl << "1 5" << endl << "1 7" << endl
           << "2 0" << endl << "2 1" << endl << "2 3" << endl << "2 4" << endl << "2 6" << endl << "2 7" << endl
           << "3 0" << endl << "3 1" << endl << "3 2" << endl << "3 5" << endl << "3 6" << endl << "3 7" << endl
           << "4 0" << endl << "4 1" << endl << "4 2" << endl << "4 5" << endl << "4 6" << endl << "4 7" << endl
           << "5 0" << endl << "5 1" << endl << "5 3" << endl << "5 4" << endl << "5 6" << endl << "5 7" << endl
           << "6 0" << endl << "6 2" << endl << "6 3" << endl << "6 4" << endl << "6 5" << endl << "6 7" << endl
           << "7 1" << endl << "7 2" << endl << "7 3" << endl << "7 4" << endl << "7 5" << endl << "7 6" << endl << endl;

    // Запись координаты ключевой клетки и её значение
    workunitFile << workunitKeyRowId << " " << workunitKeyColumnId << " " << workunitKeyCellValue << endl;
    // Запись координат текущей клетки и идентификатора текущей клетки в пути
    workunitFile << workunitKeyRowId << " " << workunitKeyColumnId << " " << 0 << endl;
    workunitFile << endl;

    // Запись состояний диагоналей, строк, столбцов и истории использования значений (копия из Write!)
      // Запись информации о значениях на главной диагонали
      for (int i = 0; i < Rank; i++)
      {
        workunitFile << primary[i] << " ";
      }
      workunitFile << endl;

      // Запись информации о значениях на побочной диагонали
      for (int i = 0; i < Rank; i++)
      {
        workunitFile << secondary[i] << " ";
      }
      workunitFile << endl;

      // Дополнительная пустая строка
      workunitFile << endl;

      // Запись информации о значениях в строках
      for (int i = 0; i < Rank; i++)
      {
        for (int j = 0; j < Rank; j++)
        {
          workunitFile << rows[i][j] << " ";
        }
        workunitFile << endl;
      }
      workunitFile << endl;

      // Запись информации о значениях в столбцах
      for (int i = 0; i < Rank; i++)
      {
        for (int j = 0; j < Rank; j++)
        {
          workunitFile << columns[i][j] << " ";
        }
        workunitFile << endl;
      }
      workunitFile << endl;

      // Запись информации об истории значений в клетках квадрата
      for (int h = 0; h < Rank; h++)
      {
        for (int i = 0; i < Rank; i++)
        {
          for (int j = 0; j < Rank; j++)
          {
            workunitFile << cellsHistory[i][j][h] << " ";
          }
          workunitFile << endl;
        }
        workunitFile << endl;
      }
      workunitFile << endl;

      // Запись информации о числе сгенерированных WU
      workunitFile << 0 << endl
             << endl;

    // Запись данных компоненты перетасовки строк
    workunitFile << "# Move search component status" << endl
           << endl
           << "0 0 0" << endl
           << "0 0" << endl;

    // Закрываем файл
    workunitFile.close();

  // Увеличиваем счётчик сгенерированных заданий 
  workunitsCount++;
}


// Формирование очередного задания для квадратов 9-го ранга
void Generator::ProcessWorkunitR9()
{
  fstream workunitFile;
  string workunitFileName;
  std::stringstream wuNameBuilder;
  int workunitNumber;
  int workunitKeyCellValue;

  // Увеличиваем счётчик номера снимка
  snapshotNumber++;

  // Формируем задание
    // Определение номера задания и стопорного значения ключевой клетки
    workunitNumber = snapshotNumber;
    workunitKeyCellValue = newSquare.Matrix[workunitKeyRowId][workunitKeyColumnId];

    // Формируем название файла
    wuNameBuilder << "wu_" << std::setw(9) << std::setfill('0') << workunitNumber << ".txt";
    workunitFileName = wuNameBuilder.str();

    // Формируем файл с заданием
    cout << "Workunit " << workunitFileName << ":" << endl;
    cout << newSquare << endl;

    // Открываем файл
    workunitFile.open(workunitFileName.c_str(), std::ios_base::out);
  
    // Записываем заголовки и квадрат, с которого надо начинать генерацию
    workunitFile << "# Move search of pairs OLDS status" << endl
          << endl
          << "# Generation of DLS status" << endl
          << endl
          << Rank << endl
          << endl
          << newSquare
          << endl;
  
    // Записываем число клеток в пути
    workunitFile << (Rank - 1)*Rank - CellsInPath << endl
           << endl;

    // Записываем путь переходов по клеткам
    workunitFile << "1 2" << endl << "1 3" << endl << "1 4" << endl << "1 5" << endl << "1 6" << endl << "1 8" << endl
           << "2 0" << endl << "2 1" << endl << "2 3" << endl << "2 4" << endl << "2 5" << endl << "2 7" << endl << "2 8" << endl
           << "3 0" << endl << "3 1" << endl << "3 2" << endl << "3 4" << endl << "3 6" << endl << "3 7" << endl << "3 8" << endl
           << "4 0" << endl << "4 1" << endl << "4 2" << endl << "4 3" << endl << "4 5" << endl << "4 6" << endl << "4 7" << endl << "4 8" << endl
           << "5 0" << endl << "5 1" << endl << "5 2" << endl << "5 4" << endl << "5 6" << endl << "5 7" << endl << "5 8" << endl
           << "6 0" << endl << "6 1" << endl << "6 3" << endl << "6 4" << endl << "6 5" << endl << "6 7" << endl << "6 8" << endl
           << "7 0" << endl << "7 2" << endl << "7 3" << endl << "7 4" << endl << "7 5" << endl << "7 6" << endl << "7 8" << endl
           << "8 1" << endl << "8 2" << endl << "8 3" << endl << "8 4" << endl << "8 5" << endl << "8 6" << endl << "8 7" << endl << endl;

    // Запись координаты ключевой клетки и её значение
    workunitFile << workunitKeyRowId << " " << workunitKeyColumnId << " " << workunitKeyCellValue << endl;
    // Запись координат текущей клетки и идентификатора текущей клетки в пути
    workunitFile << workunitKeyRowId << " " << workunitKeyColumnId << " " << 0 << endl;
    workunitFile << endl;

    // Запись состояний диагоналей, строк, столбцов и истории использования значений (копия из Write!)
      // Запись информации о значениях на главной диагонали
      for (int i = 0; i < Rank; i++)
      {
        workunitFile << primary[i] << " ";
      }
      workunitFile << endl;

      // Запись информации о значениях на побочной диагонали
      for (int i = 0; i < Rank; i++)
      {
        workunitFile << secondary[i] << " ";
      }
      workunitFile << endl;

      // Дополнительная пустая строка
      workunitFile << endl;

      // Запись информации о значениях в строках
      for (int i = 0; i < Rank; i++)
      {
        for (int j = 0; j < Rank; j++)
        {
          workunitFile << rows[i][j] << " ";
        }
        workunitFile << endl;
      }
      workunitFile << endl;

      // Запись информации о значениях в столбцах
      for (int i = 0; i < Rank; i++)
      {
        for (int j = 0; j < Rank; j++)
        {
          workunitFile << columns[i][j] << " ";
        }
        workunitFile << endl;
      }
      workunitFile << endl;

      // Запись информации об истории значений в клетках квадрата
      for (int h = 0; h < Rank; h++)
      {
        for (int i = 0; i < Rank; i++)
        {
          for (int j = 0; j < Rank; j++)
          {
            workunitFile << cellsHistory[i][j][h] << " ";
          }
          workunitFile << endl;
        }
        workunitFile << endl;
      }
      workunitFile << endl;

      // Запись информации о числе сгенерированных WU
      workunitFile << 0 << endl
             << endl;

    // Запись данных компоненты перетасовки строк
    workunitFile << "# Move search component status" << endl
           << endl
           << "0 0 0" << endl
           << "0 0" << endl;

    // Закрываем файл
    workunitFile.close();

  // Увеличиваем счётчик сгенерированных заданий 
  workunitsCount++;
}
