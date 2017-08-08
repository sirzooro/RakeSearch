// Генератор Диагональных латинских квадратов

# include "Generator.h"

using namespace std;

// Конструктор по умолчанию
Generator::Generator()
{
  // Сброс настроек
  Reset();

  // Задание текстовых констант
  generatorStateHeader = "# Generation of DLS status";
}


// Конструктор копирования
Generator::Generator(Generator& source)
{
  CopyState(source);
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
    for (int i = 0; i < MaxCellsInPath; i++)
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

    // Сброс числа сгенерированных квадратов
    squaresCount = 0;

    // Сброс флага инициализированности
    isInitialized = No;
}


// Инициализация генератора
void Generator::Initialize(string start, string result, string checkpoint, string temp)
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

    // Считывание числа клеток в пути обхода
    is >> cellsInPath;

    // Считывание из потока пути обхода клеток
    for (int i = 0; i < cellsInPath; i++)
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

    // Считываем число сгенерированных квадратов
    is >> squaresCount;

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

  // Запись числа клеток в пути обхода
  os << cellsInPath << endl;
  os << endl;
  
  // Запись в поток пути обхода клеток
  for (int i = 0; i < cellsInPath; i++)
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

  // Запись в поток информации о числе сгенерированных квадратов
  os << squaresCount << endl;
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
  for (int i = 0; i < cellsInPath; i++)
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
  squaresCount = source.squaresCount;
  rowId = source.rowId;
  columnId = source.columnId;
  cellId = source.cellId;

  // Копирование адресов текстовых констант
  generatorStateHeader = source.generatorStateHeader;
}


// Заданием имен файлов параметров и контрольной точки
void Generator::SetFileNames(string start, string result, string checkpoint, string temp)
{
  startParametersFileName = start;
  resultFileName = result;
  checkpointFileName = checkpoint;
  tempCheckpointFileName = temp;
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


// Запуск генерации квадратов
void Generator::Start()
{
  int isGet;        // Флаг получения нового значения для клетки
  int cellValue;    // Новое значение для клетки
  int oldCellValue; // Старое значение, стоявшее в клетке

  int stop = 0;     // Флаг достижения окончания расчёта

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

            // Возвращение предыдущего значения без зачистки истории (так как мы работаем с этой клеткой)
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
            if (cellId == cellsInPath - 1)
            {
              // Обрабатываем найденный квадрат
              ProcessSquare();
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

        // Проверяем условие окончания поиска
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
    }
    while (!stop);
  }
}


// Обработка квадрата
void Generator::ProcessSquare()
{
  // Увеличиваем счётчик найденных квадратов
  squaresCount++;

  // Генерируем событие
  this->SquareGenerated(newSquare);
}

