// Поиск пар диагональных латинских квадратов методом перетасовки строк

# include "MovePairSearch.h"
# include "boinc_api.h"

// Конструктор по умолчанию
MovePairSearch::MovePairSearch()
{
  Reset();
}

// Сброс настроек поиска
void MovePairSearch::Reset()
{
  // Сброс настроек генератора ДЛК
  squareAGenerator.Reset();

  // Сброс значений отдельного поиска
  ClearBeforeNextSearch();

  // Сброс значений глобальных счётчиков
  totalPairsCount = 0;
  totalSquaresWithPairs = 0;
  totalProcessedSquaresSmall = 0;
  totalProcessedSquaresLarge = 0;

  startParametersFileName = "start_parameters.txt";
  resultFileName = "result.txt";
  checkpointFileName = "checkpoint.txt";
  tempCheckpointFileName = "tmp_checkpoint.txt";

  // Задание константы - заголовка в файле параметров или контрольной точке
  moveSearchGlobalHeader = "# Move search of pairs OLDS status";
  moveSearchComponentHeader = "# Move search component status";

  // Сброс флага инициализации
  isInitialized = 0;
}


// Очистка необходимых переменных перед очередным поиском
void MovePairSearch::ClearBeforeNextSearch()
{
  // Сброс значений матриц квадратов A и B, а также матрицы использования строк 
  // при формировании квадрата B
  for (int i = 0; i < Rank; i++)
  {
    for (int j = 0; j < Rank; j++)
    {
      squareA[i][j] = -1;
      squareB[i][j] = -1;
      rowsHistory[i][j] = 1;
    }
  }

  // Сброс значений в вектрах использования строк в очередной перестановке и номеров строк, использованных для текущего квадрата
  for (int i = 0; i < Rank; i++)
  {
    rowsUsage[i] = 1;
    currentSquareRows[i] = -1;
  }

  // Сброс значения счётчиков пар найденных для заданного ДЛК
  pairsCount = 0;
}


// Инициализация поиска
void MovePairSearch::InitializeMoveSearch(string start, string result, 
                                       string checkpoint, string temp)
{
  fstream startFile;
  fstream checkpointFile;

  // Считывание названий имен файлов
  startParametersFileName = start;
  resultFileName = result;
  checkpointFileName = checkpoint;
  tempCheckpointFileName = temp;

  // Считываем состояние генератора и поиска из файла контрольной точки или начальных значений
    // Открытие файлов со стартовыми параметрами и файла контрольной точки
    startFile.open(startParametersFileName.c_str(), std::ios_base::in);
    checkpointFile.open(checkpointFileName.c_str(), std::ios_base::in);

    // Считываение состояния
	if (checkpointFile.is_open())
	{
		// Считывание состояния из файла контрольной точки
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
		// Считывание состояния из файла стартовых параметров
		Read(startFile);
		isStartFromCheckpoint = 0;
    }

    // Закрытие файлов
    startFile.close();
    checkpointFile.close();
}


// Чтение состояния поиска из потока
void MovePairSearch::Read(istream& is)
{
  string marker;

  // Сброс флага инициализированности
  isInitialized = 0;

  // Считывание состояния поиска
    // Находим маркер начала состояния
    do
	{
		std::getline(is, marker);

		if (is.eof())
		{
			throw ("Expected start marker, but EOF found.");
		}
    }
    while (marker != moveSearchGlobalHeader);
    
    // Считываем состояние генератора ДЛК
    is >> squareAGenerator;

    // Находим маркер компоненты перетасовки
	do
	{
		std::getline(is, marker);

		if (is.eof())
		{
			throw ("Expected start marker, but EOF found.");
		}
	}
    while (marker != moveSearchComponentHeader);

    // Считываем переменные поиска перетасовкой (по факту - переменные со статистикой)
    is >> pairsCount;
    is >> totalPairsCount;
    is >> totalSquaresWithPairs;
    is >> totalProcessedSquaresLarge;
    is >> totalProcessedSquaresSmall;

  // Выставление флага инициализированности
  isInitialized = 1;
}


// Запись состояния поиска в поток
void MovePairSearch::Write(ostream& os)
{
  // Запись состояния поиска
    // Запись заголовка
    os << moveSearchGlobalHeader << endl;
    os << endl;

    // Запись состояния генератора ДЛК
    os << squareAGenerator;
    os << endl;

    // Запись заголовка блока перетасовки
    os << moveSearchComponentHeader << endl;
    os << endl;

    // Запись статистических показателей
    os << pairsCount << " " << totalPairsCount << " " << totalSquaresWithPairs << endl;
    os << totalProcessedSquaresLarge << " " << totalProcessedSquaresSmall << endl;
    os << endl;
}


// Создание контрольной точки
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


// Запуск поиска ортогональных квадратов методом перестановки строк
void MovePairSearch::StartMoveSearch()
{
  // Подписываемся на событие нахождения очередного ДЛК
  squareAGenerator.Subscribe(this); 

  // Запускаем генерацию ДЛК
  squareAGenerator.Start();

  // Отписываемся от события нахождения очередного ДЛК
  squareAGenerator.Unsubscribe();

  // Вывод итогов поиска
  ShowSearchTotals();
}


// Обработчик события построения ДЛК, запускающий поиск к нему пары
void MovePairSearch::OnSquareGenerated(Square newSquare)
{
  // Очистка перед поиском квадратов, ортогональных очередному ДЛК
  ClearBeforeNextSearch();

  // Запоминание найденного квадрата
  for (int i = 0; i < Rank; i++)
  {
    for (int j = 0; j < Rank; j++)
    {
      squareA[i][j] = newSquare.Matrix[i][j];
    }
  }

  // Запуск перетасовки строк
  MoveRows();

  // Проверка взаимной ортогональности квадратов
  if (pairsCount > 0)
  {
    CheckMutualOrthogonality();
  }

  // Собирание статистики по обработанным квадратам
  totalProcessedSquaresSmall++;

  // Фиксация информации о ходе обработки
  if (totalProcessedSquaresSmall % CheckpointInterval == 0)
  {
    // Обновить прогресс выполнения для клиента BOINC
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

    boinc_fraction_done(fraction_done); // Сообщить клиенту BOINC о доле выполнения задания

    // Проверка, может ли клиент BOINC создать контрольную точку,
    // и если может, то запустить функцию её записи
    if (boinc_time_to_checkpoint()) {
      CreateCheckpoint();
      boinc_checkpoint_completed(); // BOINC знает, что контрольная точка записана
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

// Перетасовка строк заданного ДЛК в поиске ОДЛК к нему
void MovePairSearch::MoveRows()
{
  int currentRowId = 1;
  int isRowGet = 0;
  int gettingRowId = -1;
  int oldRowId = -1;

  int diagonalValues[Rank];
  int duplicationDetected = 0;

  // Записываем первую строку квадрата A в B с целью поиска нормализованных квадратов
  for (int j = 0; j < Rank; j++)
  {
    squareB[0][j] = squareA[0][j];
  }

  // Отмечаем задействование первой строки, т.к. она фиксированная
  rowsUsage[0] = 0;
  rowsHistory[0][0] = 0;
  currentSquareRows[0] = 0;

  while (currentRowId > 0)
  {
    // Подбираем строку из исходного квадрата на позицию currentRowId формируемого квадрата
    isRowGet = 0;
    gettingRowId = -1;
    for (int i = 0; i < Rank; i++)
    {
      // Проверяем i-ю строку исходного квадрата
      if (rowsUsage[i] && rowsHistory[currentRowId][i])
      {
        isRowGet = 1;
        gettingRowId = i;

        break;
      }
    }

    // Обрабатываем результат поиска
    if (isRowGet)
    {
      // Обрабатываем нахождение новой строки
        // Заносим в квадрат новую строку
          // Считываем номер строки, которая сейчас стоит в квадрате
          oldRowId = currentSquareRows[currentRowId];
          // Записываем новую строку в квадрат, массив флагов использованных строк, 
          //в историю использованных строк и массив текущих строк квадрата
            // Записываем новую строку в квадрат
            for (int j = 0; j < Rank; j++)
            {
              squareB[currentRowId][j] = squareA[gettingRowId][j];
            }
            // Отмечаем строку в массие используемых строк
            rowsUsage[gettingRowId] = 0;
            // Отмечаем строку в истории использования строки
            rowsHistory[currentRowId][gettingRowId] = 0;
            // Записываем строку в массив текущих строк квадрата
            currentSquareRows[currentRowId] = gettingRowId;

        // Очищаем для предыдущей строки флаги использования
          // Убираем отметку в массиве используемых строк
          if (oldRowId != -1)
          {
            rowsUsage[oldRowId] = 1;
          }

        // Проверяем диагональность получающейся части квадрата
          // Сбрасываем флаг сигнализирующий и дубликатах на диагоналях
          duplicationDetected = 0;
          // Проверка главной диагонали
            // Сбрасываем флаги использованных значений
            for (int i = 0; i < Rank; i++)
            {
              diagonalValues[i] = 1;
            }
            // Проверка значений главной диагонали
            for (int i = 0; i <= currentRowId; i++)
            {
              // Проверка i-го элемента главной диагонали - клетки (i, i)
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
          // Проверка побочной диагонали, если это имеет смысл
          if (!duplicationDetected)
          {
            // Проверка побочной диагонали
              // Сбрасываем флаги использованных значений
              for (int i = 0; i < Rank; i++)
              {
                diagonalValues[i] = 1;
              }
              // Проверка значений побочной диагонали начиная с "её хвоста"
              for (int i = 0; i <= currentRowId; i++)
              {

                // Проверка i-го значения побочной диагонали - элемента (i, rank - 1 - i)
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

        // Обработка итогов проверки диагональности квадрата
        if (!duplicationDetected)
        {
          // Делаем следующий шаг вперёд в зависимости от текущего положения
          if (currentRowId == Rank - 1)
          {
            // Обрабатываем найденный квадрат
            ProcessOrthoSquare();
          }
          else
          {
            // Делаем шаг вперёд
            currentRowId++;
          }
        }
    }
    else
    {
      // Обрабатываем ненахождение новой строки - делаем шаг назад, зачищая флаги задействования, 
      // историю использования, перечень текущих строк квадрата и зачищая сам квадрат
        // Считываем номер текущей строки
        oldRowId = currentSquareRows[currentRowId];
        // Зачищаем текущую строку в квадрате
        for (int j = 0; j < Rank; j++)
        {
          squareB[currentRowId][j] = -1;
        }
        // Зачищаем текущий состав квадрата
        currentSquareRows[currentRowId] = -1;
        // Зачищаем флаг возможного задействования
        rowsUsage[oldRowId] = 1;
        // Зачищаем историю работы с этой строкой
        for (int i = 0; i < Rank; i++)
        {
          rowsHistory[currentRowId][i] = 1;
        }
        // Делаем шаг назад
        currentRowId--;
    }
  }
}


// Обработка найденного ортогонального квадрата
void MovePairSearch::ProcessOrthoSquare()
{
  int isDifferent = 0;      // Число отличий в строках от исходного квадрата (для отсева формирования его копии)
  ofstream resultFile;        // Поток для I/O в файл с результатами

  Square a(squareA);        // Квадрат A как объект
  Square b(squareB);        // Квадрат B как объект

  int orthoMetric = Rank*Rank;  // Значение метрики ортогональности, говорящее о том, что квадраты - полностью ортогональны

  // Проверяем его на то, что он копия исходного
  isDifferent = 0;
  
  for (int i = 0; i < Rank; i++)
  {
    if (currentSquareRows[i] != i)
    {
      isDifferent = 1;
      break;
    }
  }

  // Обработка найденного квадрата
  if (isDifferent && Square::OrthoDegree(a, b) == orthoMetric 
      && b.IsDiagonal() && b.IsLatin() && a.IsDiagonal() && b.IsLatin())
  {
    // Запись информации о найденном квадрате
      // Увеличение счётчика квадратов
      pairsCount++;
      totalPairsCount++;

      // Запоминание базового квадрата
      if (pairsCount == 1)
      {
        orthoSquares[pairsCount - 1] = a;
        totalSquaresWithPairs++;
      }

      // Запоминание квадрата - пары
      if (pairsCount < OrhoSquaresCacheSize)
      {
        orthoSquares[pairsCount] = b;
      }

      // Вывод заголовка
      if (pairsCount == 1)
      {
        if(isDebug)
        {
          // Вывод информации о первом квадрате пары в виде заголовка
          cout << "{" << endl;
          cout << "# ------------------------" << endl;
          cout << "# Detected pair for the square: " << endl;
          cout << "# ------------------------" << endl;
          cout << a;
          cout << "# ------------------------" << endl;
        }
        // Вывод информации в файл
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

      // Вывод информации о найденной паре
        if(isDebug)
        {
          // Вывод информации в консоль
          cout << b << endl;
        }

        // Вывод информации в файл
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


// Проверка взаимной ортогональности набора квадратов, найденного в текущем поиске
void MovePairSearch::CheckMutualOrthogonality()
{
  int orthoMetric = Rank*Rank;
  int maxSquareId;
  ofstream resultFile;

  // Определение верхней границы обрабатываемых квадратов
  if (pairsCount < OrhoSquaresCacheSize)
  {
    maxSquareId = pairsCount;
  }
  else
  {
    maxSquareId = OrhoSquaresCacheSize - 1;
  }

  // Открываем файл с результатами
  resultFile.open(resultFileName.c_str(), std::ios_base::binary | std::ios_base::app);
  if (!resultFile.is_open()) { cout << "Error opening file!"; return; }

  // Проверка взаимной ортогональности набора квадратов
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

  // Выводим общее число найденых ОДЛК
  if(isDebug) cout << "# Pairs found: " << pairsCount << endl;
  resultFile << "# Pairs found: " << pairsCount << endl;

  // Ставим отметку об окончании секции результатов
  if(isDebug) cout << "}" << endl;
  resultFile << "}" << endl;

  // Закрываем файл с результатами
  resultFile.close();
}


void MovePairSearch::ShowSearchTotals()
{
  ofstream resultFile;

  if(isDebug)
  {
    // Вывод итогов в консоль
    cout << "# ------------------------" << endl;
    cout << "# Total pairs found: " << totalPairsCount << endl;
    cout << "# Total squares with pairs: " << totalSquaresWithPairs << endl;
    cout << "# ------------------------" << endl;
  }

  // Вывод итогов в файл
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
