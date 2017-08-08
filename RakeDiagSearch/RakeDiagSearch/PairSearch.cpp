// Поиск пар ортогональных диагональных латинских квадратов

# include "PairSearch.h"
# include <string>

using namespace std;

// Конструктор по умолчанию
PairSearch::PairSearch()
{
	// Сброс всех настроек
	Reset();

	// Задание текстовых констант
	pairStateHeader = "# Search of pairs ODLS status";
}

// Сброс всех значений внутренних структур
void PairSearch::Reset()
{
	// Сброс генератора первого квадрата и генератора квадрата - пары
	squareAGenerator.Reset();
	Generator::Reset();

	// Очистка массива пар значений
	for (int i = 0; i < Rank; i++)
	{
		for (int j = 0; j < Rank; j++)
		{
			pairsDictionary[i][j] = Free;
		}
	}

	// Сброс числа найденных пар и общего числа квадратов с парами
	pairsCount = 0;
	totalPairsCount = 0;
	totalSquaresWithPairs = 0;
	
	// Сброс флага инициализированности
	isInitialized = No;
}


// Инициализация поиска
void PairSearch::Initialize(string start, string result, string checkpoint, string temp)
{
	fstream startFile;
	fstream checkpointFile;
	string marker;

	// Сбрасываем флаг инициализированности в "Нет"
	isInitialized = No;

	// Сбрасываем значения внутренних структур
	Reset();

	// Запоминаем имена файлов конфигурации, контрольных точек и результатов
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
			Read(checkpointFile);
			isStartFromCheckpoint = Yes;
		}
		else
		{
			// Считывание состояния из файла стартовых параметров
			Read(startFile);
			isStartFromCheckpoint = No;
		}

		// Закрытие файлов
		startFile.close();
		checkpointFile.close();

		// Считывание стартовых параметров для перезапуска поиска ортогональных квадратов к очередному квадрату
		startFile.open(startParametersFileName.c_str(), std::ios_base::in);

		if (startFile.is_open())
		{
			// Стартовых параметров
				// Поиск заголовка
				do
				{
					std::getline(startFile, marker);
				}
				while (marker != pairStateHeader);

				// Считывание параметров
				startFile >> initialState;
				initialState.SetFileNames(start, result, checkpoint, temp);
		}

		startFile.close();
}


// Считывание состояние поиска из потока
void PairSearch::Read(std::istream& is)
{
	int result = Yes;

	// Считываем состояние генератора диагональных латинских квадратов
	is >> squareAGenerator;

	// Считываем состояние базового класса поиска
	Generator::Read(is);

	// Считываем дополнительные переменные класса поиска ортогональных квадратов
	if (isInitialized)
	{
		is >> pairsCount;
		is >> totalPairsCount;
		is >> totalSquaresWithPairs;
		
		for (int i = 0; i < Rank; i++)
		{
			for (int j = 0; j < Rank; j++)
			{
				is >> pairsDictionary[i][j];
			}
		}
	}

	// Считываем квадрат, к которому ищем ортогональные квадраты
	is >> squareA;
}


// Запись состояния поиска в поток
void PairSearch::Write(std::ostream& os)
{
	// Запись состояния генератора
	os << squareAGenerator << endl;

	// Запись строки с заголовком
	os << pairStateHeader << endl << endl;

	// Запись состояния базового класса поиска
	Generator::Write(os);
	os << endl;

	// Запись дополнительных переменных, связанных с поиском ортогональных квадратов
	os << pairsCount << " " << totalPairsCount << " " << totalSquaresWithPairs << endl;
	os << endl;

	for (int i = 0; i < Rank; i++)
	{
		for (int j = 0; j < Rank; j++)
		{
			os << pairsDictionary[i][j] << " ";
		}
		os << endl;
	}
	os << endl;

	// Записываем квадрат к которому ищем ортогональные квадраты
	os << squareA << endl;
}


// Оператор записи состояня поиска
std::ostream& operator << (std::ostream& os, PairSearch& value)
{
	value.Write(os);

return os;
}


// Оператор считывания состояния поиска
std::istream& operator >> (std::istream& is, PairSearch& value)
{
	value.Read(is);

return is;
}


// Сброс состояния поиска на необходимое для обработки очередного квадрата
void PairSearch::ResetForProcessSquare()
{
	// Восстанавливаем состояние начала поиска пару к найденному ДЛК
	Generator::CopyState(initialState);

	// Очистка массива пар значений
	for (int i = 0; i < Rank; i++)
	{
		for (int j = 0; j < Rank; j++)
		{
			pairsDictionary[i][j] = Free;
		}
	}

	// Сброс числа найденных пар
	pairsCount = 0;

	// Выставление флага инициализированности
	isInitialized = Yes;
}

// Обработчик сгенерированного диагонального латинского квадрата
void PairSearch::OnSquareGenerated(Square newSquare)
{
	// Сбрасываем состояние поиска на необходимое для обработки очередного квадрата
	ResetForProcessSquare();

	// Считываем матрицу сгенерированного квадрата к котрому ищем пару
	squareA = newSquare;

	// Запускаем поиск квадратов, ортогональных заданному
	FindPairSquare();

	// Вывод результатов поиска, если таковые были
	if (pairsCount > 0)
	{
		PrintSearchFooter();
	}
}

// Запуск поиска пар ДЛК
void PairSearch::StartPairSearch()
{
	// Подписываемся на событие генерации нового ДЛК
        squareAGenerator.SquareGenerated.connect(boost::bind(
                   &PairSearch::OnSquareGenerated, this, _1));

	// Заканчиваем обработку текущего квадрата
	if (isStartFromCheckpoint == Yes)
	{
		FindPairSquare();
	}

	// Запускаем генерацию диагональных латинских квадратов
	squareAGenerator.Start();

	// Отписываемся от события генерации нового ДЛК
        squareAGenerator.SquareGenerated.disconnect_all_slots();

	// Вывод итогов поиска пар
	PrintSearchTotals();
}


// Обработка нахождения пары к сгенерированному ДЛК
void PairSearch::ProcessPairSquare()
{
	ofstream resultFile;		// Поток для I/O в файл с результатами

	// Учитываем найденную пару
	pairsCount++;
	totalPairsCount++;

	if (pairsCount == 1)
	{
		totalSquaresWithPairs++;
	}

	// Вывод информации о найденной паре
		// Вывод заголовка
		if (pairsCount == 1)
		{
			// Вывод информации о первом квадрате пары в виде заголовка
			cout << "# ------------------------" << endl;
			cout << "# Detected orthogonal square(s) for the square: " << endl;
			cout << endl;
			cout << squareA;
			cout << "# ------------------------" << endl;

			// Вывод информации в файл
			resultFile.open(resultFileName.c_str(), std::ios_base::app);
			resultFile << "# ------------------------" << endl;
			resultFile << "# Detected orthogonal square(s) for the square: " << endl;
			resultFile << endl;
			resultFile << squareA;
			resultFile << "# ------------------------" << endl;
			resultFile.close();
		}

		// Вывод информации в консоль
		cout << newSquare;

		// Вывод информации в файл
		resultFile.open(resultFileName.c_str(), std::ios_base::app);
		resultFile << newSquare;
		resultFile.close();

	// Создание контрольной точки
	CreateCheckpoint();
}


// Вывод информации об итогах поиска квадратов, ортогональных к заданному ДЛК
void PairSearch::PrintSearchFooter()
{
	ofstream resultFile;

	// Вывод итогов в консоль
	cout << "# ------------------------" << endl;
	cout << "# Pairs found: " << pairsCount << endl;
	cout << "# ------------------------" << endl;

	// Вывод итогов в файл
	resultFile.open(resultFileName.c_str(), std::ios_base::app);
	resultFile << "# ------------------------" << endl;
	resultFile << "# Pairs found: " << pairsCount << endl;
	resultFile << "# ------------------------" << endl;
	resultFile.close();
}


// Вывод информации об итогах всего поиска в целом
void PairSearch::PrintSearchTotals()
{
	ofstream resultFile;

	// Вывод итогов в консоль
	cout << "# ------------------------" << endl;
	cout << "# Total pairs found: " << totalPairsCount << endl;
	cout << "# Total squares with pairs: " << totalSquaresWithPairs << endl;
	cout << "# ------------------------" << endl;

	// Вывод итогов в файл
	resultFile.open(resultFileName.c_str(), std::ios_base::app);
	resultFile << "# ------------------------" << endl;
	resultFile << "# Total pairs found: " << totalPairsCount << endl;
	resultFile << "# Total squares with pairs: " << totalSquaresWithPairs << endl;
	resultFile << "# ------------------------" << endl;
	resultFile.close();
}


// Создание контрольной точки
void PairSearch::CreateCheckpoint()
{
	fstream checkpointFile;

	checkpointFile.open(tempCheckpointFileName.c_str(), std::ios_base::out);
	if (checkpointFile.is_open())
	{
		Write(checkpointFile);
		checkpointFile.close();
		remove(checkpointFileName.c_str());
		rename(tempCheckpointFileName.c_str(), checkpointFileName.c_str());
	}
}


// Поиск квадрата, парного к найденному
void PairSearch::FindPairSquare()
{
	int stepCounter = 0;

	int isGet;			// Флаг получения нового значения для клетки
	int cellValue;		// Новое значение для клетки
	int oldCellValue;	// Старое значение, стоявшее в клетке

	int stop = 0;		// Флаг достижения окончания расчёта

	if (isInitialized == Yes)
	{
		// Подбор значений клеток квадрата
		do
		{
			/* ---------- */
			if (stepCounter % 100000000 == 0)
			{
				stepCounter = 0;
				cout << "Current state: " << endl;
				cout << newSquare << endl;
			}
			stepCounter++;
			/* ---------- */

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

								// Проверяем получающуюся пару на её задействование другими ячейками
								// первого квадрата - squareA и проверяемого значения для второго - i
								if (pairsDictionary[squareA.Matrix[rowId][columnId]][i] == 0)
								{
									isGet = 0;
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

						// Отмечаем значение найденной пары в использованных
							// Отмечаем использованние новой пары
							pairsDictionary[squareA.Matrix[rowId][columnId]][cellValue] = 0;

							// Снимаем отметку с использования предыдущей пары
							if (oldCellValue != -1)
							{
								pairsDictionary[squareA.Matrix[rowId][columnId]][oldCellValue] = 1;
							}

						// Обработка окончания формирования квадрата
						if (cellId == cellsInPath - 1)
						{
							// Обрабатываем найденный квадрат
							ProcessPairSquare();
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

								// Отмечаем освобождения значения в использованных парах
								pairsDictionary[squareA.Matrix[rowId][columnId]][cellValue] = 1;
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
