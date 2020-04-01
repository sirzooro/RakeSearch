// Поиск пар диагональных латинских квадратов методом "перетасовки" строк

#include "RakeSearch.h"

// Конструктор по умолчанию
RakeSearch::RakeSearch()
{
    Reset();
}

// Задание имен файлов параметров и контрольной точки
void RakeSearch::SetFileNames(const string &start, const string &result, const string &checkpoint, const string &temp)
{
    startParametersFileName = start;
    resultFileName = result;
    checkpointFileName = checkpoint;
    tempCheckpointFileName = temp;
}

// Сброс значений внутренних структур
void RakeSearch::Reset()
{
    // Очистка матриц квадратов
    for (int i = 0; i < Rank; i++)
    {
        for (int j = 0; j < Rank; j++)
        {
            squareA[i][j] = Square::Empty;
            squareB[i][j] = Square::Empty;
        }
    }

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
    flagsPrimary = (1u << Rank) - 1;
    flagsSecondary = (1u << Rank) - 1;

    // Сброс значений в матрицах использования элементов в столбцах и строках
    for (int i = 0; i < Rank; i++)
    {
        flagsColumns[i] = (1u << Rank) - 1;
        flagsRows[i] = (1u << Rank) - 1;
    }

    // Сброс значений в кубе истории использования значений в клетках
    for (int i = 0; i < Rank; i++)
    {
        for (int j = 0; j < Rank; j++)
        {
            flagsCellsHistory[i][j] = (1u << Rank) - 1;
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

    // Сброс значений, связанных с перестановкой строк
    // Сброс числа найденых пар для заданного ДЛК
    pairsCount = 0;

    // Сброс значений глобальных счётчиков
    totalPairsCount = 0;
    totalSquaresWithPairs = 0;

    // Задание имён входных файлов
    startParametersFileName = "start_parameters.txt";
    resultFileName = "result.txt";
    checkpointFileName = "checkpoint.txt";
    tempCheckpointFileName = "tmp_checkpoint.txt";

    // Задание константы - заголовка в файле параметров или контрольной точке
    workunitHeader = "# RakeSearch of diagonal Latin squares";

    // Сброс флага инициализации
    isInitialized = 0;
}

// Инициализация поиска
void RakeSearch::Initialize(const string &start, const string &result, const string &checkpoint, const string &temp)
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

    // Считывание состояния из файла контрольной точки
    if (checkpointFile.is_open())
    {
        // Считывание состояния из существующего файла контрольной точки
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

    // Считывание состояния из файла стартовых параметров
    if (isStartFromCheckpoint != 1)
    {
        // Считывание состояния из существующего файла стартовых параметров
        Read(startFile);
        isStartFromCheckpoint = 0;
    }

    // Закрытие файлов
    startFile.close();
    checkpointFile.close();
}

// Чтение состояния поиска из потока
void RakeSearch::Read(istream &is)
{
    string marker;
    int rankToVerify;
    unsigned int storedBit = 0;
    Square currentSquare;

    // Сброс флага инициализированности
    isInitialized = 0;

    // Считывание состояния поиска
    // Находим маркер начала состояния
    do
    {
        std::getline(is, marker);

        if (is.eof())
        {
            throw("Expected start marker, but EOF found.");
        }
    } while (marker != workunitHeader);

    // Считываем состояние генератора ДЛК
    // Считывание из потока ранга квадрата
    is >> rankToVerify;

    // Считывание данных поиска нужного нам ранга
    if (rankToVerify == Square::Rank)
    {
        // Считывание из потока квадрата A - первого квадрата пары
        is >> currentSquare;
        for (int i = 0; i < Rank; i++)
        {
            for (int j = 0; j < Rank; j++)
            {
                squareA[i][j] = currentSquare.Matrix[i][j];
            }
        }

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
            is >> storedBit;
            if (storedBit)
            {
                flagsPrimary |= 1u << i;
            }
            else
            {
                flagsPrimary &= ~(1u << i);
            }
        }

        // Считывание информации о значениях на побочной диагонали
        for (int i = 0; i < Rank; i++)
        {
            is >> storedBit;
            if (storedBit)
            {
                flagsSecondary |= 1u << i;
            }
            else
            {
                flagsSecondary &= ~(1u << i);
            }
        }

        // Считывание информации о значениях в строках
        for (int i = 0; i < Rank; i++)
        {
            for (int j = 0; j < Rank; j++)
            {
                is >> storedBit;
                if (storedBit)
                {
                    flagsRows[i] |= 1u << j;
                }
                else
                {
                    flagsRows[i] &= ~(1u << j);
                }
            }
        }

        // Считывание информации о значениях в столбцах
        for (int i = 0; i < Rank; i++)
        {
            for (int j = 0; j < Rank; j++)
            {
                is >> storedBit;
                if (storedBit)
                {
                    flagsColumns[i] |= 1u << j;
                }
                else
                {
                    flagsColumns[i] &= ~(1u << j);
                }
            }
        }

        // Считывание информации об истории значений в клетках квадрата
        for (int h = 0; h < Rank; h++)
        {
            for (int i = 0; i < Rank; i++)
            {
                for (int j = 0; j < Rank; j++)
                {
                    is >> storedBit;
                    if (storedBit)
                    {
                        flagsCellsHistory[i][j] |= 1u << h;
                    }
                    else
                    {
                        flagsCellsHistory[i][j] &= ~(1u << h);
                    }
                }
            }
        }

        // Считываем число сгенерированных квадратов
        is >> squaresCount;

        // Выставляем флаг инициализированности
        isInitialized = Yes;
    }

    // Считываем переменные поиска перетасовкой (по факту - переменные со статистикой)
    is >> pairsCount;
    is >> totalPairsCount;
    is >> totalSquaresWithPairs;

    // Выставление флага инициализированности
    isInitialized = 1;
}

// Запись состояния поиска в поток
void RakeSearch::Write(std::ostream &os)
{
    Square currentSquare(squareA); // Первый квадрат пары, сформированный к моменту записи

    // Запись состояния поиска
    // Запись заголовка
    os << workunitHeader << endl;
    os << endl;

    // Запись состояния генератора ДЛК
    // Запись в поток ранга квадрата
    os << Square::Rank << endl;

    // Запись в поток квадрата
    os << currentSquare;

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
    os << rowId << " " << columnId << " " << cellId << endl;

    // Записываем пустую строку для удобства
    os << endl;

    // Запись информации о задействованных значениях и истории значений
    // Запись информации о значениях на главной диагонали
    for (int i = 0; i < Rank; i++)
    {
        os << ((flagsPrimary & (1u << i)) > 0 ? 1 : 0) << " ";
    }
    os << endl;

    // Запись информации о значениях на побочной диагонали
    for (int i = 0; i < Rank; i++)
    {
        os << ((flagsSecondary & (1u << i)) > 0 ? 1 : 0) << " ";
    }
    os << endl;

    // Дополнительная пустая строка
    os << endl;

    // Запись информации о значениях в строках
    for (int i = 0; i < Rank; i++)
    {
        for (int j = 0; j < Rank; j++)
        {
            os << ((flagsRows[i] & (1u << j)) > 0 ? 1 : 0) << " ";
        }
        os << endl;
    }
    os << endl;

    // Запись информации о значениях в столбцах
    for (int i = 0; i < Rank; i++)
    {
        for (int j = 0; j < Rank; j++)
        {
            os << ((flagsColumns[i] & (1u << j)) > 0 ? 1 : 0) << " ";
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
                os << ((flagsCellsHistory[i][j] & (1u << h)) > 0 ? 1 : 0) << " ";
            }
            os << endl;
        }
        os << endl;
    }
    os << endl;

    // Запись в поток информации о числе сгенерированных квадратов
    os << squaresCount << endl;
    os << endl;

    // Запись статистических показателей
    os << pairsCount << " " << totalPairsCount << " " << totalSquaresWithPairs << endl;
    os << endl;
}

// Создание контрольной точки
void RakeSearch::CreateCheckpoint()
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
    {
        cerr << "Error opening checkpoint file!" << endl;
    }
}

// Обработка найденного, возможно что ортогонального квадрата
void RakeSearch::ProcessOrthoSquare()
{
    ofstream resultFile; // Поток для I/O в файл с результатами

    Square a(squareA); // Квадрат A как объект
    Square b(squareB); // Квадрат B как объект

    int orthoDegree = -1; // Метрика ортогональности проверяемых квадратов

    // Обработка найденного квадрата
    orthoDegree = Square::OrthoDegree(a, b);
    if (orthoDegree >= MinOrthoMetric && b.IsDiagonal() && b.IsLatin() && a.IsDiagonal() && a.IsLatin())
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
            if (isDebug)
            {
                // Вывод информации о первом квадрате пары в виде заголовка
                cout << "{" << endl;
                cout << "# ------------------------" << endl;
                cout << "# Detected pair for the square: " << endl;
                cout << "# Degree of orthogonality: " << orthoDegree << endl;
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
                resultFile << "# Degree of orthogonality: " << orthoDegree << endl;
                resultFile << "# ------------------------" << endl;
                resultFile << a;
                resultFile << "# ------------------------" << endl;
                resultFile.close();
            }
            else
            {
                std::cerr << "Error opening file!";
            }
            // Создание контрольной точки
            CreateCheckpoint();
            boinc_checkpoint_completed();
        }

        // Вывод информации о найденной паре
        if (isDebug)
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
void RakeSearch::CheckMutualOrthogonality()
{
    int orthoMetric = Rank * Rank;
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
    if (!resultFile.is_open())
    {
        cout << "Error opening file!";
        return;
    }

    // Проверка взаимной ортогональности набора квадратов
    for (int i = 0; i <= maxSquareId; i++)
    {
        for (int j = i + 1; j <= maxSquareId; j++)
        {
            if (Square::OrthoDegree(orthoSquares[i], orthoSquares[j]) == orthoMetric)
            {
                if (isDebug)
                    cout << "# Square " << i << " # " << j << endl;
                resultFile << "# Square " << i << " # " << j << endl;
            }
        }
    }
    if (isDebug)
        cout << endl;
    resultFile << endl;

    // Выводим общее число найденых ОДЛК
    if (isDebug)
        cout << "# Pairs found: " << pairsCount << endl;
    resultFile << "# Pairs found: " << pairsCount << endl;

    // Ставим отметку об окончании секции результатов
    if (isDebug)
        cout << "}" << endl;
    resultFile << "}" << endl;

    // Закрываем файл с результатами
    resultFile.close();
}

// Обработка квадрата
void RakeSearch::ProcessSquare()
{
    double fraction_done; // Доля выполнения задания

    // Увеличиваем счётчик найденных квадратов
    squaresCount++;

    // Запуск перетасовки строк
    PermuteRows();

    // Проверка взаимной ортогональности квадратов
    if (pairsCount > 0)
    {
        CheckMutualOrthogonality();
    }

    // Фиксация информации о ходе обработки
    if (squaresCount % CheckpointInterval == 0)
    {
        // Обновить прогресс выполнения для клиента BOINC
        fraction_done = squaresCount / 1000000000;

        if (fraction_done >= 1)
            fraction_done = 0.999999999;

        boinc_fraction_done(fraction_done); // Сообщить клиенту BOINC о доле выполнения задания

        // Проверка, может ли клиент BOINC создать контрольную точку,
        // и если может, то запустить функцию её записи
        if (boinc_time_to_checkpoint())
        {
            CreateCheckpoint();
            boinc_checkpoint_completed(); // BOINC знает, что контрольная точка записана
        }

        if (isDebug)
        {
            Square squareToShow(squareA);

            cout << "# ------------------------" << endl;
            cout << "# Processed " << squaresCount << " squares." << endl;
            cout << "# Last processed square:" << endl;
            cout << endl;
            cout << squareToShow;
            cout << "# ------------------------" << endl;
        }
    }
}

// Вывод итогов поиска
void RakeSearch::ShowSearchTotals()
{
    ofstream resultFile;

    if (isDebug)
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
        resultFile << "# Processed " << squaresCount << " squares" << endl;
        resultFile << "# ------------------------" << endl;
        resultFile.close();
    }
    else
        cerr << "Error opening file!" << endl;
}

// Запуск генерации квадратов
void RakeSearch::Start()
{
#ifdef __arm__
    unsigned int reversedMask; // Перевёрнутая битовая маска с флагами (для ARM, команды clz!)
#endif
    unsigned long bitIndex; // Позиция выставленного бита в "массиве" флагов
    unsigned int freeValuesMask; // Итоговая маска из битов - флагов занятости значений
    unsigned int isGet; // Флаг получения нового значения для клетки
    int cellValue;      // Новое значение для клетки
    int oldCellValue;   // Старое значение, стоявшее в клетке

    int stop = 0; // Флаг достижения окончания расчёта

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
            freeValuesMask = (1u << Rank) - 1;

            // Подбираем значение для клетки
            // Так как для 9 и 10 ранга в рамках маски workunit-а уже заданы обе диагонали, то проверку по ним - отключаем
            // // Применяем маску главной диагонали
            // if (rowId == columnId)
            // {
            // 	freeValuesMask &= flagsPrimary;
            // }
            //
            // // Применяем маску побочной диагонали
            // if (Rank - rowId - 1 == columnId)
            // {
            // 	freeValuesMask &= flagsSecondary;
            // }

            // Применяем маску строк, столбцов и истории значений
            freeValuesMask &= flagsColumns[columnId] & flagsRows[rowId] & flagsCellsHistory[rowId][columnId];

// Определяем минимально возможное для использования в клетке значение
#ifdef __arm__
            asm("rbit %1, %0" : "=r"(reversedMask) : "r"(freeValuesMask) :);
            asm("clz %1, %0" : "=r"(bitIndex) : "r"(reversedMask) :);
            if (bitIndex < Rank)
            {
                cellValue = bitIndex;
                isGet = 1;
            }
            else
            {
                isGet = 0;
            }
#else
#ifdef _MSC_VER
            isGet = _BitScanForward(&bitIndex, freeValuesMask);
            cellValue = bitIndex;
#else
            bitIndex = __builtin_ffs(freeValuesMask);
            if (bitIndex > 0)
            {
                isGet = 1;
                cellValue = bitIndex - 1;
            }
            else
            {
                isGet = 0;
            }
#endif
#endif

            // Обработка результата поиска
            if (isGet)
            {
                // Обработка найденного нового значения
                // Считывание текущего значения
                oldCellValue = squareA[rowId][columnId];
                // Запись нового значения
                // Записываем значение в квадрат
                squareA[rowId][columnId] = cellValue;
                // Отмечаем значение в столбцах
                flagsColumns[columnId] ^= 1u << cellValue;
                // Отмечаем значение в строках
                flagsRows[rowId] ^= 1u << cellValue;
                // Так как для 9 и 10 ранга в рамках маски workunit-а уже заданы обе диагонали, то проверку по ним - отключаем
                // // Отмечаем значение в диагоналях
                // if (rowId == columnId)
                // {
                // 	flagsPrimary ^= 1u << cellValue;
                // }
                // if (rowId == Rank - 1 - columnId)
                // {
                // 	flagsSecondary ^= 1u << cellValue;
                // }
                // Отмечаем значение в истории значений клетки
                flagsCellsHistory[rowId][columnId] ^= 1u << cellValue;

                // Возвращение предыдущего значения без зачистки истории (так как мы работаем с этой клеткой)
                if (oldCellValue != Square::Empty)
                {
                    // Возвращаем значение в столбцы
                    flagsColumns[columnId] |= 1u << oldCellValue;
                    // Возвращаем значение в строки
                    flagsRows[rowId] |= 1u << oldCellValue;
                    // Так как для 9 и 10 ранга в рамках маски workunit-а уже заданы обе диагонали, то проверку по ним - отключаем
                    // // Возвращаем значение в диагонали
                    // if (rowId == columnId)
                    // {
                    // 	flagsPrimary |= 1u << oldCellValue;
                    // }
                    // if (rowId == Rank - 1 - columnId)
                    // {
                    // 	flagsSecondary |= 1u << oldCellValue;
                    // }
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
                cellValue = squareA[rowId][columnId];
                // Возвращаем значение в служебные массивы
                if (cellValue != Square::Empty)
                {
                    // Возвращаем значение в столбцы
                    flagsColumns[columnId] |= 1u << cellValue;
                    // Возвращаем значение в строки
                    flagsRows[rowId] |= 1u << cellValue;
                    // Так как для 9 и 10 ранга в рамках маски workunit-а уже заданы обе диагонали, то проверку по ним - отключаем
                    // // Возвращаем значение в диагонали
                    // if (rowId == columnId)
                    // {
                    // 	flagsPrimary |= 1u << cellValue;
                    // }
                    // if (rowId == Rank - 1 - columnId)
                    // {
                    // 	flagsSecondary |= 1u << cellValue;
                    // }
                    // Сбрасываем клетку квадрата
                    squareA[rowId][columnId] = Square::Empty;
                    // Зачищаем историю клетки (rowId, columnId)
                    flagsCellsHistory[rowId][columnId] = (1u << Rank) - 1;
                }

                // Делаем шаг назад
                cellId--;
            }

            // Проверяем условие окончания поиска
            if (squareA[keyRowId][keyColumnId] == keyValue)
            {
                if (keyValue == Square::Empty)
                {
                    if (cellId < 0)
                    {
                        stop = Yes;
                    }
                }
                else
                {
                    stop = Yes;
                }
            }
        } while (!stop);
    }

    // Вывод итогов поиска
    ShowSearchTotals();
}

// Перетасовка строк заданного ДЛК в поиске ОДЛК к нему
void RakeSearch::PermuteRows()
{
    // Инициализация переменных заменила блок кода, приведённый под ней
    const unsigned int rowsHistoryMasks[Rank] = {
        1022, 1021, 1019, 1015, 1007, 991,
        959,  895,  767,  511}; // Набор битовых масок истории задействования строк, где i-е значение - это маска использования строк квадрата A на месте i-й строки квадрата B. При этом, каждый i-й бит i-й строки выставляется в 0 для того, чтобы исключить комбинации со строками, совпадающими хотя бы по одной строке с исходным квадратом. При зачистке истории строки (при возврате к предыдущей), возвращается значение из этого набора масок.
    int currentSquareRows[Rank] = {
        0,  -1, -1, -1, -1, -1,
        -1, -1, -1, -1}; // Массив с перечнем текущих строк, использованных в квадрате. На i-й позиции - номер используемой строки
    unsigned int primaryDiagonalFlags = 1022; // Битовые флаги свободных значений на основной диагонали
    unsigned int secondaryDiagonalFlags = 511; // Битовые флаги свободных значений на побочной диагонали
    unsigned int rowsUsageFlags = 1022; // Битовые флаги задействования строк в текущей комбинации
    unsigned int rowsHistoryFlags[Rank] = {
        1022, 1021, 1019, 1015, 1007,
        991,  959,  895,  767,  511}; // Битовые флаги истории задействования строк в рамках текущей комбинации
    int currentRowId = 1; // Номер обрабатываемой строки в формируемом квадрате

    int isRowGet = 0; // Флаг нахождения строки для проверки в рамках комбинации
    int oldRowId =
        -1; // Номер предыдущего вариант исходной строки строки формируемого квадрата, затираемой новым вариантом

    unsigned long insertedRowId = 0; // Номер новой строки для проверки совместимости с комбинацией
    unsigned int isRowFree =
        0; // Флаг наличия "свободной" для проверки в комбинации строки, которую потом можно проверять на диагональность
    unsigned int freeRowsMask =
        0; // Маска из соединения флагов задействования строк и истории задействования строк исходного квадрата для текущей строки новой комбинации
#ifdef __arm__
    unsigned int reversedMask; // Перевёрнутая битовая маска с флагами (для ARM, команды clz!)
#endif

    // Важная деталь! В состав битовых флагов истории зайдействования строк, уже включен фильтр, отсеивающий комбинацию, приводящую к исходному квадрату

    // Блок кода, заменённый инициализацией переменных
    //	// Сброс значений в векторах использования строк в очередной перестановке и номеров строк, использованных для текущего квадрата
    //	for (int i = 0; i < Rank; i++)
    //	{
    //		currentSquareRows[i] = -1;
    //	}
    //
    //	// Выставление флагов возможности использования значений на диагоналях
    //	primaryDiagonalFlags = (1u << Rank) - 1;
    //	secondaryDiagonalFlags = (1u << Rank) - 1;
    //	// Сброс флагов задействовая строк в новой комбинации
    //	rowsUsageFlags = (1u << Rank) - 1;
    //	// Сброс истории задействования строк
    //	for (int i = 0; i < Rank; i++)
    //	{
    //		rowsHistoryFlags[i] = (1u << Rank) - 1;
    //		rowsHistoryFlags[i] ^= 1u << i;
    //	}
    //
    //	// Выставляем флаги задействования первой (фиксированной) строки и её значений на диагоналях
    //	rowsUsageFlags ^= 1u;
    //	rowsHistoryFlags[0] ^= 1u;
    //	currentSquareRows[0] = 0;
    //	primaryDiagonalFlags ^= 1u << squareA[0][0];
    //	secondaryDiagonalFlags ^= 1u << squareA[0][Rank - 1];

    // Перебор различных комбинаций строк заданного ДЛК с целью поиска ОДЛК к нему
    // Сброс числа найденых пар для заданного ДЛК
    pairsCount = 0;

    // Инициализация переменных (отключенный блок кода)

    // Подбор всех возможных комбинаций из строк исходного квадрата
    while (currentRowId > 0)
    {
        // Очередной шаг в построении новой комбинации
        // Подбор строки из исходного квадрата на позицию currentRowId формирумого квадрата
        isRowGet = 0;

        // Поиск строки, не задействованной сейчас и ранее - для этой же строки
        do
        {
            // Проверка очередной строки на возможность добавления в новую комбинациию
            // Получение маски свободных строк
            freeRowsMask = rowsUsageFlags & rowsHistoryFlags[currentRowId];
// Получение номера свободной строки
#ifdef __arm__
            asm("rbit %1, %0" : "=r"(reversedMask) : "r"(freeRowsMask) :);
            asm("clz %1, %0" : "=r"(insertedRowId) : "r"(reversedMask) :);

            if (insertedRowId < Rank)
            {
                isRowFree = 1;
            }
            else
            {
                isRowFree = 0;
            }
#else
#ifdef _MSC_VER
            isRowFree = _BitScanForward(&insertedRowId, freeRowsMask);
#else
            insertedRowId = __builtin_ffs(freeRowsMask);
            if (insertedRowId > 0)
            {
                isRowFree = 1;
                insertedRowId--;
            }
            else
            {
                isRowFree = 0;
            }
#endif
#endif

            // Обработка найденной строки
            if (isRowFree)
            {
                // Проверка соместимости найденной строки по диагоналям с уже выбранными в комбинацию строками
                if ((primaryDiagonalFlags & (1u << squareA[insertedRowId][currentRowId])) &&
                    (secondaryDiagonalFlags & (1u << squareA[insertedRowId][Rank - currentRowId - 1])))
                {
                    isRowGet = 1;
                }

                // Отмечание строки с номером bitIndex как уже проверенной
                rowsHistoryFlags[currentRowId] ^= 1u << insertedRowId;
            }
        } while (isRowFree && !isRowGet);

        // Обработка найденной строки
        if (isRowGet)
        {
            // Добавление найденной строки в комбинацию строк формируемого квадрата
            // Считывание номера строки, уже присутствующей на этой позиции
            oldRowId = currentSquareRows[currentRowId];

            // Отмечание использования новой строки в массивах флагов
            // Сброс флага в массиве свободных/использованных строк
            rowsUsageFlags ^= 1u << insertedRowId;
            // Запись номера найденной строки в перечень строк новой комбинации
            currentSquareRows[currentRowId] = insertedRowId;
            // Отмечание диагональных значений из новой строки как использованных
            primaryDiagonalFlags ^= 1u << squareA[insertedRowId][currentRowId];
            secondaryDiagonalFlags ^= 1u << squareA[insertedRowId][Rank - currentRowId - 1];

            // Убирание флага занятости строки, убираемой из новой комбинации и значений её диагональных клеток
            if (oldRowId != -1)
            {
                rowsUsageFlags |= 1u << oldRowId;
                primaryDiagonalFlags |= 1u << squareA[oldRowId][currentRowId];
                secondaryDiagonalFlags |= 1u << squareA[oldRowId][Rank - currentRowId - 1];
            }

            // Обработка полученной комбинации
            if (currentRowId == Rank - 1)
            {
                // Обработка завершённой комбинации
                // Формирование квадрата по найденной комбинации
                for (int i = 0; i < Rank; i++)
                {
                    for (int j = 0; j < Rank; j++)
                    {
                        squareB[i][j] = squareA[currentSquareRows[i]][j];
                    }
                }

                // Обработка полученного квадрата
                ProcessOrthoSquare();
            }
            else
            {
                // Шаг вперёд - к следующей строк формируемого квадрата
                currentRowId++;
            }
        }
        else
        {
            // Шаг назад в поиске
            // Считывание номера строки исходного квадрата высвобождаемой, при уходе из текущей строки создаваемой комбинации
            oldRowId = currentSquareRows[currentRowId];
            // Отмечание текущей строки новой комбинации как незадействованной
            currentSquareRows[currentRowId] = -1;
            // Убирание флагов задействования строки из исходного квадрата и задействования её диагональных значений в новой комбинации
            if (oldRowId != -1)
            {
                rowsUsageFlags |= 1u << oldRowId;
                primaryDiagonalFlags |= 1u << squareA[oldRowId][currentRowId];
                secondaryDiagonalFlags |= 1u << squareA[oldRowId][Rank - currentRowId - 1];
            }
            // Очистка истории оставлямой строки комбинации
            // Первоначальная "зачистка" (выставление единиц во всех битах):
            // rowsHistoryFlags[currentRowId] = (1u << Rank) - 1;
            // заменена на восстановление маски с "обнулённой" строкой, соответствующей первоначальному квадрату
            rowsHistoryFlags[currentRowId] = rowsHistoryMasks[currentRowId];
            // Переход к предыдущей строке создаваемой комбинации
            currentRowId--;
        }
    }
}
