// Поиск пар диагональных латинских квадратов методом "перетасовки" строк

#include "RakeSearch.h"
#include <string.h>
#include <type_traits>

#ifdef HAS_SIMD
#ifdef __SSE2__
#include "immintrin.h"
#endif
#ifdef __ARM_NEON
#include "arm_neon.h"
#endif
#endif // HAS_SIMD

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

// Square::Empty is equal -1, all other values and non-negative.
// CPU sets sign bit in status register automatically when executing instructions,
// so sign check instead of value check can give faster code.
#define IsCellEmpty(val) ((val) < 0)

// Конструктор по умолчанию
RakeSearch::RakeSearch()
{
    Reset();
}

// Задание имен файлов параметров и контрольной точки
void RakeSearch::SetFileNames(const string& start, const string& result, const string& checkpoint, const string& temp)
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

    memset(squareA_Mask, 0, sizeof(squareA_Mask));
#if defined(HAS_SIMD) || defined(UT_BUILD)
    memset(squareA_MaskT, 0, sizeof(squareA_MaskT));
#endif

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
    flagsPrimary = AllFree;
    flagsSecondary = AllFree;

    // Сброс значений в матрицах использования элементов в столбцах и строках
    for (int i = 0; i < Rank; i++)
    {
        flagsColumns[i] = AllFree;
        flagsRows[i] = AllFree;
    }

    // Сброс значений в кубе истории использования значений в клетках
    for (int i = 0; i < Rank; i++)
    {
        for (int j = 0; j < Rank; j++)
        {
            flagsCellsHistory[i][j] = AllFree;
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
void RakeSearch::Initialize(const string& start, const string& result, const string& checkpoint, const string& temp)
{
    ifstream startFile;
    ifstream checkpointFile;

    // Считывание названий имен файлов
    startParametersFileName = start;
    resultFileName = result;
    checkpointFileName = checkpoint;
    tempCheckpointFileName = temp;

    // Считываем состояние генератора и поиска из файла контрольной точки или начальных значений
    // Открытие файлов со стартовыми параметрами и файла контрольной точки
    startFile.open(startParametersFileName.c_str(), std::ios_base::in);
    checkpointFile.open(checkpointFileName.c_str(), std::ios_base::in);

    Read(startFile);
    array<int, MaxPathPrefixes> tmpPrefixes;
    GeneratePathPrefixes(tmpPrefixes, 0);
    startFile.seekg(0);

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
    else
    {
        isStartFromCheckpoint = 0;
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

void RakeSearch::GeneratePathPrefixes(array<int, MaxPathPrefixes>& tmp, int pathPos)
{
    if (MaxPathPrefixes == pathPos)
    {
        pathPrefixes.emplace_back(tmp);
    }
    else
    {
        int r = path[pathPos][0];
        int c = path[pathPos][1];

        int rh = flagsRows[r];
        int ch = flagsColumns[c];
        int hh = flagsCellsHistory[r][c];

        int freeVals = rh & ch & hh;

        for (int n = 0; n < Rank; ++n)
        {
            int m = 1 << n;
            if (0 != (freeVals & m))
            {
                flagsRows[r] &= ~m;
                flagsColumns[c] &= ~m;
                flagsCellsHistory[r][c] &= ~m;

                tmp[pathPos] = __builtin_ctz(m);

                GeneratePathPrefixes(tmp, pathPos + 1);

                flagsRows[r] |= m;
                flagsColumns[c] |= m;
            }
        }

        flagsRows[r] = rh;
        flagsColumns[c] = ch;
        flagsCellsHistory[r][c] = hh;
    }
}

// Чтение состояния поиска из потока
void RakeSearch::Read(istream& is)
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
        flagsPrimary = 0;
        for (int i = 0; i < Rank; i++)
        {
            is >> storedBit;
            if (storedBit)
            {
                SetBit(flagsPrimary, i);
            }
        }

        // Считывание информации о значениях на побочной диагонали
        flagsSecondary = 0;
        for (int i = 0; i < Rank; i++)
        {
            is >> storedBit;
            if (storedBit)
            {
                SetBit(flagsSecondary, i);
            }
        }

        // Считывание информации о значениях в строках
        memset(flagsRows, 0, sizeof(flagsRows));
        for (int i = 0; i < Rank; i++)
        {
            for (int j = 0; j < Rank; j++)
            {
                is >> storedBit;
                if (storedBit)
                {
                    SetBit(flagsRows[i], j);
                }
            }
        }

        // Считывание информации о значениях в столбцах
        memset(flagsColumns, 0, sizeof(flagsColumns));
        for (int i = 0; i < Rank; i++)
        {
            for (int j = 0; j < Rank; j++)
            {
                is >> storedBit;
                if (storedBit)
                {
                    SetBit(flagsColumns[i], j);
                }
            }
        }

        // Считывание информации об истории значений в клетках квадрата
        memset(flagsCellsHistory, 0, sizeof(flagsCellsHistory));
        for (int h = 0; h < Rank; h++)
        {
            for (int i = 0; i < Rank; i++)
            {
                for (int j = 0; j < Rank; j++)
                {
                    is >> storedBit;
                    if (storedBit)
                    {
                        SetBit(flagsCellsHistory[i][j], h);
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

    // Data loaded. Perform necessary post-loading tasks.
    if (cellId == cellsInPath - 1)
    {
        // Start from WU
        // Convert old checkpoint format to new one if used
        int row = path[cellsInPath - 2][0], col = path[cellsInPath - 2][1];
        if (0 != flagsCellsHistory[row][col])
        {
            int tmpColumns[Rank];
            int tmpRows[Rank];
            memcpy(tmpColumns, flagsColumns, sizeof(flagsColumns));
            memcpy(tmpRows, flagsRows, sizeof(flagsRows));

            // Convert cellsHistory into candidates
            for (int i = cellsInPath - 1; i >= 0; --i)
            {
                row = path[i][0];
                col = path[i][1];
                int bit = 1 << squareA[row][col];
                tmpColumns[col] |= bit;
                tmpRows[row] |= bit;
                flagsCellsHistory[row][col] &= tmpColumns[col] & tmpRows[row];

                // Update rows/cols data for last cell in path, it is no longer set
                if (i == cellsInPath - 1)
                {
                    flagsColumns[col] = tmpColumns[col];
                    flagsRows[row] = tmpRows[row];
                }
            }
        }
    }
    else
    {
        // Start from checkpoint
        // Check if there are no cells on diagonals in path
        for (int i = 0; i < cellsInPath; i++)
        {
            int row = path[i][0], col = path[i][1];
            if ((row == col) || (row == Rank - 1 - col))
            {
                std::cerr << "Error: Cell on diagonal in path! R=" << row << " C=" << col << std::endl;
                return;
            }
        }
    }
}

// Запись состояния поиска в поток
void RakeSearch::Write(std::ostream& os)
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
        os << GetBit01(flagsPrimary, i) << " ";
    }
    os << endl;

    // Запись информации о значениях на побочной диагонали
    for (int i = 0; i < Rank; i++)
    {
        os << GetBit01(flagsSecondary, i) << " ";
    }
    os << endl;

    // Дополнительная пустая строка
    os << endl;

    // Запись информации о значениях в строках
    for (int i = 0; i < Rank; i++)
    {
        for (int j = 0; j < Rank; j++)
        {
            os << GetBit01(flagsRows[i], j) << " ";
        }
        os << endl;
    }
    os << endl;

    // Запись информации о значениях в столбцах
    for (int i = 0; i < Rank; i++)
    {
        for (int j = 0; j < Rank; j++)
        {
            os << GetBit01(flagsColumns[i], j) << " ";
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
                os << GetBit01(flagsCellsHistory[i][j], h) << " ";
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

        // The stream for output into the results file
        // It must be here, creation and destruction of iostream is costly!
        ofstream resultFile;
        resultFile.open(resultFileName.c_str(), std::ios_base::binary | std::ios_base::app);
        if (!resultFile.is_open())
        {
            std::cerr << "Error opening file!";
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
            if (resultFile.is_open())
            {
                resultFile << "{" << endl;
                resultFile << "# ------------------------" << endl;
                resultFile << "# Detected pair for the square: " << endl;
                resultFile << "# Degree of orthogonality: " << orthoDegree << endl;
                resultFile << "# ------------------------" << endl;
                resultFile << a;
                resultFile << "# ------------------------" << endl;
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

    pairsCount = 0;

    GenerateSquareMasks();

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
        for (auto it = pathPrefixes.begin() + pathPrefixPos; it != pathPrefixes.end(); ++it)
        {
            const auto& prefix = *it;
            bool greaterOrEqual = true;
            for (int n = 0; n < MaxPathPrefixes; ++n)
            {
                if (prefix[n] != squareA[path[n][0]][path[n][1]])
                {
                    greaterOrEqual = prefix[n] >= squareA[path[n][0]][path[n][1]];
                    break;
                }
            }
            if (greaterOrEqual)
                break;
            ++pathPrefixPos;
        }
        fraction_done = pathPrefixPos / (double)pathPrefixes.size();

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
            cout << "# Done: " << pathPrefixPos << "/" << pathPrefixes.size() << " = " << fraction_done * 100.0 << "%"
                 << endl;
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

// Start the squares generation
void RakeSearch::Start()
{
    // Check value of keyValue and pass result as a type to StartImpl
    if (IsCellEmpty(keyValue))
        StartImpl<true_type>();
    else
        StartImpl<false_type>();

    // Вывод итогов поиска
    ShowSearchTotals();
}

// Actual implementation of the squares generation
// Note: values on diagonal are preset in WU, so corresponding parts of code are commented out.
// It turned out that it was quite costly to have instructions which were doing nothing.
template <typename IsKeyValueEmpty> inline void RakeSearch::StartImpl()
{
    int cellValue;           // New value for the cell
    int cellValueCandidates; // Candidates for value for the cell

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

    // Selection of the value for the next cell
    // Read coordinates of the cell
    rowId = path[cellId][0];
    columnId = path[cellId][1];

    // Generate new value for the cell (rowId, columnId)
    // Select the value for the cell
    // Check the i value for possibility to be written into the cell (rowId, columnId)
    cellValueCandidates = flagsColumns[columnId] & flagsRows[rowId];

    if (isInitialized == Yes)
    {
        // Check if there are no candidates at the beginning, or if calculations are resumed from checkpoint
        if ((cellId == cellsInPath - 1) || (0 == cellValueCandidates))
            goto StepDown;

        // Selection of the cells values
        while (1)
        {
            // Process the search result
            // 1st loop (used to be "if (cellValueCandidates)" part) - handle case when at least one cell value candidate is present
            while (1)
            {
                // Extract lowest bit set
                int bit = (-cellValueCandidates) & cellValueCandidates;

                // Write the value into the square
                squareA[rowId][columnId] = __builtin_ctz(bit);

                // Process the finish of the square generation
                if (cellId == cellsInPath - 1)
                {
                    // Process the found square
                    ProcessSquare();

                    // Check the finish condition of search
                    if (!IsKeyValueEmpty::value)
                    {
                        // Set the flag if the terminal value is other
                        if (squareA[keyRowId][keyColumnId] == keyValue)
                        {
                            break;
                        }
                    }

                    break;
                }
                else
                {
                    // Mark the value in columns
                    flagsColumns[columnId] &= ~bit;
                    // Mark the value in rows
                    flagsRows[rowId] &= ~bit;

                    // Mark the value in the history of cell values
                    flagsCellsHistory[rowId][columnId] = cellValueCandidates & ~bit;

                    // Step forward
                    cellId++;

                    // Check the finish condition of search
                    if (!IsKeyValueEmpty::value)
                    {
                        // Set the flag if the terminal value is other
                        if (squareA[keyRowId][keyColumnId] == keyValue)
                        {
                            break;
                        }
                    }

                    // Selection of the value for the next cell
                    // Read coordinates of the cell
                    rowId = path[cellId][0];
                    columnId = path[cellId][1];

                    // Generate new value for the cell (rowId, columnId)
                    // Select the value for the cell
                    // Check the i value for possibility to be written into the cell (rowId, columnId)
                    cellValueCandidates = flagsColumns[columnId] & flagsRows[rowId];

                    if (!cellValueCandidates)
                        break;
                }
            }

            // 2nd loop (used to be "else" part) - handle case when there are no cell value candidates
        StepDown:
            while (1)
            {
                // Step backward
                cellId--;

                // Check the finish condition of search
                if (IsKeyValueEmpty::value)
                {
                    // Set the flag if the terminal value is "-1" which means we must leave the cell
                    if (cellId < 0 /*&& IsCellEmpty(newSquare.Matrix[keyRowId][keyColumnId])*/)
                    {
                        return;
                    }
                }

                // Selection of the value for the next cell
                // Read coordinates of the cell
                rowId = path[cellId][0];
                columnId = path[cellId][1];

                // Process the fact of not-founding a new value in the cell (rowId; columnId)
                // Restore the previous value from the square into arrays
                // Read the current value
                cellValue = squareA[rowId][columnId];

                // Restore the value into auxilary arrays
                // Restore the value into columns
                SetFree(flagsColumns[columnId], cellValue);
                // Restore the value into rows
                SetFree(flagsRows[rowId], cellValue);

                cellValueCandidates = flagsCellsHistory[rowId][columnId];

                if (cellValueCandidates)
                    break;
            }
        }
    }
}

#if defined(__ARM_NEON) && !defined(__aarch64__) && defined(HAS_SIMD)
__attribute__((always_inline)) inline void RakeSearch::transposeMatrix4x4(int srcRow, int srcCol, int destRow,
                                                                          int destCol)
{
    uint16x4_t v1, v2;
    v1 = vld1_u16((uint16_t*)(&squareA_Mask[srcRow + 0][srcCol + 0]));
    v2 = vld1_u16((uint16_t*)(&squareA_Mask[srcRow + 0][srcCol + 2]));
    uint16x4_t v1_1 = vuzp_u16(v1, v2).val[0];
    v1 = vld1_u16((uint16_t*)(&squareA_Mask[srcRow + 1][srcCol + 0]));
    v2 = vld1_u16((uint16_t*)(&squareA_Mask[srcRow + 1][srcCol + 2]));
    uint16x4_t v2_1 = vuzp_u16(v1, v2).val[0];
    v1 = vld1_u16((uint16_t*)(&squareA_Mask[srcRow + 2][srcCol + 0]));
    v2 = vld1_u16((uint16_t*)(&squareA_Mask[srcRow + 2][srcCol + 2]));
    uint16x4_t v3_1 = vuzp_u16(v1, v2).val[0];
    v1 = vld1_u16((uint16_t*)(&squareA_Mask[srcRow + 3][srcCol + 0]));
    v2 = vld1_u16((uint16_t*)(&squareA_Mask[srcRow + 3][srcCol + 2]));
    uint16x4_t v4_1 = vuzp_u16(v1, v2).val[0];

    uint16x4x2_t v12_2 = vtrn_u16(v1_1, v2_1);
    uint16x4x2_t v34_2 = vtrn_u16(v3_1, v4_1);

    uint32x2x2_t v13_3 = vtrn_u32(vreinterpret_u32_u16(v12_2.val[0]), vreinterpret_u32_u16(v34_2.val[0]));
    uint32x2x2_t v24_3 = vtrn_u32(vreinterpret_u32_u16(v12_2.val[1]), vreinterpret_u32_u16(v34_2.val[1]));

    vst1_u32((uint32_t*)(&squareA_MaskT[destRow + 0][destCol + 0]), v13_3.val[0]);
    vst1_u32((uint32_t*)(&squareA_MaskT[destRow + 1][destCol + 0]), v24_3.val[0]);
    vst1_u32((uint32_t*)(&squareA_MaskT[destRow + 2][destCol + 0]), v13_3.val[1]);
    vst1_u32((uint32_t*)(&squareA_MaskT[destRow + 3][destCol + 0]), v24_3.val[1]);
}
#endif

void RakeSearch::GenerateSquareMasks()
{
    // Generate bitmasks
#if defined(__AVX2__) && defined(HAS_SIMD)
    // AVX2 has "shift by vector" instruction, use it here
    // Note: AVX512 instructions which use ZMM registers cause too big
    // CPU frequency throttling. It does not make sense to use them in this
    // one place only, as everything else will be slowed down too.
    int n = 0;
    for (; n < Rank * Rank - 7; n += 8)
    {
        __m256i v = _mm256_load_si256((__m256i*)(&squareA[0][0] + n));
        v = _mm256_sllv_epi32(_mm256_set1_epi32(1), v);
        _mm256_store_si256((__m256i*)(&squareA_Mask[0][0] + n), v);
    }
    // Use SSE instruction if possible at the end
    if ((Rank * Rank) % 8 >= 4)
    {
        __m128i v = _mm_load_si128((__m128i*)(&squareA[0][0] + n));
        v = _mm_sllv_epi32(_mm_set1_epi32(1), v);
        _mm_store_si128((__m128i*)(&squareA_Mask[0][0] + n), v);

        n += 4;
    }
    // Process remaining elements
    if ((Rank * Rank) % 4 > 0)
    {
        for (; n < Rank * Rank; n++)
        {
            int x = *(&squareA[0][0] + n);
            *((&squareA_Mask[0][0] + n)) = 1 << x;
        }
    }
#elif defined(__SSSE3__) && defined(HAS_SIMD)
    // SSSE3 added shuffle instruction, which can be used to build small lookup table.
    // Maximum val needs more than 8 bits, so some extra check and shift by constant is
    // required. This is still faster than unvectorized code.
    const __m128i vcLut = _mm_set_epi8(128, 64, 32, 16, 8, 4, 2, 1, 128, 64, 32, 16, 8, 4, 2, 1);
    const __m128i vc0 = _mm_setzero_si128();
    const __m128i vc8 = _mm_set1_epi16(8);
    int n = 0;
    for (; n < Rank * Rank - 7; n += 8)
    {
        // Load data
        __m128i v1 = _mm_load_si128((__m128i*)(&squareA[0][0] + n));
        __m128i v2 = _mm_load_si128((__m128i*)(&squareA[0][0] + n + 4));

        // Pack two 32x4 vectors into one 8x16
        v1 = _mm_packs_epi32(v1, v2);
        __m128i v_lut_idx = _mm_packs_epi16(v1, vc0);
        // Get mask fom LUT
        __m128i v_lo = _mm_shuffle_epi8(vcLut, v_lut_idx);
        // Convert vector 8x16 to 16x8
        v_lo = _mm_unpacklo_epi8(v_lo, vc0);

        // Check for numbers >= 8, and prepare mask for them
        __m128i v_cmp_lt8 = _mm_cmplt_epi16(v1, vc8);
        __m128i v_hi = _mm_slli_epi16(v_lo, 8);

        // Create resulting vector
#ifdef __SSE4_1__
        v1 = _mm_blendv_epi8(v_hi, v_lo, v_cmp_lt8);
#else
        v1 = _mm_or_si128(_mm_and_si128(v_cmp_lt8, v_lo), _mm_andnot_si128(v_cmp_lt8, v_hi));
#endif

        // Convert vector 16x8 into two 32x4, and store results
        v2 = _mm_unpackhi_epi16(v1, vc0);
        v1 = _mm_unpacklo_epi16(v1, vc0);

        _mm_store_si128((__m128i*)(&squareA_Mask[0][0] + n), v1);
        _mm_store_si128((__m128i*)(&squareA_Mask[0][0] + n + 4), v2);
    }
    // Process remaining elements
    for (; n < Rank * Rank; n++)
    {
        int x = *(&squareA[0][0] + n);
        *((&squareA_Mask[0][0] + n)) = 1 << x;
    }
#else
    // Default non-SIMD code
    // Note: this will be autovectorized for ARM NEON. gcc has limit how many times
    // it can unroll loop, so single loop with manual vectorization would be slower.
    // Two nested loops are below limit, so autovectorization creates expected
    // machine code here.
    for (int i = 0; i < Rank; i++)
    {
        for (int j = 0; j < Rank; j++)
        {
            squareA_Mask[i][j] = 1u << squareA[i][j];
        }
    }
#endif

    // Create transposed copy of squareA_Mask if needed
#if defined(__SSE2__) && defined(HAS_SIMD)
    __m128i v1, v2;
    v1 = _mm_loadu_si128((__m128i*)(&squareA_Mask[0][0]));
    v2 = _mm_loadu_si128((__m128i*)(&squareA_Mask[0][4]));
    __m128i v1_1 = _mm_packs_epi32(v1, v2);
    v1 = _mm_loadu_si128((__m128i*)(&squareA_Mask[1][0]));
    v2 = _mm_loadu_si128((__m128i*)(&squareA_Mask[1][4]));
    __m128i v2_1 = _mm_packs_epi32(v1, v2);
    v1 = _mm_loadu_si128((__m128i*)(&squareA_Mask[2][0]));
    v2 = _mm_loadu_si128((__m128i*)(&squareA_Mask[2][4]));
    __m128i v3_1 = _mm_packs_epi32(v1, v2);
    v1 = _mm_loadu_si128((__m128i*)(&squareA_Mask[3][0]));
    v2 = _mm_loadu_si128((__m128i*)(&squareA_Mask[3][4]));
    __m128i v4_1 = _mm_packs_epi32(v1, v2);
    v1 = _mm_loadu_si128((__m128i*)(&squareA_Mask[4][0]));
    v2 = _mm_loadu_si128((__m128i*)(&squareA_Mask[4][4]));
    __m128i v5_1 = _mm_packs_epi32(v1, v2);
    v1 = _mm_loadu_si128((__m128i*)(&squareA_Mask[5][0]));
    v2 = _mm_loadu_si128((__m128i*)(&squareA_Mask[5][4]));
    __m128i v6_1 = _mm_packs_epi32(v1, v2);
    v1 = _mm_loadu_si128((__m128i*)(&squareA_Mask[6][0]));
    v2 = _mm_loadu_si128((__m128i*)(&squareA_Mask[6][4]));
    __m128i v7_1 = _mm_packs_epi32(v1, v2);
    v1 = _mm_loadu_si128((__m128i*)(&squareA_Mask[7][0]));
    v2 = _mm_loadu_si128((__m128i*)(&squareA_Mask[7][4]));
    __m128i v8_1 = _mm_packs_epi32(v1, v2);

    __m128i v1_2 = _mm_unpacklo_epi16(v1_1, v2_1);
    __m128i v2_2 = _mm_unpackhi_epi16(v1_1, v2_1);
    __m128i v3_2 = _mm_unpacklo_epi16(v3_1, v4_1);
    __m128i v4_2 = _mm_unpackhi_epi16(v3_1, v4_1);
    __m128i v5_2 = _mm_unpacklo_epi16(v5_1, v6_1);
    __m128i v6_2 = _mm_unpackhi_epi16(v5_1, v6_1);
    __m128i v7_2 = _mm_unpacklo_epi16(v7_1, v8_1);
    __m128i v8_2 = _mm_unpackhi_epi16(v7_1, v8_1);

    __m128i v1_3 = _mm_unpacklo_epi32(v1_2, v3_2);
    __m128i v2_3 = _mm_unpackhi_epi32(v1_2, v3_2);
    __m128i v3_3 = _mm_unpacklo_epi32(v2_2, v4_2);
    __m128i v4_3 = _mm_unpackhi_epi32(v2_2, v4_2);
    __m128i v5_3 = _mm_unpacklo_epi32(v5_2, v7_2);
    __m128i v6_3 = _mm_unpackhi_epi32(v5_2, v7_2);
    __m128i v7_3 = _mm_unpacklo_epi32(v6_2, v8_2);
    __m128i v8_3 = _mm_unpackhi_epi32(v6_2, v8_2);

    __m128i v1_4 = _mm_unpacklo_epi64(v1_3, v5_3);
    __m128i v2_4 = _mm_unpackhi_epi64(v1_3, v5_3);
    __m128i v3_4 = _mm_unpacklo_epi64(v2_3, v6_3);
    __m128i v4_4 = _mm_unpackhi_epi64(v2_3, v6_3);
    __m128i v5_4 = _mm_unpacklo_epi64(v3_3, v7_3);
    __m128i v6_4 = _mm_unpackhi_epi64(v3_3, v7_3);
    __m128i v7_4 = _mm_unpacklo_epi64(v4_3, v8_3);
    __m128i v8_4 = _mm_unpackhi_epi64(v4_3, v8_3);

    _mm_store_si128((__m128i*)(&squareA_MaskT[0][0]), v1_4);
    _mm_store_si128((__m128i*)(&squareA_MaskT[1][0]), v2_4);
    _mm_store_si128((__m128i*)(&squareA_MaskT[2][0]), v3_4);
    _mm_store_si128((__m128i*)(&squareA_MaskT[3][0]), v4_4);
    _mm_store_si128((__m128i*)(&squareA_MaskT[4][0]), v5_4);
    _mm_store_si128((__m128i*)(&squareA_MaskT[5][0]), v6_4);
    _mm_store_si128((__m128i*)(&squareA_MaskT[6][0]), v7_4);
    _mm_store_si128((__m128i*)(&squareA_MaskT[7][0]), v8_4);

    // Transpose data from last columns (excluding bottom-right part)
    for (int i = 0; i < 8; i++)
    {
        for (int j = 8; j < Rank; j++)
        {
            squareA_MaskT[j][i] = squareA_Mask[i][j];
        }
    }
    // Transpose data from last rows
    for (int i = 8; i < Rank; i++)
    {
        for (int j = 0; j < Rank; j++)
        {
            squareA_MaskT[j][i] = squareA_Mask[i][j];
        }
    }
#elif defined(__ARM_NEON) && defined(HAS_SIMD)
#ifdef __aarch64__
    uint16x8_t v1, v2;
    v1 = vld1q_u16((uint16_t*)(&squareA_Mask[0][0]));
    v2 = vld1q_u16((uint16_t*)(&squareA_Mask[0][4]));
    uint16x8_t v1_1 = vuzp1q_u16(v1, v2);
    v1 = vld1q_u16((uint16_t*)(&squareA_Mask[1][0]));
    v2 = vld1q_u16((uint16_t*)(&squareA_Mask[1][4]));
    uint16x8_t v2_1 = vuzp1q_u16(v1, v2);
    v1 = vld1q_u16((uint16_t*)(&squareA_Mask[2][0]));
    v2 = vld1q_u16((uint16_t*)(&squareA_Mask[2][4]));
    uint16x8_t v3_1 = vuzp1q_u16(v1, v2);
    v1 = vld1q_u16((uint16_t*)(&squareA_Mask[3][0]));
    v2 = vld1q_u16((uint16_t*)(&squareA_Mask[3][4]));
    uint16x8_t v4_1 = vuzp1q_u16(v1, v2);
    v1 = vld1q_u16((uint16_t*)(&squareA_Mask[4][0]));
    v2 = vld1q_u16((uint16_t*)(&squareA_Mask[4][4]));
    uint16x8_t v5_1 = vuzp1q_u16(v1, v2);
    v1 = vld1q_u16((uint16_t*)(&squareA_Mask[5][0]));
    v2 = vld1q_u16((uint16_t*)(&squareA_Mask[5][4]));
    uint16x8_t v6_1 = vuzp1q_u16(v1, v2);
    v1 = vld1q_u16((uint16_t*)(&squareA_Mask[6][0]));
    v2 = vld1q_u16((uint16_t*)(&squareA_Mask[6][4]));
    uint16x8_t v7_1 = vuzp1q_u16(v1, v2);
    v1 = vld1q_u16((uint16_t*)(&squareA_Mask[7][0]));
    v2 = vld1q_u16((uint16_t*)(&squareA_Mask[7][4]));
    uint16x8_t v8_1 = vuzp1q_u16(v1, v2);

    uint16x8_t v1_2 = vtrn1q_u16(v1_1, v2_1);
    uint16x8_t v2_2 = vtrn2q_u16(v1_1, v2_1);
    uint16x8_t v3_2 = vtrn1q_u16(v3_1, v4_1);
    uint16x8_t v4_2 = vtrn2q_u16(v3_1, v4_1);
    uint16x8_t v5_2 = vtrn1q_u16(v5_1, v6_1);
    uint16x8_t v6_2 = vtrn2q_u16(v5_1, v6_1);
    uint16x8_t v7_2 = vtrn1q_u16(v7_1, v8_1);
    uint16x8_t v8_2 = vtrn2q_u16(v7_1, v8_1);

    uint32x4_t v1_3 = vtrn1q_u32(vreinterpretq_u32_u16(v1_2), vreinterpretq_u32_u16(v3_2));
    uint32x4_t v2_3 = vtrn1q_u32(vreinterpretq_u32_u16(v2_2), vreinterpretq_u32_u16(v4_2));
    uint32x4_t v3_3 = vtrn2q_u32(vreinterpretq_u32_u16(v1_2), vreinterpretq_u32_u16(v3_2));
    uint32x4_t v4_3 = vtrn2q_u32(vreinterpretq_u32_u16(v2_2), vreinterpretq_u32_u16(v4_2));
    uint32x4_t v5_3 = vtrn1q_u32(vreinterpretq_u32_u16(v5_2), vreinterpretq_u32_u16(v7_2));
    uint32x4_t v6_3 = vtrn1q_u32(vreinterpretq_u32_u16(v6_2), vreinterpretq_u32_u16(v8_2));
    uint32x4_t v7_3 = vtrn2q_u32(vreinterpretq_u32_u16(v5_2), vreinterpretq_u32_u16(v7_2));
    uint32x4_t v8_3 = vtrn2q_u32(vreinterpretq_u32_u16(v6_2), vreinterpretq_u32_u16(v8_2));

    uint64x2_t v1_4 = vtrn1q_u64(vreinterpretq_u64_u32(v1_3), vreinterpretq_u64_u32(v5_3));
    uint64x2_t v2_4 = vtrn1q_u64(vreinterpretq_u64_u32(v2_3), vreinterpretq_u64_u32(v6_3));
    uint64x2_t v3_4 = vtrn1q_u64(vreinterpretq_u64_u32(v3_3), vreinterpretq_u64_u32(v7_3));
    uint64x2_t v4_4 = vtrn1q_u64(vreinterpretq_u64_u32(v4_3), vreinterpretq_u64_u32(v8_3));
    uint64x2_t v5_4 = vtrn2q_u64(vreinterpretq_u64_u32(v1_3), vreinterpretq_u64_u32(v5_3));
    uint64x2_t v6_4 = vtrn2q_u64(vreinterpretq_u64_u32(v2_3), vreinterpretq_u64_u32(v6_3));
    uint64x2_t v7_4 = vtrn2q_u64(vreinterpretq_u64_u32(v3_3), vreinterpretq_u64_u32(v7_3));
    uint64x2_t v8_4 = vtrn2q_u64(vreinterpretq_u64_u32(v4_3), vreinterpretq_u64_u32(v8_3));

    vst1q_u64((uint64_t*)(&squareA_MaskT[0][0]), v1_4);
    vst1q_u64((uint64_t*)(&squareA_MaskT[1][0]), v2_4);
    vst1q_u64((uint64_t*)(&squareA_MaskT[2][0]), v3_4);
    vst1q_u64((uint64_t*)(&squareA_MaskT[3][0]), v4_4);
    vst1q_u64((uint64_t*)(&squareA_MaskT[4][0]), v5_4);
    vst1q_u64((uint64_t*)(&squareA_MaskT[5][0]), v6_4);
    vst1q_u64((uint64_t*)(&squareA_MaskT[6][0]), v7_4);
    vst1q_u64((uint64_t*)(&squareA_MaskT[7][0]), v8_4);

    // Transpose data from last columns (excluding bottom-right part)
    for (int i = 0; i < 8; i++)
    {
        for (int j = 8; j < Rank; j++)
        {
            squareA_MaskT[j][i] = squareA_Mask[i][j];
        }
    }
    // Transpose data from last rows
    for (int i = 8; i < Rank; i++)
    {
        for (int j = 0; j < Rank; j++)
        {
            squareA_MaskT[j][i] = squareA_Mask[i][j];
        }
    }
#else  // !__aarch64__
    transposeMatrix4x4(0, 0, 0, 0);
    transposeMatrix4x4(0, 4, 4, 0);
    transposeMatrix4x4(4, 0, 0, 4);
    transposeMatrix4x4(4, 4, 4, 4);

    static_assert(Rank < 12, "Use transposeMatrix4x4 more times to optimize transpose for Rank 12+");

    // Transpose data from last columns (excluding bottom-right part)
    for (int i = 0; i < 8; i++)
    {
        for (int j = 8; j < Rank; j++)
        {
            squareA_MaskT[j][i] = squareA_Mask[i][j];
        }
    }
    // Transpose data from last rows
    for (int i = 8; i < Rank; i++)
    {
        for (int j = 0; j < Rank; j++)
        {
            squareA_MaskT[j][i] = squareA_Mask[i][j];
        }
    }
#endif // !__aarch64__
#else
    // Non-SIMD code does not use squareA_MaskT, so nothing here

#ifdef UT_BUILD
    // Transpose data for tests
    for (int i = 0; i < Rank; i++)
    {
        for (int j = 0; j < Rank; j++)
        {
            squareA_MaskT[j][i] = squareA_Mask[i][j];
        }
    }
#endif

#endif // __SSE2__
}

// Permute the rows of the given DLS, trying to find ODLS for it
void RakeSearch::PermuteRows()
{
    static_assert(Rank <= 16, "Function needs update for Rank 17+");

    // Masks for rowUsage which exclude current row for rows beside 1st one.
    // This is done to prevent generation of squares which has some rows unpermuted.
    // Every such row reduces final OrthoDegree of square pair by Rank, so squares
    // with such rows are eliminated early. This also improves performance.
    const unsigned int rowsUsageMasks[RankAligned] = {
        0xFFFF & ~0x0001, 0xFFFF & ~0x0002, 0xFFFF & ~0x0004, 0xFFFF & ~0x0008, 0xFFFF & ~0x0010, 0xFFFF & ~0x0020,
        0xFFFF & ~0x0040, 0xFFFF & ~0x0080, 0xFFFF & ~0x0100, 0xFFFF & ~0x0200, 0xFFFF & ~0x0400, 0xFFFF & ~0x0800,
        0xFFFF & ~0x1000, 0xFFFF & ~0x2000, 0xFFFF & ~0x4000, 0xFFFF & ~0x8000};

    int currentSquareRows[Rank];
    unsigned int rowsHistoryFlags[Rank];

    int currentRowId;
    int gettingRowId = -1;

    int diagonalValues1, diagonalValues2;

    int diagonalValuesHistory[Rank][2];

    int rowsUsage; // Flags of the rows usage at the current moment; rowsUsage[number of the row] = 0 | 1, where 0 means the row is already used, 1 - not.

    int rowCandidates; // Rows which still have to be checked

    // Mark the usage of the 1st row, because it is fixed
    rowsUsage = AllFree & ~1u;
    currentSquareRows[0] = 0;

#ifndef HAS_SIMD
    // For non-vectorized builds start from 1st row, and mark all rows except 0th as candidates
    currentRowId = 1;
    rowCandidates = AllFree & ~1u;
#else
    // SSE2/AVX2 version performs duplicate check when it steps forward to next row
    // and saves result as a candidates. So we have to start from 0th row and
    // mark it as the only candidate in order to perform duplicate check for 1st row.
    currentRowId = 0;
    rowCandidates = 1;
#endif

    // Set bits for diagonal values in 1st row. 1st row always has values 0,1,2,3...
    diagonalValues1 = 1;
    diagonalValues2 = 1u << (Rank - 1);

    diagonalValuesHistory[0][0] = diagonalValues1;
    diagonalValuesHistory[0][1] = diagonalValues2;

#if defined(__ARM_NEON) && defined(HAS_SIMD)
    // Set the powers of 2
    const uint16_t powersOf2[16] = {0x0001, 0x0002, 0x0004, 0x0008, 0x0010, 0x0020, 0x0040, 0x0080,
                                    0x0100, 0x0200, 0x0400, 0x0800, 0x1000, 0x2000, 0x4000, 0x8000};
#ifdef __aarch64__
    const uint16x8_t vPowersOf2_1 = vld1q_u16(powersOf2);
    const uint16x8_t vPowersOf2_2 = vld1q_u16(powersOf2 + 8);
#else
    const uint16x4_t vPowersOf2_1 = vld1_u16(powersOf2);
    const uint16x4_t vPowersOf2_2 = vld1_u16(powersOf2 + 4);
    const uint16x4_t vPowersOf2_3 = vld1_u16(powersOf2 + 8);
    static_assert(Rank <= 12, "ARM NEON code needs more vector(s) for Ranks 13+");
#endif
#endif

    // Note: code below assumes that rowCandidates is non-zero at beginning
    // so 1st nested while loop should execute first. If it may not be the case,
    // change code to handle this.
    // 1st loop (used to be "if (rowCandidates)" part) - handle case when at least one row candidate is present
    while (1)
    {
        // Select a row from the initial square for the position currentRowId of the generated square
        // Process the search result
        while (1)
        {
#ifndef HAS_SIMD
            int rowCandidatesMasked = rowCandidates & rowsUsageMasks[currentRowId];
            int bit1, bit2;
            bool foundCandidate = false;
            while (rowCandidatesMasked)
            {
                gettingRowId = __builtin_ctz(rowCandidatesMasked);
                // Process the new found row

                // Check diagonality of the generated part of the square
                // Check the main diagonal and secondary diagonal
                // Get bits for current row
                bit1 = squareA_Mask[gettingRowId][currentRowId];
                bit2 = squareA_Mask[gettingRowId][Rank - 1 - currentRowId];

                // Duplicate check
                int duplicationDetected = (diagonalValues1 & bit1) | (diagonalValues2 & bit2);

                // Process the results of checking the square for diagonality
                if (duplicationDetected)
                {
                    // Mark the row in the history of the used rows
                    ClearBit(rowCandidates, gettingRowId);
                    ClearBit(rowCandidatesMasked, gettingRowId);
                }
                else
                {
                    foundCandidate = true;
                    break;
                }
            }
            if (!foundCandidate)
                break;
#else // HAS_SIMD
            gettingRowId = __builtin_ctz(rowCandidates);
            // Process the new found row

            // Check diagonality of the generated part of the square
            // Check the main diagonal and secondary diagonal
            // Get bits for current row
            int bit1 = squareA_Mask[gettingRowId][currentRowId];
            int bit2 = squareA_Mask[gettingRowId][Rank - 1 - currentRowId];
#endif

            // Mark the row in the history of the used rows
            ClearBit(rowCandidates, gettingRowId);

            // Write the row into the array of the current rows
            currentSquareRows[currentRowId] = gettingRowId;

            // Step forward depending on the current position
            if (currentRowId == Rank - 1)
            {
                // Write rows into the square in correct order
                for (int n = 0; n < Rank; ++n)
                {
                    memcpy(&squareB[n][0], &squareA[currentSquareRows[n]][0], Rank * sizeof(squareB[n][0]));
                }

                // Process the found square
                ProcessOrthoSquare();
                break;
            }
            else
            {
                // Save new bitmasks for diagonal values for further use
                diagonalValues1 |= bit1;
                diagonalValues2 |= bit2;
                diagonalValuesHistory[currentRowId][0] = diagonalValues1;
                diagonalValuesHistory[currentRowId][1] = diagonalValues2;

                // Save remaining candidates in row history
                rowsHistoryFlags[currentRowId] = rowCandidates;

                // Mark the row in the array of the used rows
                ClearBit(rowsUsage, gettingRowId);

#ifndef HAS_SIMD
                // Set new row candidates
                rowCandidates = rowsUsage;
#endif

                // Step forward
                currentRowId++;

#ifdef HAS_SIMD
#ifdef __AVX2__
                // load bitmasks for columns which will be on diagonals
                __m256i vCol1 = _mm256_load_si256((const __m256i*)&squareA_MaskT[currentRowId][0]);
                __m256i vCol2 = _mm256_load_si256((const __m256i*)&squareA_MaskT[Rank - 1 - currentRowId][0]);

                // AND loaded values with diagonal masks
                __m256i vDiagMask1 = _mm256_set1_epi16(diagonalValues1);
                __m256i vDiagMask2 = _mm256_set1_epi16(diagonalValues2);

                vCol1 = _mm256_and_si256(vCol1, vDiagMask1);
                vCol2 = _mm256_and_si256(vCol2, vDiagMask2);

                // non-zero means that number is duplicated, zero means that it is unique
                // OR these values together first
                vCol1 = _mm256_or_si256(vCol1, vCol2);

#if defined(__AVX512F__) && defined(__AVX512VL__)
                // check if result is zero and get result as a bitmask
                __mmask16 resultMask = _mm256_testn_epi16_mask(vCol1, vCol1);

                // AND result with masked rowsUsage
                rowCandidates = resultMask & rowsUsage & rowsUsageMasks[currentRowId];
#else  /* !AVX512 */
                // check if result is zero
                vCol1 = _mm256_cmpeq_epi16(vCol1, _mm256_setzero_si256());

                // there are 2 bits per result, so we need to pack int16 to int8 first
                vCol1 = _mm256_packs_epi16(vCol1, _mm256_setzero_si256());
                unsigned int mask = _mm256_movemask_epi8(vCol1);

                // AVX internally has two separate lanes :(
                mask = (mask & 0x000000FF) | ((mask & 0x00FF0000) >> 8);

                // AND result with masked rowsUsage
                rowCandidates = mask & rowsUsage & rowsUsageMasks[currentRowId];
#endif // !AVX512

// AVX causes too much CPU throttling, so any benefits from longer vectors are lost
// and app is slower than SSE ones. Fortunately AVX2 app does not have this problem!
/*#elif defined(__AVX__)
                // AVX has floating-point support only, but still is useful here

                // load bitmasks for columns which will be on diagonals
                __m256 vCol1 = _mm256_load_ps((const float*)&squareA_MaskT[currentRowId][0]);
                __m256 vCol2 = _mm256_load_ps((const float*)&squareA_MaskT[Rank - 1 - currentRowId][0]);

                // AND loaded values with diagonal masks
                __m256 vDiagMask1 = _mm256_castsi256_ps(_mm256_set1_epi16(diagonalValues1));
                __m256 vDiagMask2 = _mm256_castsi256_ps(_mm256_set1_epi16(diagonalValues2));

                vCol1 = _mm256_and_ps(vCol1, vDiagMask1);
                vCol2 = _mm256_and_ps(vCol2, vDiagMask2);

                // non-zero means that number is duplicated, zero means that it is unique
                // OR these values together first
                vCol1 = _mm256_or_ps(vCol1, vCol2);

                // check if result is zero
                __m128i vCol1a = _mm256_castsi256_si128(_mm256_castps_si256(vCol1));
                __m128i vCol1b = _mm_castps_si128(_mm256_extractf128_ps(vCol1, 1));

                vCol1a = _mm_cmpeq_epi16(vCol1a, _mm_setzero_si128());
                vCol1b = _mm_cmpeq_epi16(vCol1b, _mm_setzero_si128());

                // create mask from vector
                // there are 2 bits per result, so we need to pack int16 to int8 first
                vCol1a = _mm_packs_epi16(vCol1a, _mm_setzero_si128());
                vCol1b = _mm_packs_epi16(vCol1b, _mm_setzero_si128());
                unsigned int maska = _mm_movemask_epi8(vCol1a);
                unsigned int maskb = _mm_movemask_epi8(vCol1b);

                // combine masks together, and AND result with masked rowsUsage
                rowCandidates = (maska | (maskb << 8)) & rowsUsage & rowsUsageMasks[currentRowId];
*/
#elif defined(__SSE2__)
                // load bitmasks for columns which will be on diagonals
                __m128i vCol1a = _mm_load_si128((const __m128i*)&squareA_MaskT[currentRowId][0]);
                __m128i vCol1b = _mm_load_si128((const __m128i*)&squareA_MaskT[currentRowId][8]);
                __m128i vCol2a = _mm_load_si128((const __m128i*)&squareA_MaskT[Rank - 1 - currentRowId][0]);
                __m128i vCol2b = _mm_load_si128((const __m128i*)&squareA_MaskT[Rank - 1 - currentRowId][8]);

                // AND loaded values with diagnonal masks
                __m128i vDiagMask1 = _mm_set1_epi16(diagonalValues1);
                __m128i vDiagMask2 = _mm_set1_epi16(diagonalValues2);

                vCol1a = _mm_and_si128(vCol1a, vDiagMask1);
                vCol1b = _mm_and_si128(vCol1b, vDiagMask1);
                vCol2a = _mm_and_si128(vCol2a, vDiagMask2);
                vCol2b = _mm_and_si128(vCol2b, vDiagMask2);

                // non-zero means that number is duplicated, zero means that it is unique
                // OR these values together first
                vCol1a = _mm_or_si128(vCol1a, vCol2a);
                vCol1b = _mm_or_si128(vCol1b, vCol2b);

                // check if result is zero
                vCol1a = _mm_cmpeq_epi16(vCol1a, _mm_setzero_si128());
                vCol1b = _mm_cmpeq_epi16(vCol1b, _mm_setzero_si128());

                // create mask from vector
                // there are 2 bits per result, so we need to pack int16 to int8 first
                vCol1a = _mm_packs_epi16(vCol1a, _mm_setzero_si128());
                vCol1b = _mm_packs_epi16(vCol1b, _mm_setzero_si128());
                unsigned int maska = _mm_movemask_epi8(vCol1a);
                unsigned int maskb = _mm_movemask_epi8(vCol1b);

                // combine masks together, and AND result with masked rowsUsage
                rowCandidates = (maska | (maskb << 8)) & rowsUsage & rowsUsageMasks[currentRowId];

#elif defined(__ARM_NEON)
#ifdef __aarch64__
                // load bitmasks for columns which will be on diagonals
                // for performance reasons load this as a row from transposed square
                uint16x8_t vCol1a = vld1q_u16((const uint16_t*)&squareA_MaskT[currentRowId][0]);
                uint16x8_t vCol1b = vld1q_u16((const uint16_t*)&squareA_MaskT[currentRowId][8]);

                uint16x8_t vCol2a = vld1q_u16((const uint16_t*)&squareA_MaskT[Rank - 1 - currentRowId][0]);
                uint16x8_t vCol2b = vld1q_u16((const uint16_t*)&squareA_MaskT[Rank - 1 - currentRowId][8]);

                // AND loaded values with diagnonal masks
                uint16x8_t vDiagMask1 = vdupq_n_u16(diagonalValues1);
                uint16x8_t vDiagMask2 = vdupq_n_u16(diagonalValues2);

                vCol1a = vandq_u16(vCol1a, vDiagMask1);
                vCol1b = vandq_u16(vCol1b, vDiagMask1);

                vCol2a = vandq_u16(vCol2a, vDiagMask2);
                vCol2b = vandq_u16(vCol2b, vDiagMask2);

                // non-zero means that number is duplicated, zero means that it is unique
                // OR these values together first
                vCol1a = vorrq_u16(vCol1a, vCol2a);
                vCol1b = vorrq_u16(vCol1b, vCol2b);

                // check if result is zero
                vCol1a = vceqq_u16(vCol1a, vdupq_n_u16(0));
                vCol1b = vceqq_u16(vCol1b, vdupq_n_u16(0));

                // create mask from vector
                vCol1a = vandq_u16(vCol1a, vPowersOf2_1);
                vCol1b = vandq_u16(vCol1b, vPowersOf2_2);

                vCol1a = vorrq_u16(vCol1a, vCol1b);

                uint32_t mask = vaddvq_u64(vpaddlq_u32(vpaddlq_u16(vCol1a)));

                // AND result with masked rowsUsage
                rowCandidates = mask & rowsUsage & rowsUsageMasks[currentRowId];
#else /* !__aarch64__ */
                // load bitmasks for columns which will be on diagonals
                // for performance reasons load this as a row from transposed square
                uint16x4_t vCol1a = vld1_u16((const uint16_t*)&squareA_MaskT[currentRowId][0]);
                uint16x4_t vCol1b = vld1_u16((const uint16_t*)&squareA_MaskT[currentRowId][4]);
                uint16x4_t vCol1c = vld1_u16((const uint16_t*)&squareA_MaskT[currentRowId][8]);

                uint16x4_t vCol2a = vld1_u16((const uint16_t*)&squareA_MaskT[Rank - 1 - currentRowId][0]);
                uint16x4_t vCol2b = vld1_u16((const uint16_t*)&squareA_MaskT[Rank - 1 - currentRowId][4]);
                uint16x4_t vCol2c = vld1_u16((const uint16_t*)&squareA_MaskT[Rank - 1 - currentRowId][8]);

                // AND loaded values with diagnonal masks
                uint16x4_t vDiagMask1 = vdup_n_u16(diagonalValues1);
                uint16x4_t vDiagMask2 = vdup_n_u16(diagonalValues2);

                vCol1a = vand_u16(vCol1a, vDiagMask1);
                vCol1b = vand_u16(vCol1b, vDiagMask1);
                vCol1c = vand_u16(vCol1c, vDiagMask1);

                vCol2a = vand_u16(vCol2a, vDiagMask2);
                vCol2b = vand_u16(vCol2b, vDiagMask2);
                vCol2c = vand_u16(vCol2c, vDiagMask2);

                // non-zero means that number is duplicated, zero means that it is unique
                // OR these values together first
                vCol1a = vorr_u16(vCol1a, vCol2a);
                vCol1b = vorr_u16(vCol1b, vCol2b);
                vCol1c = vorr_u16(vCol1c, vCol2c);

                // check if result is zero
                vCol1a = vceq_u16(vCol1a, vdup_n_u16(0));
                vCol1b = vceq_u16(vCol1b, vdup_n_u16(0));
                vCol1c = vceq_u16(vCol1c, vdup_n_u16(0));

                // create mask from vector
                vCol1a = vand_u16(vCol1a, vPowersOf2_1);
                vCol1b = vand_u16(vCol1b, vPowersOf2_2);
                vCol1c = vand_u16(vCol1c, vPowersOf2_3);

                vCol1a = vorr_u16(vCol1a, vCol1b);
                vCol1a = vorr_u16(vCol1a, vCol1c);

                uint64x1_t v = vpaddl_u32(vpaddl_u16(vCol1a));
                uint32_t mask = vget_lane_u32(vreinterpret_u32_u64(v), 0);

                // AND result with masked rowsUsage
                rowCandidates = mask & rowsUsage & rowsUsageMasks[currentRowId];
#endif
#endif // AVX2/SSE2
#endif // HAS_SIMD
                if (!rowCandidates)
                    break;
            }
        }

        // 2nd loop (used to be "else" part) - handle case when there are no row candidates
        while (1)
        {
            // Process not-founding of the new row: step backward, clear the flags of usage,
            // the history of usage, the list of current rows and clear the square itself

            // Step backward
            currentRowId--;
            // Check if we are done
            if (0 == currentRowId)
                return;
            // Get saved values for previous row
            diagonalValues1 = diagonalValuesHistory[currentRowId - 1][0];
            diagonalValues2 = diagonalValuesHistory[currentRowId - 1][1];
            // Clear the flag of row usage
            SetBit(rowsUsage, currentSquareRows[currentRowId]);
            // Get saved candidates
            rowCandidates = rowsHistoryFlags[currentRowId];
            if (rowCandidates)
                break;
        }
    }
}
