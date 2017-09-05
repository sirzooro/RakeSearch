// Поиск пар диагональных латинских квадратов методом перетасовки строк

# if !defined MovePairSearch_h
# define MovePairSearch_h

# include <iostream>
# include <string>

# include "Generator.h"

using namespace std;

class MovePairSearch
{
public:
  static const int Rank = Square::Rank;

  MovePairSearch();              // Конструктор по умолчанию
  void Reset();                  // Сброс настроек поиска
  void ClearBeforeNextSearch();  // Очистка необходимых переменных перед очередным поиском
  void InitializeMoveSearch(string start, string result, string checkpoint, string temp);  // Инициализация поиска
  void StartMoveSearch();        // Запуск поиска ортогональных квадратов методом перестановки строк
  void OnSquareGenerated(Square newSquare);  // Обработчик события построения ДЛК, запускающий поиск к нему пары

private:                              
  static const int CheckpointInterval = 1000000;  // Интервал создания контрольных точек
  static const int OrhoSquaresCacheSize = 32;     // Размер кэша для хранения квадратов, ортогональных обрабатываемому

  void MoveRows();                  // Перетасовка строк заданного ДЛК в поиске ОДЛК к нему
  void ProcessOrthoSquare();        // Обработка найденного ортогонального квадрата
  void CheckMutualOrthogonality();  // Проверка взаимной ортогональности квадратов
  void CreateCheckpoint();          // Создание контрольной точки
  void Read(std::istream& is);      // Чтение состояния поиска из потока
  void Write(std::ostream& os);     // Запись состояния поиска в поток
  void ShowSearchTotals();          // Отображение общих итогов поиска

  Generator squareAGenerator;       // Генератор ДЛК
  int squareA[Rank][Rank];          // Исходный ДЛК, строки в котором будут переставляться
  int squareB[Rank][Rank];          // Подбираемый ДЛК, строки внутри которого и переставляются
  int rowsUsage[Rank];              // Массив флагов возможно задействования строк на данный момент - rowsUsage[номер строки] = 0 | 1, где 0 - строка уже используется, 1 - нет.
  int rowsHistory[Rank][Rank];          // Массив истории возможного использования строк - rowsHistory[номер строки][значение] = 0 | 1, 
                          // где 0 - строка с номером "значение" использовалась для строки "номер строки" формируемого квадрата. 1 - строку ещё можно использовать
  int currentSquareRows[Rank];          // Массив с перечнем текущих строк, использованных в квадрате. На i-й позиции - номер используемой строки
  int pairsCount;                  // Число обнаруженных диагональных квадратов в перестановках строк из найдеенного squareA
  int totalPairsCount;              // Общее число обнаруженных диагональных квадратов - в рамках всего поиска
  int totalSquaresWithPairs;            // Общее число квадратов, к которым найден хотя бы один ортогональный
  int totalProcessedSquaresLarge;          // Число обработанных ДЛК, поступивших от генератора - в миллиардах
  int totalProcessedSquaresSmall;          // Число обработанных ДЛК, поступивших от генератора в пределах до миллиарда

  Square orthoSquares[OrhoSquaresCacheSize];    // Кэш для хранения квадратов, ортогональных обрабатываемому

  string startParametersFileName;   // Название файла с параметрами запуска расчёта
  string resultFileName;            // Название файла с результатами
  string checkpointFileName;        // Название файла контрольной точки
  string tempCheckpointFileName;    // Временное название файла новой контрольной точки

  int isInitialized;                // Флаг инициализированности поиска
  int isStartFromCheckpoint;        // Флаг запуска с контрольной точки

  string moveSearchGlobalHeader;    // Заголовок, после которого в файле идут данные о состоянии поиска методом перетасовки строк
  string moveSearchComponentHeader; // Заголовок, после которого в файле идут данные о состоянии компоненты перетасовки строк
  static const bool isDebug = false; // Флаг вывода отладочной информации
};

# endif
