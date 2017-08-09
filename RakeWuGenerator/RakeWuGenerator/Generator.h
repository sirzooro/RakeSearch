// Генератор Диагональных латинских квадратов

# if !defined Generator_h
# define Generator_h

# include <iostream>
# include <fstream>
# include "Square.h"

using namespace std;

class Generator
{
public:
  Generator();                  // Конструктор по умолчанию
  Generator(Generator& source); // Конструктор копировния
  void GenerateWorkunits(string start, string result, string checkpoint, string temp, int workunits);  // Выполнение генерации заданий

  Generator& operator = (Generator&  value);                  // Оператор копирования
  friend std::ostream& operator << (std::ostream& os, Generator& value);    // Оператор записи состояние генератора
  friend std::istream& operator >> (std::istream& is, Generator& value);    // Оператор считывания состояния генератора

private:
  static const int Rank = Square::Rank;  // Ранг квадрата (для удобства)
  static const int Free = 1;        // Флаг свободного для использования значения
  static const int Used = 0;        // Флаг задействования значения в какой-либо клетке диагонали | строки | столбца
  static const int CellsInPath = 16;    // Число обрабатываемых клеток. Для диагоналей 9 ранга + 1 клетка в первой строке - 16 клеток
  //static const int CellsInPath = 7;    // Число обрабатываемых клеток. Для главной диагонали 8 ранга - 7 клеток
  static const int Yes = 1;       // Флаг "Да"
  static const int No = 0;        // Флаг "Нет"

  void Initialize(string start, string result, string checkpoint, string temp, int workunits);  // Инициализация поиска
  void Start();              // Функция генерации заданий обходом клеток
  void Reset();              // Сброс всех значений внутренних структур

  void CopyState(Generator& source); // Копирование состояния
  void Read(std::istream& is);       // Считывание состояние генератора из потока
  void Write(std::ostream& os);      // Запись состояния генератора в поток

  Square newSquare;          // Генерируемый квадрат

  int path[CellsInPath][2];  // Путь заполнения матрицы квадрата - path[i][0] - строка на шаге i, path[i][1] - столбец
  int keyRowId;              // Идентификатор строки ключевой клетки - по значению которой расчёт будет останавливаться
  int keyColumnId;           // Идентификатор столбца ключевой клетки
  int keyValue;              // Значение ключевой клетки, по достижению которого расчёт будет останавливаться

  int primary[Rank];         // Содержимое главной диагонали
  int secondary[Rank];       // Содержимое побочной диагонали
  int columns[Rank][Rank];   // Матрица значений, использовавшихся в столбцах - columns[значение][столбец] = 0|1. 0 - значение занято. 1 - свободно.
  int rows[Rank][Rank];          // Матрица значений, использовавшихся в строках - rows[строка][значение] = 0|1
  int cellsHistory[Rank][Rank][Rank];    // Куб значений, которые использовались для формирования построенной части квадрата - cellsHistory[строка][столбец][значение]

  string startParametersFileName; // Название файла с параметрами запуска расчёта
  string resultFileName;          // Название файла с результатами
  string checkpointFileName;      // Название файла контрольной точки
  string tempCheckpointFileName;  // Временное название файла новой контрольной точки

  int isInitialized;        // Флаг успешной инициализации поиска
  int rowId;                // Идентификатор строки обратываемой клетки
  int columnId;             // Идентификатор столбца обрабатываемой клетки
  int cellId;               // Идентификатор клетки в перечне шагов обхода квадрата
  
  int snapshotNumber;       // Номер снимка по которым генерируются задания (от i-1 до i)
  int workunitKeyRowId;     // Строка ключевой клетки, прописываемая в задании
  int workunitKeyColumnId;  // Столбец ключевой клетки, приписываемый в задании

  int workunitsCount;       // Число заданий, созданных в рамках этого прогона (не входит в контрольную точку!)
  int workunitsToGenerate;  // Число заданий, которые необходимо сгенерировать (не входит в контрольную точку!)

  void CreateCheckpoint();  // Создание контрольной точки
  void ProcessWorkunit();   // Формирование очередного задания
  void ProcessWorkunitR8(); // Формирование очередного задания для квадратов 8-го ранга
  void ProcessWorkunitR9(); // Формирование очередного задания для квадратов 9-го ранга

  string generatorStateHeader;  // Заголовок, после которого в файле параметров или контрольной точки идёт состояние генератора диагональных квадратов
};

# endif
