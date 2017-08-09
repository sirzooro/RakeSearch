// Диагональный латинский квадрат

# if !defined Square_h
# define Square_h

# include <iostream>

using namespace std;

class Square
{
public:
  static const int Rank = 9;          // Ранг квадрата
  static const int Empty = -1;        // Пустое, не заданное значение

  static int OrthoDegree(Square a, Square b);  // Степень ортогональности квадратов a и b

  Square();                       // Конструктор по умолчанию
  Square(int source[Rank][Rank]); // Конструктор создания квадрата по матрице
  Square(Square& source);         // Конструктор копирования

  int operator == (Square& value);    // Перегрузка оператора сравнения - сравниваются компоненты матрицы
  Square& operator = (Square& value); // Перегрузка оператора присвоения
  friend std::ostream& operator << (std::ostream& os, Square& value);    // Перегрузка оператора вывода данных квадрата
  friend std::istream& operator >> (std::istream& is, Square& value);    // Перегрузка оператора считывания данных квадрата

  int IsDiagonal(); // Проверка квадрата на то, что он является диагональным латинским квадратом
  int IsLatin();    // Проверка квадрата на то, что он является латинским квадратом
  void Initialize(int source[Rank][Rank]);    // Инициализация компонентов квадрата
  void Reset();                    // Сброс всех значеий перемененных
  void Read(std::istream& is);     // Чтение квадрата из потока
  void Write(std::ostream& os);    // Запись квадрата в поток

  int Matrix[Rank][Rank];          // Матрица квадрата

protected:
private:
};

# endif
