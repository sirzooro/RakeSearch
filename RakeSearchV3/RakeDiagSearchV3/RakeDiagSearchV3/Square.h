// Диагональный латинский квадрат

#if !defined Square_h
#define Square_h

#include <iostream>
#include "Helpers.h"

using namespace std;

class Square
{
public:
    static const int Rank = 10;        // Ранг квадрата
    static const int Empty = -1;       // Пустое, не заданное значение
    static const char HeadToken = '{'; // Символ начала информации о квадрате в потоке
    static const char TailToken = '}'; // Символ окончания информации о квадрате в потоке

    static int OrthoDegree(const Square &a, const Square &b); // Степень ортогональности квадратов a и b

    Square();                       // Конструктор по умолчанию
    Square(int source[Rank][Rank]); // Конструктор создания квадрата по матрице
    Square(const Square &source);         // Конструктор копирования

    int operator==(const Square &value) const; // Перегрузка оператора сравнения - сравниваются компоненты матрицы
    Square &operator=(const Square &value); // Перегрузка оператора присвоения
    friend std::ostream &operator<<(std::ostream &os, const Square &value); // Перегрузка оператора вывода данных квадрата
    friend std::istream &operator>>(std::istream &is, Square &value); // Перегрузка оператора считывания данных квадрата

    int IsDiagonal() const; // Проверка квадрата на то, что он является диагональным латинским квадратом
    int IsLatin() const; // Проверка квадрата на то, что он является латинским квадратом
    void Initialize(const int source[Rank][Rank]); // Инициализация компонентов квадрата
    void Reset();                            // Сброс всех значеий перемененных
    void Read(std::istream &is);             // Чтение квадрата из потока
    void Write(std::ostream &os) const;            // Запись квадрата в поток

    int Matrix[Rank][Rank] ALIGNED; // Матрица квадрата

protected:
private:
};

#endif
