// Диагональный латинский квадрат

# include "Square.h"

using namespace std;

// Конструктор по умолчанию. "Зануление" всех значений
Square::Square()
{
	Reset();
}


// Создание квадрата по заданной матрице
Square::Square(int source[Rank][Rank])
{
	Initialize(source);
}


// Конструктор копирования
Square::Square(Square& source)
{
	Initialize(source.Matrix);
}


// Инициализация внутренних структур
void Square::Initialize(int source[Rank][Rank])
{
	for (int rowId = 0; rowId < Rank; rowId++)
	{
		for (int columnId = 0; columnId < Rank; columnId++)
		{
			Matrix[rowId][columnId] = source[rowId][columnId];
		}
	}
}


// Сброс значений внутренних переменных
void Square::Reset()
{
	for (int rowId = 0; rowId < Rank; rowId++)
	{
		for (int columnId = 0; columnId < Rank; columnId++)
		{
			Matrix[rowId][columnId] = Empty;
		}
	}
}


// Оператор сравнения
int Square::operator == (Square& value)
{
	int isEqual = 1;

	for (int rowId = 0; rowId < Rank && isEqual; rowId++)
	{
		for (int columnId = 0; columnId < Rank && isEqual; columnId++)
		{
			if (Matrix[rowId][columnId] != value.Matrix[rowId][columnId])
			{
				isEqual = 0;
			}
		}
	}

return isEqual;
}


// Оператор присвоения
Square& Square::operator = (Square& value)
{
	Initialize(value.Matrix);

return *this;
}


// Оператор вывода данных квадрата
std::ostream& operator << (std::ostream& os, Square& value)
{
	value.Write(os);

return os;
}


// Оператор считывания данных квадрата
std::istream& operator >> (std::istream& is, Square& value)
{
	value.Read(is);

return is;
}


// Чтение квадрата из потока
void Square::Read(std::istream& is)
{
	char readedChar;	// "Буфер", используемый для посимволного чтения

	// Считывание символа начала блока информации
	do
	{
		is >> readedChar;
	}
	while (readedChar != HeadToken);

	// Считывание матрицы квадрата
	for (int rowId = 0; rowId < Rank; rowId++)
	{
		for (int columnId = 0; columnId < Rank; columnId++)
		{
			is >> Matrix[rowId][columnId];
		}
	}

	// Считывание символа окончания блока информации
	do
	{
		is >> readedChar;
	}
	while (readedChar != TailToken);
}


// Запись квадрата в поток
void Square::Write(std::ostream& os)
{
	// Запись символа начала блока информации
	os << HeadToken << endl;

	// Запись матрицы квадрата
	for (int rowId = 0; rowId < Rank; rowId++)
	{
		for (int columnId = 0; columnId < Rank; columnId++)
		{
			os << Matrix[rowId][columnId] << " ";
		}
		os << endl;
	}

	// Запись символа окончания блока информации
	os << TailToken << endl;
}

// Проверка квадрата на то, что он является диагональным латинским квадратом
int Square::IsDiagonal()
{
	int isDiagonal = 1;

	// Проверка первой диагонали - диагонали [0;0] - [rank;rank]
	for (int itemId = 0; itemId < Rank && isDiagonal; itemId++)
	{
		// Проверка на совпадение элемента [itemId; itemId] со всеми остальными элементами диагонали
		for (int comparedId = itemId + 1; comparedId < Rank && isDiagonal; comparedId++)
		{
			if (Matrix[itemId][itemId] == Matrix[comparedId][comparedId])
			{
				isDiagonal = 0;
			}
		}
	}

	// Проверка второй диагонали - диагонали [rank - itemId - 1; itemId]
	for (int itemId = 0; itemId < Rank && isDiagonal; itemId ++)
	{
		// Проверка на совпадение элемента [rank - itemId - 1; itemId] со всеми остальными элементами диагонали
		for (int comparedId = itemId + 1; comparedId < Rank && isDiagonal; comparedId++)
		{
			if (Matrix[(Rank - itemId - 1)][itemId] == Matrix[(Rank - comparedId - 1)][comparedId])
			{
				isDiagonal = 0;
			}
		}
	}

	return isDiagonal;
}


// Проверка квадрата на то, что он является латинским квадратом
int Square::IsLatin()
{
	int isLatin = 1;

	// Проверка столбцов квадрата
	for (int columnId = 0; columnId < Rank && isLatin; columnId++)
	{
		// Проверка корректности столбца columnId
		for (int rowId = 0; rowId < Rank && isLatin; rowId++)
		{
			// Проверка на совпадение элемента [rowId; columnId] со всеми остальными элементами столбца columndId
			for (int comparedRowId = rowId + 1; comparedRowId < Rank && isLatin; comparedRowId++)
			{
				if (Matrix[comparedRowId][columnId] == Matrix[rowId][columnId])
				{
					isLatin = 0;
				}
			}
		}
	}

	// Проверка строк квадрата
	for (int rowId = 0; rowId < Rank && isLatin; rowId++)
	{
		// Проверка корректности строки rowId
		for (int columnId = 0; columnId < Rank && isLatin; columnId++)
		{
			// Проверка на совпадение элемента [rowId; columnId] со всеми остальными элементами строки rowId
			for (int comparedColumnId = columnId + 1; comparedColumnId < Rank && isLatin; comparedColumnId++)
			{
				if (Matrix[rowId][columnId] == Matrix[rowId][comparedColumnId])
				{
					isLatin = 0;
				}
			}
		}
	}

	return isLatin;
}


// Проверка ортогональности квадратов a и b
int Square::OrthoDegree(Square a, Square b)
{
	int degree = 0;				// Степерь ортогональности
	int freePair[Rank][Rank];	// Массив использования пар значений в получающемся греко-латинском квадрате

	// Инциализируем все пары как свободные
	for (int i = 0; i < Rank; i++)
	{
		for (int j = 0; j < Rank; j++)
		{
			freePair[i][j] = 1;
		}
	}

	// Отмечаем пары, использованные в греко-латинском квадрате
	for (int rowId = 0; rowId < Rank; rowId++)
	{
		for (int columnId = 0; columnId < Rank; columnId++)
		{
			if (freePair[a.Matrix[rowId][columnId]][b.Matrix[rowId][columnId]])
			{
				freePair[a.Matrix[rowId][columnId]][b.Matrix[rowId][columnId]] = 0;
				degree++;
			}
		}
	}

return degree;
}
