#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

#include <iomanip>

#include <functional>

#include <ciso646>

#include <algorithm>

#include <vector>

#include <fstream>

#include <cstdlib>

template <class Type>
class Matrix
{
	template <typename Type>
	friend std::istream& operator >> (std::istream&, Matrix<Type>&);
	template <typename Type>
	friend std::ostream& operator << (std::ostream&, const Matrix<Type>&);
public:
	Matrix() = default;
	Matrix(const char*);
	Matrix(const Matrix<Type>&);
	Matrix(std::initializer_list<Type>);
	Matrix(std::initializer_list<std::initializer_list<Type>>);
	Matrix(const size_t& m, const size_t& n = ULLONG_MAX) { helperConstructor(m, n); }
	~Matrix() { helperDestructor(); }

	size_t sizeRow() const noexcept { return row; }
	size_t sizeCol() const noexcept { return col; }
	bool isNull() const noexcept;
	bool isUnit() const noexcept;
	bool isEmpty() const noexcept { return matrixPtr == nullptr; }
	bool isSquare() const noexcept { return row == col; }
	bool isDiagonal() const noexcept;
	bool isSymmetric() const noexcept;
	void exportInFile(const char*) const;
	void swap(Matrix<Type>&) noexcept;
	void clear() noexcept { helperDestructor(); }
	void transpose() noexcept;
	void fill(const Type&) noexcept;
	void deleteRow(const size_t&);
	void deleteCol(const size_t&);
	void swapRows(const size_t&, const size_t&);
	void swapCols(const size_t&, const size_t&);
	Type summa() const noexcept;
	Type maxElement() const noexcept;
	Type minElement() const noexcept;
	std::vector<Type> getRow(const size_t&) const ;
	std::vector<Type> getCol(const size_t&) const;

	const Matrix<Type>& operator =  (const Matrix<Type>&) noexcept;
	const Matrix<Type>& operator += (const Matrix<Type>&);
	const Matrix<Type>& operator -= (const Matrix<Type>&);
	const Matrix<Type>& operator *= (const Matrix<Type>&);
	const Matrix<Type>& operator += (const Type&);
	const Matrix<Type>& operator -= (const Type&);
	const Matrix<Type>& operator *= (const Type&);
	Matrix<Type> operator + (const Matrix<Type>&) const;
	Matrix<Type> operator - (const Matrix<Type>&) const;
	Matrix<Type> operator * (const Matrix<Type>&) const;
	Matrix<Type> operator + (const Type&) const;
	Matrix<Type> operator - (const Type&) const;
	Matrix<Type> operator * (const Type&) const;
	Matrix<Type>& operator ++ ();
	Matrix<Type> operator  ++ (int);
	Matrix<Type>& operator -- ();
	Matrix<Type> operator  -- (int);
	bool operator == (const Matrix<Type>&) const noexcept;
	bool operator != (const Matrix<Type>& other) const noexcept { return (*this == other); }
	Type& operator () (const size_t&, const size_t&);
	const Type& operator () (const size_t&, const size_t&) const;
private:
	Type** matrixPtr;
	size_t row, 
		   col;

	void helperConstructor(const size_t&, const size_t&);
	void helperDestructor();
	void helperCrement(const bool&);
	bool examinationOfFile(const char*);
};

template <class Type>
inline Matrix<Type>::Matrix(const char* fileName)
{
	std::ifstream objectFile(fileName, std::ios::in);
	if (!objectFile) 
	{
		std::cerr << "Файл " << fileName << " не может быть открыт\n";
		exit(EXIT_FAILURE);
	}

	objectFile.close();

	if (!examinationOfFile(fileName)) 
	{
		std::cerr << "Содержимое файла некорректно\n";
		exit(EXIT_FAILURE);
	}

	objectFile.open(fileName, std::ios::in);

	objectFile >> row >> col;

	helperConstructor(row, col);

	for (size_t i = 0; i < row; ++i)
	{
		for (size_t j = 0; j < col; ++j)
		{
			objectFile >> matrixPtr[i][j];
		}
	}

	objectFile.close();
}

template <class Type>
inline Matrix<Type>::Matrix(const Matrix<Type>& other)
{
	helperConstructor(other.row, other.col);

	for (size_t i = 0; i < this->row; ++i)
	{
		for (size_t j = 0; j < this->col; ++j)
		{
			this->matrixPtr[i][j] = other.matrixPtr[i][j];
		}
	}
}

template<class Type>
inline Matrix<Type>::Matrix(std::initializer_list<Type> list)
{
	helperConstructor(1, list.size());

	for (size_t i = 0; i < list.size(); ++i)
	{
		matrixPtr[0][i] = *(list.begin() + i);
	}
}

template<class Type>
inline Matrix<Type>::Matrix(std::initializer_list<std::initializer_list<Type>> list)
{
	size_t maxSizeOfColumn = 0;

	for (const std::initializer_list<Type>& subList : list)
	{
		maxSizeOfColumn = std::max(maxSizeOfColumn, subList.size());
	}

	helperConstructor(list.size(), maxSizeOfColumn);

	size_t i = 0;

	for (const std::initializer_list<Type>& subList : list)
	{
		size_t j = 0;

		for (const Type& value : subList)
		{
			matrixPtr[i][j] = value;
			++j;
		}
		while (j < maxSizeOfColumn)
		{
			matrixPtr[i][j] = 0;
			++j;
		}
		++i;
	}
}

template <class Type>
inline bool Matrix<Type>::isNull() const noexcept
{
	if (isEmpty())
	{
		return false;
	}

	for (size_t i = 0; i < row; ++i)
	{
		for (size_t j = 0; j < col; ++j)
		{
			if (i == j)
			{
				if (matrixPtr[i][j] != 1) return false;
			}
			else
			{
				if (matrixPtr[i][j] != 0) return false;
			}
		}
	}

	return true;
}

template <class Type>
inline bool Matrix<Type>::isUnit() const noexcept
{
	if (isEmpty() or row != col)
	{
		return false;
	}

	for (size_t i = 0; i < row; ++i)
	{
		for (size_t j = 0; j < col; ++j)
		{
			if (i == j)
			{
				if (matrixPtr[i][j] != 1)
				{
					return false;
				}
			}
			else
			{
				if (matrixPtr[i][j] != 0)
				{
					return false;
				}
			}
		}
	}

	return true;
}

template <class Type>
inline bool Matrix<Type>::isDiagonal() const noexcept
{
	if (isEmpty())
	{
		return false;
	}

	for (size_t i = 0; i < row; ++i)
	{
		for (size_t j = 0; j < col; ++j)
		{
			if (i != j && matrixPtr[i][j] != 0)
			{
				return false;
			}
		}
	}

	return true;
}

template <class Type>
inline bool Matrix<Type>::isSymmetric() const noexcept
{
	Matrix<Type> temp(*this);

	return temp == temp.transpose();
}

template <class Type>
inline void Matrix<Type>::exportInFile(const char* fileName) const
{
	std::ofstream output(fileName, std::ios::out);
	if (!output)
	{
		throw std::invalid_argument("Не удалось открыть файл для записи!");
	}

	output << this->row << ' ' << this->col << ' ';

	for (size_t i = 0; i < row; ++i)
	{
		for (size_t j = 0; j < col; ++j)
		{
			output << this->matrixPtr[i][j] << ' ';
		}
	}

	output.close();
}

template <class Type>
inline void Matrix<Type>::swap(Matrix<Type>& other) noexcept
{
	std::swap(matrixPtr, other.matrixPtr);
	std::swap(row, other.row);
	std::swap(col, other.col);
}

template <class Type>
inline void Matrix<Type>::transpose() noexcept
{
	Matrix<Type> temp(col, row);

	for (size_t i = 0; i < row; ++i)
	{
		for (size_t j = 0; j < col; ++j)
		{
			temp.matrixPtr[j][i] = matrixPtr[i][j];
		}
	}

	*this = temp;
}

template <class Type>
inline void Matrix<Type>::fill(const Type& value) noexcept
{
	for (size_t i = 0; i < row; ++i)
	{
		for (size_t j = 0; j < col; ++j)
		{
			matrixPtr[i][j] = value;
		}
	}
}

template <class Type>
inline void Matrix<Type>::deleteRow(const size_t& numOfRow)
{
	if (numOfRow < 1 || numOfRow > this->row)
	{
		throw std::length_error("Попытка удалить несуществующую строку!");
	}

	if (row == 1)
	{
		helperDestructor();
		return;
	}

	size_t indexOfNumOfRow = 0;

	Matrix<Type> newMatrix(row - 1, col);

	for (size_t i = 0; i < newMatrix.row; ++i, ++indexOfNumOfRow)
	{
		if (i + 1 == numOfRow)
		{
			++indexOfNumOfRow;
		}

		for (size_t j = 0; j < newMatrix.col; ++j)
		{
			newMatrix.matrixPtr[i][j] = matrixPtr[indexOfNumOfRow][j];
		}
	}

	*this = newMatrix;
}

template <class Type>
inline void Matrix<Type>::deleteCol(const size_t& numOfCol)
{
	if (numOfCol < 1 || numOfCol > this->col)
	{
		throw std::length_error("Попытка удалить несуществующий столбец!");
	}

	if (col == 1)
	{
		helperDestructor();
		return;
	}

	Matrix<Type> newMatrix(row, col - 1);

	for (size_t i = 0; i < newMatrix.row; ++i)
	{
		for (size_t j = 0, indexOfNumOfCol = 0; j < newMatrix.col; ++j, ++indexOfNumOfCol)
		{
			if (j + 1 == numOfCol)
			{
				++indexOfNumOfCol;
			}

			newMatrix.matrixPtr[i][j] = matrixPtr[i][indexOfNumOfCol];
		}
	}

	*this = newMatrix;
}

template <class Type>
inline void Matrix<Type>::swapRows(const size_t& firstIndex, const size_t& secondIndex)
{
	if (firstIndex < 1 or secondIndex < 1 or firstIndex > row or secondIndex > row)
	{
		throw std::out_of_range("Обращение к несуществующей строке!");
	}

	if (firstIndex == secondIndex)
	{
		return;
	}

	for (size_t i = 0; i < col; ++i)
	{
		std::swap(matrixPtr[firstIndex - 1][i], matrixPtr[secondIndex - 1][i]);
	}
}

template <class Type>
inline void Matrix<Type>::swapCols(const size_t& firstIndex, const size_t& secondIndex)
{
	if (firstIndex < 1 or secondIndex < 1 or firstIndex > col or secondIndex > col)
	{
		throw std::out_of_range("Обращение к несуществующему столбцу!");
	}

	if (firstIndex == secondIndex)
	{
		return;
	}

	for (size_t i = 0; i < row; ++i)
	{
		std::swap(matrixPtr[i][firstIndex - 1], matrixPtr[i][secondIndex - 1]);
	}
}

template <class Type>
inline Type Matrix<Type>::summa() const noexcept
{
	Type sum = 0;

	for (size_t i = 0; i < row; ++i)
	{
		for (size_t j = 0; j < col; ++j)
		{
			sum += matrixPtr[i][j];
		}
	}

	return sum;
}

template <class Type>
inline Type Matrix<Type>::maxElement() const noexcept
{
	Type maximum = matrixPtr[0][0];

	for (size_t i = 0; i < row; ++i)
	{
		for (size_t j = 0; j < col; ++j)
		{
			if (matrixPtr[i][j] > maximum)
			{
				maximum = matrixPtr[i][j];
			}
		}
	}

	return maximum;
}

template <class Type>
inline Type Matrix<Type>::minElement() const noexcept
{
	Type minimum = matrixPtr[0][0];

	for (size_t i = 0; i < row; ++i)
	{
		for (size_t j = 0; j < col; ++j)
		{
			if (matrixPtr[i][j] < minimum)
			{
				minimum = matrixPtr[i][j];
			}
		}
	}

	return minimum;
}

template <class Type>
inline std::vector<Type> Matrix<Type>::getRow(const size_t& numOfRow) const
{
	if (numOfRow < 1 || numOfRow > this->row)
	{
		throw std::length_error("Попытка получить несуществующую строку!");
	}

	if (isEmpty())
	{
		return std::vector<Type>();
	}

	std::vector<Type> ofRow(col);

	for (size_t index = 0; index < col; ++index)
	{
		ofRow[index] = matrixPtr[numOfRow - 1][index];
	}

	return ofRow;
}

template <class Type>
inline std::vector<Type> Matrix<Type>::getCol(const size_t& numOfCol) const
{
	if (numOfCol < 1 || numOfCol > this->col)
	{
		throw std::length_error("Попытка получить несуществующую строку!");
	}

	if (isEmpty())
	{
		return std::vector<Type>();
	}

	std::vector<Type> ofCol(row);

	for (size_t index = 0; index < row; ++index)
	{
		ofCol[index] = matrixPtr[index][numOfCol - 1];
	}

	return ofCol;
}

template <class Type>
inline const Matrix<Type>& Matrix<Type>::operator =  (const Matrix<Type>& other) noexcept
{
	if (row != other.row or col != other.col)
	{
		helperDestructor();
		helperConstructor(other.row, other.col);
	}

	for (size_t i = 0; i < row; ++i)
	{
		for (size_t j = 0; j < col; ++j)
		{
			matrixPtr[i][j] = other.matrixPtr[i][j];
		}
	}

	return *this;
}

template <class Type>
inline const Matrix<Type>& Matrix<Type>::operator += (const Matrix<Type>& other)
{
	*this = *this + other;
	
	return *this;
}

template <class Type>
inline const Matrix<Type>& Matrix<Type>::operator -= (const Matrix<Type>& other)
{
	*this = *this - other;

	return *this;
}

template <class Type>
inline const Matrix<Type>& Matrix<Type>::operator *= (const Matrix<Type>& other)
{
	*this = *this * other;
	
	return *this;
}

template <class Type>
inline const Matrix<Type>& Matrix<Type>::operator += (const Type& value)
{
	*this = *this + value;
	
	return *this;
}

template <class Type>
inline const Matrix<Type>& Matrix<Type>::operator -= (const Type& value)
{
	*this = *this - value;
	
	return *this;
}

template <class Type>
inline const Matrix<Type>& Matrix<Type>::operator *= (const Type& value)
{
	*this = *this * value;
	
	return *this;
}

template <class Type>
inline Matrix<Type> Matrix<Type>::operator + (const Matrix<Type>& other) const
{
	if (row != other.row or col != other.col)
	{
		throw std::logic_error("Матрицы имеют не одинаковые размеры для их сложения!");
	}

	Matrix<Type> result(*this);

	for (size_t i = 0; i < this->row; ++i)
	{
		for (size_t j = 0; j < this->col; ++j)
		{
			result.matrixPtr[i][j] += other.matrixPtr[i][j];
		}
	}

	return result;
}

template <class Type>
inline Matrix<Type> Matrix<Type>::operator - (const Matrix<Type>& other) const
{
	if (row != other.row or col != other.col)
	{
		throw std::logic_error("Матрицы имеют не одинаковые размеры для их вычитания!");
	}

	Matrix<Type> result(*this);

	for (size_t i = 0; i < this->row; ++i)
	{
		for (size_t j = 0; j < this->col; ++j)
		{
			result.matrixPtr[i][j] -= other.matrixPtr[i][j];
		}
	}

	return result;
}

template <class Type>
inline Matrix<Type> Matrix<Type>::operator * (const Matrix<Type>& other) const
{
	if (col != other.row)
	{
		throw std::logic_error("Матрицы не соответствуют требованиям для их умножения!");
	}

	Matrix<Type> result(this->row, other.col);

	for (size_t i = 0; i < this->row; ++i)
	{
		for (size_t j = 0; j < other.col; ++j)
		{
			result.matrixPtr[i][j] = 0;

			for (size_t k = 0; k < this->col; ++k)
			{
				result.matrixPtr[i][j] += this->matrixPtr[i][k] * other.matrixPtr[k][j];
			}
		}
	}

	return result;
}

template <class Type>
inline Matrix<Type> Matrix<Type>::operator + (const Type& value) const
{
	Matrix<Type> result(this->row, this->col);

	for (size_t i = 0; i < this->row; ++i)
	{
		for (size_t j = 0; j < this->col; ++j)
		{
			result.matrixPtr[i][j] = this->matrixPtr[i][j] + value;
		}
	}

	return result;
}

template <class Type>
inline Matrix<Type> Matrix<Type>::operator - (const Type& value) const
{
	Matrix<Type> result(this->row, this->col);

	for (size_t i = 0; i < this->row; ++i)
	{
		for (size_t j = 0; j < this->col; ++j)
		{
			result.matrixPtr[i][j] = this->matrixPtr[i][j] - value;
		}
	}

	return result;
}

template <class Type>
inline Matrix<Type> Matrix<Type>::operator * (const Type& value) const
{
	Matrix<Type> result(this->row, this->col);

	for (size_t i = 0; i < this->row; ++i)
	{
		for (size_t j = 0; j < this->col; ++j)
		{
			result.matrixPtr[i][j] = this->matrixPtr[i][j] * value;
		}
	}

	return result;
}

template <class Type>
inline Matrix<Type>& Matrix<Type>::operator ++ ()
{
	helperCrement(true);

	return *this;
}

template <class Type>
inline Matrix<Type> Matrix<Type>::operator ++ (int)
{
	Matrix<Type> temp(*this);

	helperCrement(true);

	return temp;
}

template <class Type>
inline Matrix<Type>& Matrix<Type>::operator -- ()
{
	helperCrement(false);

	return *this;
}

template <class Type>
inline Matrix<Type> Matrix<Type>::operator -- (int)
{
	Matrix<Type> temp(*this);

	helperCrement(false);

	return temp;
}

template <class Type>
inline bool Matrix<Type>::operator == (const Matrix<Type>& other) const noexcept
{
	if (row != other.row or col != other.col)
	{
		return false;
	}

	for (size_t i = 0; i < row; ++i)
	{
		for (size_t j = 0; j < col; ++j)
		{
			if (matrixPtr[i][j] != other.matrixPtr[i][j])
			{
				return false;
			}
		}
	}

	return true;
}

template <class Type>
inline Type& Matrix<Type>::operator () (const size_t& a, const size_t& b)
{
	if (a < 0 or a >= row or b < 0 or b >= col)
	{
		throw std::out_of_range("Индекс выходит за пределы матрицы!");
	}

	return matrixPtr[a][b];
}

template <class Type>
inline const Type& Matrix<Type>::operator () (const size_t& a, const size_t& b) const
{
	if (a < 0 or a >= row or b < 0 or b >= col)
	{
		throw std::out_of_range("Индекс выходит за пределы матрицы!");
	}

	return matrixPtr[a][b];
}

template <class Type>
inline void Matrix<Type>::helperConstructor(const size_t& m, const size_t& n)
{
	if (m < 0 or n < 0)
	{
		throw std::logic_error("Попытка создать матрицу отрицательного размера!");
	}

	row = m;
	col = (n == ULLONG_MAX ? m : n);

	matrixPtr = new Type * [row];

	for (size_t i = 0; i < row; ++i)
	{
		matrixPtr[i] = new Type[col];
	}
}

template<class Type>
inline void Matrix<Type>::helperDestructor()
{
	for (size_t i = 0; i < row; ++i)
	{
		delete[] matrixPtr[i];
	}

	delete[] matrixPtr;

	row = col = 0;
}

template <class Type>
inline void Matrix<Type>::helperCrement(const bool& is)
{
	std::function<void(Type&, const bool&)> help = [](Type& value, const bool& is) -> void
		{
			is ? ++value : --value;
		};

	for (size_t i = 0; i < row; ++i)
	{
		for (size_t j = 0; j < col; ++j)
		{
			help(matrixPtr[i][j], is);
		}
	}
}

template <class Type>
bool Matrix<Type>::examinationOfFile(const char* file)
{

	std::ifstream examFile(file, std::ios::in);

	int examRow, examCol;

	examFile >> examRow >> examCol;

	if (examRow <= 0 or examCol <= 0) 
	{
		// Некорректные размеры матрицы
		return false;
	}

	if (examFile.fail()) 
	{
		// Ошибка при чтении размеров матрицы
		return false;
	}

	for (size_t i = 0; i < examRow; ++i) 
	{
		for (size_t j = 0; j < examCol; ++j) 
		{
			Type value;

			examFile >> value;

			if (examFile.fail())
			{
				// Ошибка при чтении элемента матрицы
				return false;
			}
		}
	}

	char nextChar;

	if (examFile >> nextChar) 
	{
		// Символы после матрицы
		return false;
	}

	examFile.close();

	return true;
}

template <typename MatrixType>
inline std::istream& operator >> (std::istream& input, Matrix<MatrixType>& matrix)
{
	if (matrix.isEmpty())
	{
		throw std::invalid_argument("Попытка инициализировать несуществующую область памяти!");
	}

	std::cout << "Введите значения матрицы размера : " << matrix.row <<
		" на " << matrix.col << " --->\n\n";

	for (size_t i = 0; i < matrix.row; ++i)
	{
		std::cout << i + 1 << " | ";

		for (size_t j = 0; j < matrix.col; ++j)
		{
			input >> matrix.matrixPtr[i][j];
		}
	}

	return input;
}

template<typename MatrixType>
inline std::ostream& operator << (std::ostream& output, const Matrix<MatrixType>& matrix)
{
	if (matrix.isEmpty())
	{
		output << "Пустая матрица!";
	}
	else
	{
		output << std::setiosflags(std::ios::fixed | std::ios::showpoint) << std::setprecision(1);

		size_t sizeOfColumnMax = matrix.col;

		long double* columnMax = new long double[sizeOfColumnMax];

		for (size_t i = 0; i < sizeOfColumnMax; ++i)
		{
			columnMax[i] = 0.0;
		}

		for (size_t i = 0; i < sizeOfColumnMax; ++i)
		{
			long double maximum = (long double)LLONG_MIN;

			for (size_t k = 0; k < matrix.row; ++k)
			{
				if (matrix.matrixPtr[k][i] > maximum)
				{
					maximum = (long double)matrix.matrixPtr[k][i];
				}
			}

			columnMax[i] = maximum;
		}

		output << "     ";

		std::function<size_t(size_t)> longs = [](size_t number) -> size_t
			{
				if (number == 0)
				{
					return 1;
				}

				size_t result = 0;

				while (number > 0)
				{
					++result;

					number /= 10;
				}

				return result;
			};

		for (size_t i = 0; i < matrix.col; ++i)
		{
			output << '[' << i + 1 << ']' << "  ";

			size_t space = longs((size_t)columnMax[i]);

			for (size_t j = 0; j < space - 1; ++j)
			{
				output << ' ';
			}
		}

		for (size_t l = 0; l < sizeOfColumnMax; ++l)
		{
			columnMax[l] = longs((size_t)columnMax[l]);
		}

		output << '\n';

		for (size_t i = 0; i < matrix.row; ++i)
		{
			output << ((i + 1) > 9 ? "[" : " [") << i + 1 << "]  ";

			for (size_t j = 0; j < matrix.col; ++j)
			{
				output << matrix.matrixPtr[i][j] << "  ";

				size_t plus;

				if (matrix.matrixPtr[i][j] != 0)
				{
					plus = (size_t)columnMax[j] - longs((size_t)matrix.matrixPtr[i][j]);
				}
				else
				{
					plus = longs((size_t)columnMax[j]) - 1;

				}
				for (size_t k = 0; k < plus; ++k)
				{
					output << ' ';
				}
			}
			output << '\n';
		}

		delete[] columnMax;
	}

	return output;
}

#endif // !MATRIX_H
