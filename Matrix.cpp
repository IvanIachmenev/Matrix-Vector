#include "Matrix.h"

using namespace mat_vec;

Matrix::Matrix(size_t size, double value):wigth(size), higth(size), matrix(new double*[size])
{
    for(size_t i = 0; i < size; i++)
    {
        matrix[i] = new double[size];
        for(size_t j = 0; j < size; j++)
        {
            matrix[i][j] = value;
        }
    }
}

Matrix Matrix::eye(size_t size)
{
    Matrix m(size, 0.0);
    for(size_t i = 0; i < size; i++)
    {
        for(size_t j = 0; j < size; j++)
        {
            if(i == j)
            {
                m.matrix[i][j] = 1;
            }
        }
    }
    return m;
}

Matrix::Matrix(size_t rows, size_t cols, double value):higth(rows), wigth(cols), matrix(new double*[rows])
{
    for(size_t i = 0; i < rows; i++)
    {
        matrix[i] = new double[cols];
        for(size_t j = 0; j < cols; j++)
        {
            matrix[i][j] = value;
        }
    }
}

Matrix::Matrix(const Matrix &src)
{
    std::pair<size_t, size_t> hw = src.shape();
    higth = hw.first;
    wigth = hw.second;
    matrix = new double*[higth];
    for(size_t i = 0; i < higth; i++)
    {
        matrix[i] = new double[wigth];
        for(size_t j = 0; j < wigth; j++)
        {
            matrix[i][j] = src.get(i, j);
        }
    }
}

Matrix::~Matrix()
{
    for(size_t i = 0; i < higth; i++)
    {
        delete []matrix[i];
    }
    delete []matrix;
}

double Matrix::get(size_t row, size_t col) const
{
    return matrix[row][col];
}

double Matrix::get(size_t row, size_t col, double k) const
{
    return matrix[row][col] = k;
}

std::pair<size_t, size_t> Matrix::shape() const
{
    std::pair<size_t, size_t> hw;
    hw.first = higth;
    hw.second = wigth;
    return hw;
}

Matrix &Matrix::operator+=(const Matrix &rhs)
{
    std::pair<size_t, size_t> hw = rhs.shape();
    if(higth != hw.first && wigth != hw.second)
    {
        Matrix k(0, 0, 0);
        std::cout << "Опеация невозможна!" << std::endl;
        return k;
    }

    Matrix m(higth, wigth, 0);
    for(size_t i = 0; i < higth; i++)
    {
        for(size_t j = 0; j < wigth; j++)
        {
            matrix[i][j] += rhs.get(i, j);
            m.matrix[i][j] = matrix[i][j];
        }
    }
    return m;
}

Matrix Matrix::operator+(const Matrix &rhs) const
{
    std::pair<size_t, size_t> hw = rhs.shape();
    if(higth != hw.first && wigth != hw.second)
    {
        Matrix k(0, 0, 0);
        std::cout << "Опеация невозможна!" << std::endl;
        return k;
    }

    Matrix m(higth, wigth, 0);
    for(size_t i = 0; i < higth; i++)
    {
        for(size_t j = 0; j < wigth; j++)
        {
            matrix[i][j] += rhs.get(i, j);
            m.matrix[i][j] = matrix[i][j];
        }
    }
    return m;
}

Matrix Matrix::operator-(const Matrix &rhs) const
{
    std::pair<size_t, size_t> hw = rhs.shape();
    if(higth != hw.first && wigth != hw.second)
    {
        Matrix k(0, 0, 0);
        std::cout << "Опеация невозможна!" << std::endl;
        return k;
    }

    Matrix m(higth, wigth, 0);
    for(size_t i = 0; i < higth; i++)
    {
        for(size_t j = 0; j < wigth; j++)
        {
            matrix[i][j] -= rhs.get(i, j);
            m.matrix[i][j] = matrix[i][j];
        }
    }
    return m;
}

Matrix &Matrix::operator-=(const Matrix &rhs)
{
    std::pair<size_t, size_t> hw = rhs.shape();
    if(higth != hw.first && wigth != hw.second)
    {
        Matrix k(0, 0, 0);
        std::cout << "Опеация невозможна!" << std::endl;
        return k;
    }

    Matrix m(higth, wigth, 0);
    for(size_t i = 0; i < higth; i++)
    {
        for(size_t j = 0; j < wigth; j++)
        {
            matrix[i][j] -= rhs.get(i, j);
            m.matrix[i][j] = matrix[i][j];
        }
    }
    return m;
}

Matrix Matrix::operator*(double k) const
{
    Matrix m(higth, wigth, 0);
    for(size_t i = 0; i < higth; i++)
    {
        for(size_t j = 0; j < wigth; j++)
        {
            matrix[i][j] *= k;
            m.matrix[i][j] = matrix[i][j];
        }
    }
    return m;
}

Matrix &Matrix::operator*=(double k)
{
    Matrix m(higth, wigth, 0);
    for(size_t i = 0; i < higth; i++)
    {
        for(size_t j = 0; j < wigth; j++)
        {
            matrix[i][j] *= k;
            m.matrix[i][j] = matrix[i][j];
        }
    }
    return m;
}

Matrix Matrix::operator/(double k) const
{
    Matrix m(higth, wigth, 0);
    for(size_t i = 0; i < higth; i++)
    {
        for(size_t j = 0; j < wigth; j++)
        {
            matrix[i][j] /= k;
            m.matrix[i][j] = matrix[i][j];
        }
    }
    return m;
}

Matrix &Matrix::operator/=(double k)
{
    Matrix m(higth, wigth, 0);
    for(size_t i = 0; i < higth; i++)
    {
        for(size_t j = 0; j < wigth; j++)
        {
            matrix[i][j] /= k;
            m.matrix[i][j] = matrix[i][j];
        }
    }
    return m;
}

Matrix Matrix::operator*(const Matrix &rhs) const
{
    std::pair<size_t, size_t> hw = rhs.shape();
    if(higth != hw.second)
    {
        Matrix k(0, 0, 0);
        std::cout << "Операция невозможна!" << std::endl;
        return k;
    }

    Matrix m(hw.second, higth, 0);
    for(size_t i = 0; i < higth; i++)
    {
        for(size_t j = 0; j < hw.second; j++)
        {
            for(size_t k = 0; k < wigth; k++)
            {
                m.matrix[i][j] += matrix[i][k] * rhs.matrix[k][j];
            }
        }
    }
    return m;
}

Matrix &Matrix::operator*=(const Matrix &rhs)
{
    std::pair<size_t, size_t> hw = rhs.shape();
    if(higth != hw.second)
    {
        Matrix k(0, 0, 0);
        std::cout << "Операция невозможна!" << std::endl;
        return k;
    }

    Matrix m(hw.second, higth, 0);
    for(size_t i = 0; i < higth; i++)
    {
        for(size_t j = 0; j < hw.second; j++)
        {
            for(size_t k = 0; k < wigth; k++)
            {
                m.matrix[i][j] += matrix[i][k] * rhs.matrix[k][j];
            }
        }
    }
    return m;
}

bool Matrix::operator==(const Matrix &rhs) const
{
    std::pair<size_t, size_t> hw = rhs.shape();
    if(higth != hw.first || wigth != hw.second)
    {
        return false;
    }

    for(size_t i = 0; i < higth; i++)
    {
        for(size_t j = 0; j < wigth; j++)
        {
            if(matrix[i][j] != rhs.matrix[i][j])
            {
                return false;
            }
        }
    }

    return true;
}

bool Matrix::operator!=(const Matrix &rhs) const
{
    std::pair<size_t, size_t> hw = rhs.shape();
    size_t cnt = 0;
    if(higth == hw.first && wigth == hw.second)
    {
        for(size_t i = 0; i < higth; i++)
        {
            for(size_t j = 0; j < wigth; j++)
            {
                if(matrix[i][j] == rhs.matrix[i][j])
                {
                    cnt++;
                }
            }
        }
        if(cnt != hw.first * hw.second)
        {
            return true;
        }else
        {
            return false;
        }
    }


    return true;
}

Vector Matrix::operator*(const Vector &vec) const
{
    size_t length = vec.size();
    if(wigth != length)
    {
        Vector k(0, 0);
        std::cout << "Операция невозможна!" << std::endl;
        return k;
    }

    Vector v(higth, 0);
    double k = 0;
    for(size_t i = 0; i < higth; i++)
    {
        for(size_t j = 0; j < length; j++)
        {
            k += matrix[i][j] * vec.operator[](j);
        }
        v.operator[](i) = k;
        k = 0;
    }
    return v;
}

double Matrix::det() const
{
    if(higth != wigth)
    {
        std::cout << "Операция невозможна!" << std::endl;
        return 0;
    }

    double **matrix = new double*[higth];
    for(size_t i = 0; i < higth; i++)
    {
        matrix[i] = new double[higth];
        for(size_t j = 0; j < higth; j++)
        {
            matrix[i][j] = this->matrix[i][j];
        }
    } 

    const double EPS = 1E-9;
    double det = 1;
    for (size_t i = 0; i < higth; ++i)
    {
        int k = i;
        for (size_t j = i+1; j < higth; ++j)
        {
            if (abs(matrix[j][i]) > abs(matrix[k][i]))
            {
                k = j;
            }
        }
	    if (abs(matrix[k][i]) < EPS)
        {
		    det = 0;
		    break;
	    }
        for(size_t j = 0; j < higth; j++)
        {
            double a = matrix[i][j];
            matrix[i][j] = matrix[k][j];
            matrix[k][j] = a;
        }
        if (i != k)
        {
            det = -det;
        }
        det *= matrix[i][i];
        for (int j = i+1; j < higth; ++j)
        {
            matrix[i][j] /= matrix[i][i];
        }
        for (int j = 0; j < higth; ++j)
        {
            if (j != i && abs (matrix[j][i]) > EPS)
            {
                for (int k = i+1; k < higth; ++k)
                {
                    matrix[j][k] -= matrix[i][k] * matrix[j][i];
                }
            }
        }
    }

    for(size_t i = 0; i < higth; i++)
    {
        delete []matrix[i];
    }
    delete []matrix;

    return det;
}

Matrix Matrix::inv() const
{
    Matrix m(higth, wigth, 0.0);
    for(size_t i = 0; i < higth; i++)
    {
        for(size_t j = 0; j < higth; j++)
        {
            m.matrix[i][j] = matrix[i][j];
        }
    }
    double det = m.det();
    if(det == 0)
    {
        Matrix k(0, 0, 0);
        std::cout << "Операция невозможна!" << std::endl;
        return k;
    }
    
    for (size_t i = 0; i< higth; i++)
	{
		if (matrix[i][i] == 0)
			for (size_t j = i+1; j< higth; j++)
			{
				if (matrix[j][i]==1)
				{
					for (size_t k =0; k< 2*higth; k++)
					{
						double c = matrix[j][k];
						matrix[j][k] = matrix[i][k];
						matrix[i][k] = c;
					}
					break;
				}
			}
		for (size_t k = i+1; k < higth; k++)
		{
			if (matrix[k][i] == 1)
			{
				for (size_t j = 0; j < 2*higth; j++)
				{
					matrix[k][j] *= matrix[i][j];
				}
			}
		}
	}
    for (size_t i = higth-1; i >= 0; i--)
    {
        for (size_t k = i-1; k >= 0; k--)
        {
            if (matrix[k][i] == 1)
            {
                for (size_t j = 0; j < 2*higth; j++)
                {
                    matrix[k][j] *= matrix[i][j];
                }
            }
        }
    }
	
	//копируем в обратную матрицу В
	for (size_t i = 0; i < higth; i++)
    {
		for (size_t j = 0, k = higth; j < higth; j++, k++)
        {
            m.matrix[i][j] = matrix[i][k];
        }
    }
	return m;
}

Matrix& Matrix::operator=(const mat_vec::Matrix &rhs)
{
    for(size_t i = 0; i < higth; i++)
    {
        delete []matrix[i];
    }
    delete[] this->matrix;
    this->higth = rhs.higth;
    this->wigth = rhs.wigth;

    matrix = new double*[higth];
    for(size_t i = 0; i < higth; i++)
    {
        matrix[i] = new double[wigth];
    }

    for(size_t i = 0; i < higth; i++)
    {
        for(size_t j = 0; j < wigth; j++)
        {
            matrix[i][j] = rhs.matrix[i][j];
        }
    }
    return *this;
}

void Matrix::get(size_t row, size_t col, double k)
{
    matrix[row][col] = k;
}

void Matrix::reshape(size_t rows, size_t cols)
{
    double *m = new double[rows*cols];
    size_t k = 0;
    for(size_t i = 0; i < higth; i++)
    {
        for(size_t j = 0; j < wigth; j++)
        {
            m[k] = matrix[i][j];
            k++;
        }
    }
    for(size_t i = 0; i < higth; i++)
    {
        delete []matrix[i];
    }
    delete []matrix;

    higth = rows;
    wigth = cols;
    k = 0;

    matrix = new double*[higth];
    for(size_t i = 0; i < higth; i++)
    {
        matrix[i] = new double[wigth];
        for(size_t j = 0; j < wigth; j++)
        {
            matrix[i][j] = m[k];
            k++;
        }
    }
}

Matrix Matrix::transposed() const
{
    Matrix m(wigth, higth);
    for(size_t i = 0; i < wigth; i++)
    {
        m.matrix[i] = new double[higth];
        for(size_t j = 0; j < higth; j++)
        {
            m.matrix[i][j] = matrix[j][i];
        }
    }
    
    return m;
}

void Matrix::transpose()
{
    double **m = new double*[wigth]; 
    for(size_t i = 0; i < wigth; i++)
    {
        m[i] = new double[higth];
        for(size_t j = 0; j < higth; j++)
        {
            m[i][j] = matrix[j][i];
        }
    }
    
    for(size_t i = 0; i < higth; i++)
    {
        delete []matrix[i];
    }
    delete []matrix;
    
    double a = higth;
    higth = wigth;
    wigth = a;
    matrix = new double*[higth];
    for(size_t i = 0; i < higth; i++)
    {
        matrix[i] = new double[wigth];
        for(size_t j = 0; j < wigth; j++)
        {
            matrix[i][j] = m[i][j];
        }
    }
    for(size_t i = 0; i < higth; i++)
    {
        delete []m[i];
    }
    delete []m;
}
