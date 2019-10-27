#pragma once

#include "Vector.h"
#include "Matrix.h"

using namespace mat_vec;

Vector operator*(double k, const Vector &v)
{
    Vector V(v.size());
    for(size_t i = 0; i < v.size(); i++)
    {
        V[i] = k * v.operator[](i);
    }
    return V;
}

Vector::Vector(size_t size, double value):vector(new double[size]), length(size)
{
    for(size_t i = 0; i < size; i++)
    {
        vector[i] = value;
    }
}

Vector::~Vector()
{
    delete[] this->vector;
}

Vector::Vector(const Vector &src):length(src.length), vector(new double[src.length])
{
    for(size_t i = 0; i < length; i++)
    {
        this->vector[i] = src.vector[i];
    }
}

size_t Vector::size() const
{
    return length;
}

double &Vector::operator[](size_t n)
{
    return vector[n];
}

double Vector::operator[](size_t n) const
{
    return this->vector[n];
}

double Vector::norm() const
{
    double inv_lenth;
    for(size_t i = 0; i < length; i++)
    {
        inv_lenth += vector[i] * vector[i];
    }

    inv_lenth = sqrt(inv_lenth);
    return inv_lenth;
}

void Vector::normalize()
{
    double inv_lenth;
    for(size_t i = 0; i < length; i++)
    {
        inv_lenth += vector[i] * vector[i];
    }
    inv_lenth = sqrt(inv_lenth);


    for(size_t i = 0; i < length; i++)
    {
        vector[i] = vector[i]/inv_lenth;
    }
}

Vector Vector::normalized() const
{
    Vector v(length);
    double inv_lenth = 0;
    for(size_t i = 0; i < length; i++)
    {
        inv_lenth += vector[i] * vector[i];
    }
    inv_lenth = sqrt(inv_lenth);

    for(size_t i = 0; i < this->length; i++)
    {
        v[i] = (this->vector[i])/inv_lenth;
    }

    return v;
}


Vector Vector::operator+(const Vector &rhs) const
{
    Vector v(rhs.size());
    for(size_t i = 0; i < v.size(); i++)
    {
        v[i] = rhs.operator[](i);
    }
    if(length != v.size())
    {
        std::cout << "Операция невозможна!" << std::endl;
        return v;
    }else
    {
        for(size_t i = 0; i < v.size(); i++)
        {
            v[i] += vector[i];
        }
        return v;
    }
}

Vector &Vector::operator+=(const Vector &rhs)
{
    Vector v(rhs.size());
    if(length != v.size())
    {
        std::cout << "Операция невозможна!" << std::endl;
        return v;
    }else
    {
        for(size_t i = 0; i < length; i++)
        {
            v[i] = this->vector[i] + rhs.operator[](i);
        }
    }
    return v;
}

Vector Vector::operator-(const Vector &rhs) const
{
    Vector v(rhs.size());
    for(size_t i = 0; i < v.size(); i++)
    {
        v[i] = rhs.operator[](i);
    }
    if(length != v.size())
    {
        std::cout << "Операция невозможна!" << std::endl;
        return v;
    }else
    {
        for(size_t i = 0; i < v.size(); i++)
        {
            v[i] -= vector[i];
        }
        return v;
    }
}

Vector &Vector::operator-=(const Vector &rhs)
{
    Vector v(rhs.length);
    if(v.length != rhs.length)
    {
        std::cout << "Операция невозможна!" << std::endl;
        return v;
    }else
    {
        for(size_t i = 0; i < rhs.length; i++)
        {
            v.vector[i] -= rhs.vector[i];
        }
    }
    return v;
}

Vector Vector::operator^(const Vector &rhs) const
{
    Vector v(rhs.length);
    for(size_t i = 0; i < v.length; i++)
    {
        v[i] = rhs.vector[i];
    }
    if(length != v.length)
    {
        std::cout << "Операция невозможна!" << std::endl;
        return v;
    }else
    {
        for(size_t i = 0; i < v.size(); i++)
        {
            v[i] *= vector[i];
        }
        return v;
    }
}

Vector &Vector::operator^=(const Vector &rhs)
{
    if(length != rhs.length)
    {
        std::cout << "Операция невозможна!" << std::endl;
        return *this;
    }else
    {
        for(size_t i = 0; i < length; i++)
        {
            this->vector[i] = this->vector[i] * rhs.operator[](i);
        }
    }
    return *this;
}

double Vector::operator*(const Vector &rhs) const
{
    double cos = 0;
    for(size_t i = 0; i < length; i++)
    {
        cos += vector[i] * rhs.operator[](i);
    }

    return cos;
}

Vector Vector::operator*(double k) const
{
    Vector v(length);
    for(size_t i = 0; i < length; i++)
    {
        v[i] = vector[i] * k;
    }
    return v;
}

Vector &Vector::operator*=(double k)
{
    Vector v(length);
    for(size_t i = 0; i < length; i++)
    {
        v.vector[i] = this->vector[i] * k;
    }
    return v;
}

Vector &Vector::operator/=(double k)
{
    for(size_t i = 0; i < length; i++)
    {
        this->vector[i] = (this->vector[i]) / k;
    }
    return *this;
}

Vector Vector::operator/(double k) const
{
    Vector v(length);
    for(size_t i = 0; i < length; i++)
    {
        v.vector[i] = vector[i] / k;
    }
    return v;
}

bool Vector::operator==(const Vector &rhs) const
{
    if(length != rhs.size())
    {
        return false;
    }
    for(size_t i = 0; i < length; i++)
    {
        if(vector[i] != rhs.operator[](i))
        {
            return false;
        }
    }
    return true;
}

bool Vector::operator!=(const Vector &rhs) const
{
    size_t cnt = 0;
    if(length == rhs.length)
    {
        for(size_t i = 0; i < length; i++)
        {
            if(vector[i] == rhs.vector[i])
            {
                cnt++;
            }
        }
        if(cnt == length)
        {
            return false;
        }else
        {
            return true;
        }
    }else
    {
        return true;
    }
    
}


Vector &Vector::operator=(const Vector &rhs)
{
    delete []vector;
    vector = new double[rhs.length];
    length = rhs.length;
    for(size_t i = 0; i < length; i++)
    {
        vector[i] = rhs.vector[i];
    }

    return *this;
}

Vector &Vector::operator*=(const Matrix &mat)
{
    double *v = new double[mat.shape().second];
    for(size_t i = 0; i < length; i++)
    {
        for(size_t j = 0; j < mat.shape().second; j++)
        {
            v[i] += vector[i] * mat.get(i, j);
        }
    }
    delete []vector;
    vector = new double[mat.shape().second];
    for(size_t i = 0; i < mat.shape().second; i++)
    {
        vector[i] = v[i];
    }
    delete []v;
    return *this;
}

Vector Vector::operator*(const Matrix &mat) const
{
    Vector v(mat.shape().second, 0);
    for(size_t i = 0; i < length; i++)
    {
        for(size_t j = 0; j < mat.shape().second; j++)
        {
            v.vector[i] += vector[i] * mat.get(i, j);
        }
    }
    return v;
}
