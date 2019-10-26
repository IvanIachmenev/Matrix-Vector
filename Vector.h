#pragma once

#include "Base.h"

namespace mat_vec {

// Умножение всех элементов вектора на число слева (k * v) 1
Vector operator*(double k, const Vector &v);

class Vector {
public:
  // Конструирует вектор размера size со значениями value 2
  explicit Vector(size_t size, double value = 0);
  
  // Конструктор копирования 3
  Vector(const Vector &src);
  
  // Оператор присваивания 4
  Vector &operator=(const Vector &rhs);
  
  // Деструктор 5
  ~Vector();

  // Возвращает размер вектора 6
  size_t size() const;
  
  // Доступ к n-му элементу вектора 7
  double operator[](size_t n) const;
  double &operator[](size_t n);

  // L2 норма вектора 8
  double norm() const;
  
  // Возвращает новый вектор, полученный нормализацией текущего (this) 9
  Vector normalized() const;
  
  // Нормализует текущий вектор 10
  void normalize();

  // Поэлементное сложение векторов 11
  Vector operator+(const Vector &rhs) const;
  Vector &operator+=(const Vector &rhs);

  // Поэлементное вычитание векторов 12
  Vector operator-(const Vector &rhs) const;
  Vector &operator-=(const Vector &rhs);

  // Поэлементное умножение векторов 13
  Vector operator^(const Vector &rhs) const;
  Vector &operator^=(const Vector &rhs);

  // Скалярное произведение 14
  double operator*(const Vector &rhs) const;

  // Умножение всех элементов вектора на скаляр справа (v * k) 15
  Vector operator*(double k) const;
  Vector &operator*=(double k);

  // Деление всех элементов вектора на скаляр 16
  Vector operator/(double k) const;
  Vector &operator/=(double k);

  // Умножение вектора на матрицу 17
  Vector operator*(const Matrix &mat) const;
  Vector &operator*=(const Matrix &mat);

  // Поэлементное сравнение18
  bool operator==(const Vector &rhs) const;
  bool operator!=(const Vector &rhs) const;


private:
  double *vector, **matrix;
  size_t length, wigth;
};
} // namespace mat_vec