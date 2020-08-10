#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include "Point.h"

/// Scalar product of two vectors
double dot_product(Vector v1, Vector v2);

/// Vector product of two vectors
double cross_product(Vector v1, Vector v2);

/// Angle between two vectors
double vectors_angle(Vector v1, Vector v2);

class Vector
{
public:
    double x, y;

    Vector(double x, double y)
    {
        this->x = x;
        this->y = y;
    }

    Vector(Point p) : Vector(p.x, p.y) { }

    Vector(double x_begin, double y_begin, double x_end, double y_end) : Vector(x_end - x_begin, y_end - y_begin) { }

    Vector(Point begin, Point end) : Vector(end.x - begin.x, end.y - begin.y) { }

    Vector() : Vector(0, 0) { }

    /// Length of the vector
    double length() const;

    /// Polar angle of the vector
    double polar_angle() const;

    /// Normalizes the vector (divides by it's length)
    Vector normalize() const;

    /// Rotates the vector on angle around the point p
    Vector rotate(double angle, Point p = Point(0, 0)) const;

    friend std::istream &operator >>(std::istream &in, Vector &v);

    friend std::ostream &operator <<(std::ostream &out, const Vector v);

    Vector operator-() const;

    Vector operator+(const Vector& v) const;

    Vector& operator+=(const Vector& v);

    Vector operator-(const Vector& v) const;

    Vector& operator-=(const Vector& v);

    Vector operator*(const double& k) const;

    Vector& operator*=(const double& k);

    Vector operator/(const double& k) const;

    Vector& operator/=(const double& k);

    bool operator==(const Vector& v) const;

    bool operator!=(const Vector& v) const;
};

#endif // VECTOR_H
