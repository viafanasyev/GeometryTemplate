#include "Vector.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

/// Length of the vector
double Vector::length() const
{
    return hypot(x, y);
}

/// Polar angle of the vector
double Vector::polar_angle() const
{
    double angle = atan2(y, x);
    if (angle < 0)
    {
        angle += 2 * M_PI;
    }
    return angle;
}

/// Normalizes the vector (divides by it's length)
Vector Vector::normalize() const
{
    if ((fabs(x) < eps) && (fabs(y) < eps)) // If Vector = {0, 0}
        return *this;

    double l = length();
    return Vector(x / l, y / l);
}

/// Rotates the vector on angle around the point p
Vector Vector::rotate(double angle, Point p) const
{
    Vector tmp;
    tmp.x = (x - p.x) * cos(angle) - (y - p.y) * sin(angle) + p.x;
    tmp.y = (y - p.y) * cos(angle) + (x - p.x) * sin(angle) + p.y;
    return tmp;
}

istream &operator >>(istream &in, Vector &v)
{
    in >> v.x >> v.y;
    return in;
}

ostream &operator <<(ostream &out, const Vector v)
{
    out << v.x << ' ' << v.y;
    return out;
}

Vector Vector::operator-() const
{
    return Vector(-x, -y);
}

Vector Vector::operator+(const Vector& v) const
{
    return Vector(x + v.x, y + v.y);
}

Vector& Vector::operator+=(const Vector& v)
{
    this->x += v.x;
    this->y += v.y;
    return *this;
}

Vector Vector::operator-(const Vector& v) const
{
    return Vector(x - v.x, y - v.y);
}

Vector& Vector::operator-=(const Vector& v)
{
    this->x -= v.x;
    this->y -= v.y;
    return *this;
}

Vector Vector::operator*(const double& k) const
{
    return Vector(x * k, y * k);
}

Vector& Vector::operator*=(const double& k)
{
    this->x *= k;
    this->y *= k;
    return *this;
}

Vector Vector::operator/(const double& k) const
{
    if (fabs(k) < eps)
        throw overflow_error("Division by zero exception");
    return Vector(x / k, y / k);
}

Vector& Vector::operator/=(const double& k)
{
    if (fabs(k) < eps)
        throw overflow_error("Division by zero exception");
    this->x /= k;
    this->y /= k;
    return *this;
}

bool Vector::operator==(const Vector& v) const
{
    return (fabs(x - v.x) < eps) && (fabs(y - v.y) < eps);
}

bool Vector::operator!=(const Vector& v) const
{
    return (fabs(x - v.x) > eps) || (fabs(y - v.y) > eps);
}

/// Scalar product of two vectors
double dot_product(Vector v1, Vector v2)
{
    return v1.x * v2.x + v1.y * v2.y;
}

/// Vector product of two vectors
double cross_product(Vector v1, Vector v2)
{
    return v1.x * v2.y - v1.y * v2.x;
}

/// Angle between two vectors
double vectors_angle(Vector v1, Vector v2)
{
    if ((v1 == Vector(0, 0)) || (v2 == Vector(0, 0)))
        return M_PI / 2.0;
    return acos(dot_product(v1, v2) / (v1.length() * v2.length()));
}

Point Point::operator+(const Vector& v) const
{
    return Point(x + v.x, y + v.y);
}

Point Point::operator-(const Vector& v) const
{
    return Point(x - v.x, y - v.y);
}

Point& Point::operator+=(const Vector& v)
{
    x += v.x;
    y += v.y;
    return *this;
}

Point& Point::operator-=(const Vector& v)
{
    x -= v.x;
    y -= v.y;
    return *this;
}
