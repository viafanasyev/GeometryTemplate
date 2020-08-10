#include "Point.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

/// Polar angle of the point
double Point::polar_angle() const
{
    double angle = atan2(y, x);
    if (angle < 0)
    {
        angle += 2 * M_PI;
    }
    return angle;
}

/// Rotates the point on angle around the point p
Point Point::rotate(double angle, Point p) const
{
    Point tmp;
    tmp.x = (x - p.x) * cos(angle) - (y - p.y) * sin(angle) + p.x;
    tmp.y = (y - p.y) * cos(angle) + (x - p.x) * sin(angle) + p.y;
    return tmp;
}


/// Distance between two points
double Point::distance(Point p) const
{
    return hypot(fabs(p.x - x), fabs(p.y - y));
}

istream &operator >>(istream &in, Point &p)
{
    in >> p.x >> p.y;
    return in;
}

ostream &operator <<(ostream &out, const Point p)
{
    out << p.x << ' ' << p.y;
    return out;
}

Point Point::operator+(const Point& p) const
{
    return Point(x + p.x, y + p.y);
}

Point& Point::operator+=(const Point& p)
{
    this->x += p.x;
    this->y += p.y;
    return *this;
}

Point Point::operator-(const Point& p) const
{
    return Point(x - p.x, y - p.y);
}

Point& Point::operator-=(const Point& p)
{
    this->x -= p.x;
    this->y -= p.y;
    return *this;
}

Point Point::operator*(const double& k) const
{
    return Point(x * k, y * k);
}

Point& Point::operator*=(const double& k)
{
    this->x *= k;
    this->y *= k;
    return *this;
}

Point Point::operator/(const double& k) const
{
    if (fabs(k) < eps)
        throw overflow_error("Division by zero exception");
    return Point(x / k, y / k);
}

Point& Point::operator/=(const double& k)
{
    if (fabs(k) < eps)
        throw overflow_error("Division by zero exception");
    this->x /= k;
    this->y /= k;
    return *this;
}

bool Point::operator==(const Point& p) const
{
    return (fabs(x - p.x) < eps) && (fabs(y - p.y) < eps);
}

bool Point::operator!=(const Point& p) const
{
    return (fabs(x - p.x) > eps) || (fabs(y - p.y) > eps);
}
