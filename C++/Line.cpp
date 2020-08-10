#include "Line.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

/// Normal vector of the line
Vector Line::normal_vector() const
{
    return Vector(a, b);
}

/// Guiding vector of the line
Vector Line::guiding_vector() const
{
    return Vector(b, -a);
}

/// Normalizes the line
Line Line::normalize() const
{
    double l = normal_vector().length();
    return Line(a / l, b / l, c / l);
}

istream &operator >>(istream &in, Line &l)
{
    in >> l.a >> l.b >> l.c;
    return in;
}

ostream &operator <<(ostream &out, const Line l)
{
    out << l.a << ' ' << l.b << ' ' << l.c;
    return out;
}

Line Line::operator+(const double& k) const
{
    return Line(a, b, c + k);
}

Line& Line::operator+=(const double& k)
{
    this->c += k;
    return *this;
}

Line Line::operator-(const double& k) const
{
    return Line(a, b, c - k);
}

Line& Line::operator-=(const double& k)
{
    this->c -= k;
    return *this;
}

bool Line::operator==(const Line& l) const
{
    Line l1 = this->normalize(), l2 = l.normalize();
    return ((fabs(l1.a - l2.a) < eps) && (fabs(l1.b - l2.b) < eps) && (fabs(l1.c - l2.c) < eps))
            || ((fabs(l1.a + l2.a) < eps) && (fabs(l1.b + l2.b) < eps) && (fabs(l1.c + l2.c) < eps));
}

bool Line::operator!=(const Line& l) const
{
    Line l1 = this->normalize(), l2 = l.normalize();
    return ((fabs(l1.a - l2.a) > eps) || (fabs(l1.b - l2.b) > eps) || (fabs(l1.c - l2.c) > eps))
            && ((fabs(l1.a + l2.a) > eps) || (fabs(l1.b + l2.b) > eps) || (fabs(l1.c + l2.c) > eps));
}

bool Point::belongs(Line l) const
{
    return fabs(l.a * x + l.b * y + l.c) < eps;
}

/// True, if lines are parallel
bool parallel(Line l1, Line l2)
{
    return (fabs(l1.a - l2.a) < eps) && (fabs(l1.b - l2.b) < eps);
}

/// Distance from the point to the line
double Point::distance(Line l) const
{
    return fabs(l.a * x + l.b * y + l.c) / l.normal_vector().length();
}

/// Sign distance from the point to the line
double Point::sign_distance(Line l) const
{
    return (l.a * x + l.b * y + l.c) / l.normal_vector().length();
}

/// Projection of the point on the line
Point Point::projection(Line l) const
{
    Vector n = (-l.normal_vector().normalize()) * sign_distance(l);
        return (*this + n);
}

/// Finds intersection point P of the two lines L1 and L2.
/// Returns -1, if lines are equal; 0, if lines are parallel; 1, if lines intersect.
int intersect(Line l1, Line l2, Point &p)
{
    if (l1 == l2) return -1;
    if (parallel(l1, l2)) return 0;

    p.x = (l2.b * l1.c - l1.b * l2.c) / (l2.a * l1.b - l1.a * l2.b);
    p.y = (l1.b != 0 ? (- l1.c - l1.a * p.x) / l1.b : (- l2.c - l2.a * p.x) / l2.b);
    return 1;
}
