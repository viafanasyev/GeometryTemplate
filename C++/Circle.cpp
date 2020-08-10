#include "Circle.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

istream &operator >>(istream &in, Circle &c)
{
    in >> c.center >> c.radius;
    return in;
}

ostream &operator <<(ostream &out, const Circle c)
{
    out << c.center << ' ' << c.radius;
    return out;
}

Circle Circle::operator+(const Vector& v) const
{
    return Circle(center + v, radius);
}

Circle& Circle::operator+=(const Vector& v)
{
    this->center += v;
    return *this;
}

Circle Circle::operator-(const Vector& v) const
{
    return Circle(center - v, radius);
}

Circle& Circle::operator-=(const Vector& v)
{
    this->center -= v;
    return *this;
}

Circle Circle::operator*(const double& k) const
{
    return Circle(center, radius * k);
}

Circle& Circle::operator*=(const double& k)
{
    this->radius *= k;
    return *this;
}

bool Circle::operator==(const Circle& c) const
{
    return (center == c.center) && (fabs(radius - c.radius) < eps);
}

bool Circle::operator!=(const Circle& c) const
{
    return (center != c.center) || (fabs(radius - c.radius) > eps);
}

/// Distance from the point to the circle
double Point::distance(Circle c) const
{
    return sqrt(pow(x - c.center.x, 2) + pow(y - c.center.y, 2));
}

/// Distance from the circle to the point
double Circle::distance(Point p) const
{
    return this->center.distance(p);
}

/// Distance from the circle to the line
double Circle::distance(Line l) const
{
    return this->center.distance(l);
}

/// Distance from the circle to the circle
double Circle::distance(Circle c) const
{
    return this->center.distance(c.center);
}

/// Finds intersection points P1 and P2 of the circle C and the line L.
/// Returns number of intersections.
int intersect(Circle c, Line l, Point &p1, Point &p2)
{
    double d = c.distance(l);

    // No intersection
    if (d > c.radius)
        return 0;

    // One intersection
    if (fabs(d - c.radius) < eps)
    {
        p1 = p2 = c.center.projection(l);
        return 1;
    }

    // Two intersections
    Point pr = c.center.projection(l);
    Vector g = l.guiding_vector().normalize() * sqrt(pow(c.radius, 2) - pow(d, 2));
    p1 = pr + g;
    p2 = pr - g;
    return 2;
}

int intersect(Line l, Circle c, Point &p1, Point &p2)
{
    return intersect(c, l, p1, p2);
}

/// Finds intersection points P1 and P2 of circles C1 and C2.
/// Returns -1, if circles are equal, or number of intersections.
int intersect(Circle c1, Circle c2, Point &p1, Point &p2)
{
    // Both circles are points
    if ((c1.radius < eps) && (c2.radius < eps))
    {
        if (c1.center == c2.center)
        {
            p1 = p2 = c1.center;
            return 1;
        }
        return 0;
    }

    double d = c1.distance(c2);

    // One of the circles is a point
    if ((c1.radius < eps) || (c2.radius < eps))
    {
        if (fabs(d - max(c1.radius, c2.radius)) < eps)
        {
            p1 = p2 = c1.radius < eps ? c1.center : c2.center;
            return 1;
        }
        return 0;
    }

    // Circles are equal
    if (c1 == c2)
        return -1;

    // Centers are equal
    if (c1.center == c2.center)
        return 0;

    Line l = Line(2 * (c2.center.x - c1.center.x),
                  2 * (c2.center.y - c1.center.y),
                  pow(c1.center.x, 2) + pow(c1.center.y, 2) - pow(c1.radius, 2) - (pow(c2.center.x, 2) + pow(c2.center.y, 2) - pow(c2.radius, 2)));
    return intersect(c1, l, p1, p2);
}

/// Finds intersection points P1 and P2 of the tangents of the circle C that pass through the point P.
/// Returns number of tangents.
int tangent(Point p, Circle c, Point &p1, Point &p2)
{
    double d = c.distance(p);

    // Point is inside the circle
    if (d < c.radius)
        return 0;

    // Point is on the circle
    if (fabs(d - c.radius) < eps)
    {
        p1 = p2 = p;
        return 1;
    }

    // Point is outside of the circle
    double angle = asin(c.radius / d);
    Point p0 = c.center.rotate(angle, p);
    p1 = p + (Vector(p, p0) / d * sqrt(pow(d, 2) - pow(c.radius, 2)));
    p2 = p1.rotate(-2.0 * angle, p);
    return 2;
}

