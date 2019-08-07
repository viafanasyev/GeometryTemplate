#include "Ray.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

Vector Ray::guiding_vector() const
{
    return Vector(origin, guiding_point);
}

istream &operator >>(istream &in, Ray &r)
{
    in >> r.origin >> r.guiding_point;
    return in;
}

ostream &operator <<(ostream &out, const Ray r)
{
    out << r.origin << ' ' << r.guiding_point;
    return out;
}

bool Ray::operator==(const Ray& r) const
{
    return (guiding_vector().normalize() == r.guiding_vector().normalize()) && (origin == r.origin);
}

bool Ray::operator!=(const Ray& r) const
{
    return (guiding_vector().normalize() != r.guiding_vector().normalize()) || (origin != r.origin);
}

bool Point::belongs(Ray r) const
{
    bool b = belongs(Line(r.origin, r.guiding_point));
    b &= dot_product(r.guiding_vector(), Vector(r.origin, *this)) >= 0;
    return b;
}

/// Distance from the point to the ray
double Point::distance(Ray r) const
{
    if (dot_product(r.guiding_vector(), Vector(r.origin, *this)) >= 0)
        return distance(Line(r.origin, r.guiding_point));
    return distance(r.origin);
}

/// Finds intersection segment of two rays R1 and R2.
/// Returns -2, if rays intersection is a segment; -1, if rays intersection is a line; 0, if rays doesn't intersect; 1, if rays intersect in one point.
int intersect(Ray r1, Ray r2, Segment &intersection)
{
    int k = intersect(Line(r1.origin, r1.guiding_point), Line(r2.origin, r2.guiding_point), intersection.p1);
    intersection.p2 = intersection.p1;
    if (k == 1)
    {
        if (!intersection.p1.belongs(r1) || !intersection.p1.belongs(r2))
            k = 0;
    }
    else if (k == -1)
    {
        if (dot_product(r1.guiding_vector(), r2.guiding_vector()) < 0)
        {
            intersection.p1 = r1.origin;
            intersection.p2 = r2.origin;

            if (r1.origin == r2.origin)
                k = 1;
            else
                k = -2;
        }
    }
    return k;
}
