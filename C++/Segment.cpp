#include "Segment.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

Vector Segment::guiding_vector() const
{
    return Vector(p1, p2);
}

/// Distance from the segment to the segment
double Segment::distance(Segment s) const
{
    Segment tmp;
    if (intersect(*this, s, tmp) > 0)
        return 0;
    return min(min(p1.distance(s), p2.distance(s)), min(s.p1.distance(*this), s.p2.distance(*this)));
}

istream &operator >>(istream &in, Segment &s)
{
    in >> s.p1 >> s.p2;
    return in;
}

ostream &operator <<(ostream &out, const Segment s)
{
    out << s.p1 << ' ' << s.p2;
    return out;
}

bool Segment::operator==(const Segment& s) const
{
    return ((p1 == s.p1) && (p2 == s.p2)) || ((p1 == s.p2) && (p2 == s.p1));
}

bool Segment::operator!=(const Segment& s) const
{
    return ((p1 != s.p1) || (p2 != s.p2)) && ((p1 != s.p2) || (p2 != s.p1));
}

bool Point::belongs(Segment s) const
{
    if (s.p1 == s.p2)
        return s.p1 == *this;
    return belongs(Line(s.p1, s.p2)) && (dot_product(Vector(s.p1, s.p2), Vector(s.p1, *this)) >= 0) && (dot_product(Vector(s.p2, s.p1), Vector(s.p2, *this)) >= 0);
}

/// Distance from the point to the segment
double Point::distance(Segment s) const
{
    if (s.p1 == s.p2)
        return distance(s.p1);
    if ((dot_product(Vector(s.p1, s.p2), Vector(s.p1, *this)) >= 0) && (dot_product(Vector(s.p2, s.p1), Vector(s.p2, *this)) >= 0))
    {
        return distance(Line(s.p1, s.p2));
    }
    return min(distance(s.p1), distance(s.p2));
}

/// Finds intersection segment of two segments S1 and S2.
/// Returns -2, if segment intersection is a segment; 0, if segments doesn't intersect; 1, if segments intersect in one point.
int intersect(Segment s1, Segment s2, Segment &intersection)
{
    // One of the segments is a point
    if (s1.p1 == s1.p2)
    {
        if (s1.p1.belongs(s2))
        {
            intersection.p1 = intersection.p2 = s1.p1;
            return 1;
        }
        return 0;
    }
    if (s2.p1 == s2.p2)
    {
        if (s2.p1.belongs(s1))
        {
            intersection.p1 = intersection.p2 = s2.p1;
            return 1;
        }
        return 0;
    }

    int k = intersect(Line(s1.p1, s1.p2), Line(s2.p1, s2.p2), intersection.p1);
    intersection.p2 = intersection.p1;

    if (k == 1)
    {
        if (!intersection.p1.belongs(s1) || !intersection.p1.belongs(s2))
            k = 0;
    }
    else if (k == -1)
    {
        // S1 is in S2
        if (s1.p1.belongs(s2) && s1.p2.belongs(s2))
        {
            intersection = s1;
            k = -2;
        }
        // S2 is in S1
        else if (s2.p1.belongs(s1) && s2.p2.belongs(s1))
        {
            intersection = s2;
            k = -2;
        }

        // S1 and S2 intersect by their points
        else if (s1.p1.belongs(s2))
        {
            if (s2.p1.belongs(s1))
                intersection = Segment(s1.p1, s2.p1);
            else
                intersection = Segment(s1.p1, s2.p2);

            // Single point intersection
            if (intersection.p1 == intersection.p2)
                k = 1;
            else
                k = -2;
        }
        else if (s1.p2.belongs(s2))
        {
            if (s2.p1.belongs(s1))
                intersection = Segment(s1.p2, s2.p1);
            else
                intersection = Segment(s1.p2, s2.p2);

            // Single point intersection
            if (intersection.p1 == intersection.p2)
                k = 1;
            else
                k = -2;
        }
        else
            k = 0;
    }
    return k;
}
