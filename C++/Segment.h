#ifndef SEGMENT_H
#define SEGMENT_H

#include <iostream>
#include "Point.h"
#include "Vector.h"
#include "Line.h"

class Segment
{
public:
    Point p1, p2;

    Segment(Point p1, Point p2)
    {
        this->p1 = p1;
        this->p2 = p2;
    }

    Segment(double x1, double y1, double x2, double y2) : Segment(Point(x1, y1), Point(x2, y2)) { }

    Segment() : Segment(0, 0, 0, 0) { }

    Vector guiding_vector() const;

    /// Distance from the segment to the segment
    double distance(Segment s) const;

    friend std::istream &operator >>(std::istream &in, Segment &s);

    friend std::ostream &operator <<(std::ostream &out, const Segment s);

    bool operator==(const Segment& s) const;

    bool operator!=(const Segment& s) const;
};

#endif // SEGMENT_H
