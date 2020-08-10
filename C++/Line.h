#ifndef LINE_H
#define LINE_H

#include <iostream>
#include <cmath>
#include "Point.h"
#include "Vector.h"

/// True, if lines are parallel
bool parallel(Line l1, Line l2);

class Line
{
public:
    double a, b, c;

    Line(double a, double b, double c)
    {
        if ((fabs(a) < eps) && (fabs(b) < eps))
            throw std::invalid_argument("Invalid line exception");
        this->a = a;
        this->b = b;
        this->c = c;
    }

    Line(Point p1, Point p2)
    {
        if (p1 == p2)
            throw std::invalid_argument("Invalid line exception");
        a = p1.y - p2.y;
        b = p2.x - p1.x;
        c = -(a * p1.x + b * p1.y);
    }

    Line(double x1, double y1, double x2, double y2) : Line(Point(x1, y1), Point(x2, y2)) { }

    Line(Vector n, Point p)
    {
        a = n.x;
        b = n.y;
        c = -(a * p.x + b * p.y);
    }

    Line() : Line(0, 1, 0) { }

    /// Normal vector of the line
    Vector normal_vector() const;

    /// Guiding vector of the line
    Vector guiding_vector() const;

    /// Normalizes the line
    Line normalize() const;

    friend std::istream &operator >>(std::istream &in, Line &l);

    friend std::ostream &operator <<(std::ostream &out, const Line l);

    Line operator+(const double& k) const;

    Line& operator+=(const double& k);

    Line operator-(const double& k) const;

    Line& operator-=(const double& k);

    bool operator==(const Line& l) const;

    bool operator!=(const Line& l) const;
};

#endif // LINE_H
