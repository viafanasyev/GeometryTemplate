#ifndef CIRCLE_H
#define CIRCLE_H

#include <iostream>
#include "Point.h"
#include "Vector.h"
#include "Line.h"

class Circle
{
public:
    Point center;
    double radius;

    Circle(Point center, double radius)
    {
        this->center = center;
        this->radius = radius;
    }

    Circle(double x, double y, double radius) : Circle(Point(x, y), radius) { }

    Circle() : Circle(Point(0, 0), 0) { }

    double distance(Point p) const;
    double distance(Line l) const;
    double distance(Circle c) const;

    friend std::istream &operator >>(std::istream &in, Circle &c);

    friend std::ostream &operator <<(std::ostream &out, const Circle c);

    Circle operator+(const Vector& v) const;

    Circle& operator+=(const Vector& v);

    Circle operator-(const Vector& v) const;

    Circle& operator-=(const Vector& v);

    Circle operator*(const double& k) const;

    Circle& operator*=(const double& k);

    bool operator==(const Circle& c) const;

    bool operator!=(const Circle& c) const;
};

#endif // CIRCLE_H
