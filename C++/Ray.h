#ifndef RAY_H
#define RAY_H

#include <iostream>
#include <cmath>
#include "Point.h"
#include "Vector.h"
#include "Line.h"
#include "Segment.h"

class Ray
{
public:
    Point origin;
    Point guiding_point;

    Ray(Point origin, Point guiding_point)
    {
        if (origin == guiding_point)
            throw std::invalid_argument("Invalid guiding point");
        this->origin = origin;
        this->guiding_point = guiding_point;
    }

    Ray(Point origin, Vector guiding_vector) : Ray(origin, origin + guiding_vector) { }

    Ray (double x1, double y1, double x2, double y2) : Ray(Point(x1, y1), Point(x2, y2)) { }

    Ray() : Ray(0, 0, 1, 0) { }

    Vector guiding_vector() const;

    friend std::istream &operator >>(std::istream &in, Ray &r);

    friend std::ostream &operator <<(std::ostream &out, const Ray r);

    bool operator==(const Ray& r) const;

    bool operator!=(const Ray& r) const;
};

#endif // RAY_H
