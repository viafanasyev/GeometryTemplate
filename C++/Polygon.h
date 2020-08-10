#ifndef POLYGON_H
#define POLYGON_H

#include <iostream>
#include <vector>
#include "Point.h"
#include "Vector.h"

Polygon convex_hull(std::vector<Point> points);

class Polygon
{
private:
    std::vector<Point> vertices;

public:
    Polygon(std::vector<Point> vertices)
    {
        this->vertices = std::vector<Point>(vertices);
    }

    Polygon() : Polygon(std::vector<Point>(0)) { }

    void add(Point p);

    void insert(Point p, size_t index);

    void erase(size_t index);

    void pop();

    void clear();

    bool empty() const;

    size_t size() const;

    void resize(size_t new_size);

    void reverse();

    bool convex() const;

    double area() const;

    friend std::istream &operator >>(std::istream &in, Polygon &poly);

    friend std::ostream &operator <<(std::ostream &out, const Polygon poly);

    Point& operator[](int index);
};

#endif // POLYGON_H
