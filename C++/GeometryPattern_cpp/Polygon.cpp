#include "Polygon.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

void Polygon::add(Point p)
{
    vertices.push_back(p);
}

void Polygon::insert(Point p, size_t index)
{
    if (index > vertices.size())
        throw out_of_range("Invalid index");
    else if (index == vertices.size())
        add(p);
    else
        vertices.insert(vertices.begin() + index, p);
}

void Polygon::erase(size_t index)
{
    if (index >= vertices.size())
        throw out_of_range("Invalid index");
    vertices.erase(vertices.begin() + index);
}

void Polygon::pop()
{
    if (size() == 0)
        throw out_of_range("Polygon is empty");
    erase(size() - 1);
}

void Polygon::clear()
{
    vertices.clear();
}

bool Polygon::empty() const
{
    return vertices.size() == 0;
}

size_t Polygon::size() const
{
    return vertices.size();
}

void Polygon::resize(size_t new_size)
{
    vertices.resize(new_size);
}

void Polygon::reverse()
{
    std::reverse(vertices.begin(), vertices.end());
}

bool Polygon::convex() const
{
    int n = size();
    if (n <= 3)
        return true;

    double prev = cross_product(Vector(vertices[0], vertices[1]), Vector(vertices[1], vertices[2]));
    double now;
    for (int i = 1; i < n; ++i)
    {
        now = cross_product(Vector(vertices[i % n], vertices[(i + 1) % n]), Vector(vertices[(i + 1) % n], vertices[(i + 2) % n]));
        if (prev * now < 0)
            return false;
        prev = now;
    }
    return true;
}

double Polygon::area() const
{
    double sum = 0;
    for (size_t i = 1; i < size() - 1; ++i)
        sum += cross_product(Vector(vertices[0], vertices[i]), Vector(vertices[0], vertices[i + 1])) / 2.0;
    return fabs(sum);
}

istream &operator >>(istream &in, Polygon &poly)
{
    size_t n;
    in >> n;
    poly.resize(n);
    for (size_t i = 0; i < n; ++i)
        in >> poly[i];
    return in;
}

ostream &operator <<(ostream &out, const Polygon poly)
{
    out << poly.size();
    for (size_t i = 0; i < poly.size(); ++i)
        out << endl << poly.vertices[i];
    return out;
}

Point& Polygon::operator[](int index)
{
    return vertices[index];
}

Polygon convex_hull(vector<Point> points)
{
    points.erase(unique(points.begin(), points.end()), points.end());
    if (points.size() <= 2)
        return Polygon(points);

    Polygon convex = Polygon();
    sort(points.begin(), points.end(), [] (Point a, Point b)
    {
        return (a.y < b.y) || ((fabs(a.y - b.y) < eps) && (a.x < b.x));
    });
    convex.add(points[0]);
    Point p0 = convex[0];
    sort(points.begin(), points.end(), [p0] (Point a, Point b)
    {
        Point p1 = a - p0, p2 = b - p0;
        return (p1.polar_angle() < p2.polar_angle()) || ((fabs(p1.polar_angle() - p2.polar_angle()) < eps) && (p1.distance() < p2.distance()));
    });
    convex.add(points[1]);
    for (size_t i = 2; i < points.size(); ++i)
    {
        while ((convex.size() >= 2) && (cross_product(Vector(convex[convex.size() - 2], convex[convex.size() - 1]), Vector(convex[convex.size() - 1], points[i])) <= 0))
        {
            convex.pop();
        }
        convex.add(points[i]);
    }
    return convex;
}
