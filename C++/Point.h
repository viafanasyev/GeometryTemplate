#ifndef POINT_H
#define POINT_H

#include <iostream>

const double eps = 1e-9;

class Point;
class Vector;
class Line;
class Ray;
class Segment;
class Circle;
class Polygon;

int intersect(Line l1, Line l2, Point &p);
int intersect(Ray r1, Ray r2, Segment &intersection);
int intersect(Segment s1, Segment s2, Segment &intersection);
int intersect(Line l, Circle c, Point &p1, Point &p2);
int intersect(Circle c, Line l, Point &p1, Point &p2);
int intersect(Circle c1, Circle c2, Point &p1, Point &p2);

class Point
{
public:
    double x, y;

    Point(double x, double y)
    {
        this->x = x;
        this->y = y;
    }

    Point() : Point(0, 0) { }

    /// Polar angle of the point
    double polar_angle() const;

    /// Rotates the point on angle around the point p
    Point rotate(double angle, Point p = Point(0, 0)) const;


    /// Distance between two points
    double distance(Point p = Point(0, 0)) const;

    double distance(Line l) const;
    double distance(Ray r) const;
    double distance(Segment s) const;
    double distance(Circle c) const;

    double sign_distance(Line l) const;

    Point projection(Line l) const;

    bool belongs(Line l) const;
    bool belongs(Ray r) const;
    bool belongs(Segment s) const;

    friend std::istream &operator >>(std::istream &in, Point &p);

    friend std::ostream &operator <<(std::ostream &out, const Point p);

    Point operator+(const Point& p) const;

    Point& operator+=(const Point& p);

    Point operator+(const Vector& v) const;

    Point& operator+=(const Vector& v);

    Point operator-(const Point& p) const;

    Point& operator-=(const Point& p);

    Point operator-(const Vector& v) const;

    Point& operator-=(const Vector& v);

    Point operator*(const double& k) const;

    Point& operator*=(const double& k);

    Point operator/(const double& k) const;

    Point& operator/=(const double& k);

    bool operator==(const Point& p) const;

    bool operator!=(const Point& p) const;
};

#endif // POINT_H
