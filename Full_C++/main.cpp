#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

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
    double polar_angle() const
    {
        double angle = atan2(y, x);
        if (angle < 0)
        {
            angle += 2 * M_PI;
        }
        return angle;
    }

    /// Rotates the point on angle around the point p
    Point rotate(double angle, Point p = Point(0, 0)) const
    {
        Point tmp;
        tmp.x = (x - p.x) * cos(angle) - (y - p.y) * sin(angle) + p.x;
        tmp.y = (y - p.y) * cos(angle) + (x - p.x) * sin(angle) + p.y;
        return tmp;
    }


    /// Distance between two points
    double distance(Point p = Point(0, 0)) const
    {
        return hypot(fabs(p.x - x), fabs(p.y - y));
    }

    double distance(Line l) const;
    double distance(Ray r) const;
    double distance(Segment s) const;
    double distance(Circle c) const;

    double sign_distance(Line l) const;

    Point projection(Line l) const;

    bool belongs(Line l) const;
    bool belongs(Ray r) const;
    bool belongs(Segment s) const;

    friend istream &operator >>(istream &in, Point &p)
    {
        in >> p.x >> p.y;
        return in;
    }

    friend ostream &operator <<(ostream &out, const Point p)
    {
        out << p.x << ' ' << p.y;
        return out;
    }

    Point operator+(const Point& p) const
    {
        return Point(x + p.x, y + p.y);
    }

    Point& operator+=(const Point& p)
    {
        this->x += p.x;
        this->y += p.y;
        return *this;
    }

    friend Point& operator+=(Point& p, const Vector& v);

    Point operator-(const Point& p) const
    {
        return Point(x - p.x, y - p.y);
    }

    Point& operator-=(const Point& p)
    {
        this->x -= p.x;
        this->y -= p.y;
        return *this;
    }

    friend Point& operator-=(Point& p, const Vector& v);

    Point operator*(const double& k) const
    {
        return Point(x * k, y * k);
    }

    Point& operator*=(const double& k)
    {
        this->x *= k;
        this->y *= k;
        return *this;
    }

    Point operator/(const double& k) const
    {
        if (fabs(k) < eps)
            throw overflow_error("Division by zero exception");
        return Point(x / k, y / k);
    }

    Point& operator/=(const double& k)
    {
        if (fabs(k) < eps)
            throw overflow_error("Division by zero exception");
        this->x /= k;
        this->y /= k;
        return *this;
    }

    bool operator==(const Point& p) const
    {
        return (fabs(x - p.x) < eps) && (fabs(y - p.y) < eps);
    }

    bool operator!=(const Point& p) const
    {
        return (fabs(x - p.x) > eps) || (fabs(y - p.y) > eps);
    }
};

class Vector
{
public:
    double x, y;

    Vector(double x, double y)
    {
        this->x = x;
        this->y = y;
    }

    Vector(Point p) : Vector(p.x, p.y) { }

    Vector(double x_begin, double y_begin, double x_end, double y_end) : Vector(x_end - x_begin, y_end - y_begin) { }

    Vector(Point begin, Point end) : Vector(end - begin) { }

    Vector() : Vector(0, 0) { }

    /// Length of the vector
    double length() const
    {
        return hypot(x, y);
    }

    /// Polar angle of the vector
    double polar_angle() const
    {
        double angle = atan2(y, x);
        if (angle < 0)
        {
            angle += 2 * M_PI;
        }
        return angle;
    }

    /// Normalizes the vector (divides by it's length)
    Vector normalize() const
    {
        if ((fabs(x) < eps) && (fabs(y) < eps)) // If Vector = {0, 0}
            return *this;

        double l = length();
        return Vector(x / l, y / l);
    }

    /// Rotates the vector on angle around the point p
    Vector rotate(double angle, Point p = Point(0, 0)) const
    {
        Vector tmp;
        tmp.x = (x - p.x) * cos(angle) - (y - p.y) * sin(angle) + p.x;
        tmp.y = (y - p.y) * cos(angle) + (x - p.x) * sin(angle) + p.y;
        return tmp;
    }

    friend istream &operator >>(istream &in, Vector &v)
    {
        in >> v.x >> v.y;
        return in;
    }

    friend ostream &operator <<(ostream &out, const Vector v)
    {
        out << v.x << ' ' << v.y;
        return out;
    }

    Vector operator-() const
    {
        return Vector(-x, -y);
    }

    Vector operator+(const Vector& v) const
    {
        return Vector(x + v.x, y + v.y);
    }

    Vector& operator+=(const Vector& v)
    {
        this->x += v.x;
        this->y += v.y;
        return *this;
    }

    Vector operator-(const Vector& v) const
    {
        return Vector(x - v.x, y - v.y);
    }

    Vector& operator-=(const Vector& v)
    {
        this->x -= v.x;
        this->y -= v.y;
        return *this;
    }

    Vector operator*(const double& k) const
    {
        return Vector(x * k, y * k);
    }

    Vector& operator*=(const double& k)
    {
        this->x *= k;
        this->y *= k;
        return *this;
    }

    Vector operator/(const double& k) const
    {
        if (fabs(k) < eps)
            throw overflow_error("Division by zero exception");
        return Vector(x / k, y / k);
    }

    Vector& operator/=(const double& k)
    {
        if (fabs(k) < eps)
            throw overflow_error("Division by zero exception");
        this->x /= k;
        this->y /= k;
        return *this;
    }

    bool operator==(const Vector& v) const
    {
        return (fabs(x - v.x) < eps) && (fabs(y - v.y) < eps);
    }

    bool operator!=(const Vector& v) const
    {
        return (fabs(x - v.x) > eps) || (fabs(y - v.y) > eps);
    }
};

/// Scalar product of two vectors
double dot_product(Vector v1, Vector v2)
{
    return v1.x * v2.x + v1.y * v2.y;
}

/// Vector product of two vectors
double cross_product(Vector v1, Vector v2)
{
    return v1.x * v2.y - v1.y * v2.x;
}

/// Angle between two vectors
double vectors_angle(Vector v1, Vector v2)
{
    if ((v1 == Vector(0, 0)) || (v2 == Vector(0, 0)))
        return M_PI / 2.0;
    return acos(dot_product(v1, v2) / (v1.length() * v2.length()));
}

Point operator+(const Point& p, const Vector& v)
{
    return Point(p.x + v.x, p.y + v.y);
}

Point& operator+=(Point& p, const Vector& v)
{
    p.x += v.x;
    p.y += v.y;
    return p;
}

Point& operator-=(Point& p, const Vector& v)
{
    p.x -= v.x;
    p.y -= v.y;
    return p;
}

Point operator-(const Point& p, const Vector& v)
{
    return Point(p.x - v.x, p.y - v.y);
}

class Line
{
public:
    double a, b, c;

    Line(double a, double b, double c)
    {
        if ((fabs(a) < eps) && (fabs(b) < eps))
            throw invalid_argument("Invalid line exception");
        this->a = a;
        this->b = b;
        this->c = c;
    }

    Line(Point p1, Point p2)
    {
        if (p1 == p2)
            throw invalid_argument("Invalid line exception");
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
    Vector normal_vector() const
    {
        return Vector(a, b);
    }

    /// Guiding vector of the line
    Vector guiding_vector() const
    {
        return Vector(b, -a);
    }

    /// Normalizes the line
    Line normalize() const
    {
        double l = normal_vector().length();
        return Line(a / l, b / l, c / l);
    }

    friend istream &operator >>(istream &in, Line &l)
    {
        in >> l.a >> l.b >> l.c;
        return in;
    }

    friend ostream &operator <<(ostream &out, const Line l)
    {
        out << l.a << ' ' << l.b << ' ' << l.c;
        return out;
    }

    Line operator+(const double& k) const
    {
        return Line(a, b, c + k);
    }

    Line& operator+=(const double& k)
    {
        this->c += k;
        return *this;
    }

    Line operator-(const double& k) const
    {
        return Line(a, b, c - k);
    }

    Line& operator-=(const double& k)
    {
        this->c -= k;
        return *this;
    }

    bool operator==(const Line& l) const
    {
        Line l1 = this->normalize(), l2 = l.normalize();
        return ((fabs(l1.a - l2.a) < eps) && (fabs(l1.b - l2.b) < eps) && (fabs(l1.c - l2.c) < eps))
                || ((fabs(l1.a + l2.a) < eps) && (fabs(l1.b + l2.b) < eps) && (fabs(l1.c + l2.c) < eps));
    }

    bool operator!=(const Line& l) const
    {
        Line l1 = this->normalize(), l2 = l.normalize();
        return ((fabs(l1.a - l2.a) > eps) || (fabs(l1.b - l2.b) > eps) || (fabs(l1.c - l2.c) > eps))
                && ((fabs(l1.a + l2.a) > eps) || (fabs(l1.b + l2.b) > eps) || (fabs(l1.c + l2.c) > eps));
    }
};

/// True, if lines are parallel
bool parallel(Line l1, Line l2)
{
    return (fabs(l1.a - l2.a) < eps) && (fabs(l1.b - l2.b) < eps);
}

bool Point::belongs(Line l) const
{
    return fabs(l.a * x + l.b * y + l.c) < eps;
}

class Ray
{
public:
    Point origin;
    Point guiding_point;

    Ray(Point origin, Point guiding_point)
    {
        if (origin == guiding_point)
            throw invalid_argument("Invalid guiding point");
        this->origin = origin;
        this->guiding_point = guiding_point;
    }

    Ray(Point origin, Vector guiding_vector) : Ray(origin, origin + guiding_vector) { }

    Ray (double x1, double y1, double x2, double y2) : Ray(Point(x1, y1), Point(x2, y2)) { }

    Ray() : Ray(0, 0, 1, 0) { }

    Vector guiding_vector() const
    {
        return Vector(origin, guiding_point);
    }

    friend istream &operator >>(istream &in, Ray &r)
    {
        in >> r.origin >> r.guiding_point;
        return in;
    }

    friend ostream &operator <<(ostream &out, const Ray r)
    {
        out << r.origin << ' ' << r.guiding_point;
        return out;
    }

    bool operator==(const Ray& r) const
    {
        return (guiding_vector().normalize() == r.guiding_vector().normalize()) && (origin == r.origin);
    }

    bool operator!=(const Ray& r) const
    {
        return (guiding_vector().normalize() != r.guiding_vector().normalize()) || (origin != r.origin);
    }
};

bool Point::belongs(Ray r) const
{
    bool b = belongs(Line(r.origin, r.guiding_point));
    b &= dot_product(r.guiding_vector(), Vector(r.origin, *this)) >= 0;
    return b;
}

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

    Vector guiding_vector() const
    {
        return Vector(p1, p2);
    }

    /// Distance from the segment to the segment
    double distance(Segment s) const;

    friend istream &operator >>(istream &in, Segment &s)
    {
        in >> s.p1 >> s.p2;
        return in;
    }

    friend ostream &operator <<(ostream &out, const Segment s)
    {
        out << s.p1 << ' ' << s.p2;
        return out;
    }

    bool operator==(const Segment& s) const
    {
        return ((p1 == s.p1) && (p2 == s.p2)) || ((p1 == s.p2) && (p2 == s.p1));
    }

    bool operator!=(const Segment& s) const
    {
        return ((p1 != s.p1) || (p2 != s.p2)) && ((p1 != s.p2) || (p2 != s.p1));
    }
};

bool Point::belongs(Segment s) const
{
    if (s.p1 == s.p2)
        return s.p1 == *this;
    return belongs(Line(s.p1, s.p2)) && (dot_product(Vector(s.p1, s.p2), Vector(s.p1, *this)) >= 0) && (dot_product(Vector(s.p2, s.p1), Vector(s.p2, *this)) >= 0);
}

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

    friend istream &operator >>(istream &in, Circle &c)
    {
        in >> c.center >> c.radius;
        return in;
    }

    friend ostream &operator <<(ostream &out, const Circle c)
    {
        out << c.center << ' ' << c.radius;
        return out;
    }

    Circle operator+(const Vector& v) const
    {
        return Circle(center + v, radius);
    }

    Circle& operator+=(const Vector& v)
    {
        this->center += v;
        return *this;
    }

    Circle operator-(const Vector& v) const
    {
        return Circle(center - v, radius);
    }

    Circle& operator-=(const Vector& v)
    {
        this->center -= v;
        return *this;
    }

    Circle operator*(const double& k) const
    {
        return Circle(center, radius * k);
    }

    Circle& operator*=(const double& k)
    {
        this->radius *= k;
        return *this;
    }

    bool operator==(const Circle& c) const
    {
        return (center == c.center) && (fabs(radius - c.radius) < eps);
    }

    bool operator!=(const Circle& c) const
    {
        return (center != c.center) || (fabs(radius - c.radius) > eps);
    }
};

/// Distance from the point to the line
double Point::distance(Line l) const
{
    return fabs(l.a * x + l.b * y + l.c) / l.normal_vector().length();
}

/// Sign distance from the point to the line
double Point::sign_distance(Line l) const
{
    return (l.a * x + l.b * y + l.c) / l.normal_vector().length();
}

/// Distance from the point to the ray
double Point::distance(Ray r) const
{
    if (dot_product(r.guiding_vector(), Vector(r.origin, *this)) >= 0)
        return distance(Line(r.origin, r.guiding_point));
    return distance(r.origin);
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

/// Distance from the point to the circle
double Point::distance(Circle c) const
{
    return sqrt(pow(x - c.center.x, 2) + pow(y - c.center.y, 2));
}

/// Distance from the segment to the segment
double Segment::distance(Segment s) const
{
    Segment tmp;
    if (intersect(*this, s, tmp) > 0)
        return 0;
    return min(min(p1.distance(s), p2.distance(s)), min(s.p1.distance(*this), s.p2.distance(*this)));
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

/// Projection of the point on the line
Point Point::projection(Line l) const
{
    Vector n = (-l.normal_vector().normalize()) * sign_distance(l);
    return (*this + n);
}

/// Finds intersection point P of the two lines L1 and L2.
/// Returns -1, if lines are equal; 0, if lines are parallel; 1, if lines intersect.
int intersect(Line l1, Line l2, Point &p)
{
    if (l1 == l2) return -1;
    if (parallel(l1, l2)) return 0;

    p.x = (l2.b * l1.c - l1.b * l2.c) / (l2.a * l1.b - l1.a * l2.b);
    p.y = (l1.b != 0 ? (- l1.c - l1.a * p.x) / l1.b : (- l2.c - l2.a * p.x) / l2.b);
    return 1;
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

class Polygon
{
private:
    vector<Point> vertices;

public:
    Polygon(vector<Point> vertices)
    {
        this->vertices = vector<Point>(vertices);
    }

    Polygon() : Polygon(vector<Point>(0)) { }

    void add(Point p)
    {
        vertices.push_back(p);
    }

    void insert(Point p, size_t index)
    {
        if (index > vertices.size())
            throw out_of_range("Invalid index");
        else if (index == vertices.size())
            add(p);
        else
            vertices.insert(vertices.begin() + index, p);
    }

    void erase(size_t index)
    {
        if (index >= vertices.size())
            throw out_of_range("Invalid index");
        vertices.erase(vertices.begin() + index);
    }

    void pop()
    {
        if (size() == 0)
            throw out_of_range("Polygon is empty");
        erase(size() - 1);
    }

    void clear()
    {
        vertices.clear();
    }

    bool empty() const
    {
        return vertices.size() == 0;
    }

    size_t size() const
    {
        return vertices.size();
    }

    void resize(size_t new_size)
    {
        vertices.resize(new_size);
    }

    void reverse()
    {
        std::reverse(vertices.begin(), vertices.end());
    }

    bool convex() const
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

    double area() const
    {
        double sum = 0;
        for (size_t i = 1; i < size() - 1; ++i)
            sum += cross_product(Vector(vertices[0], vertices[i]), Vector(vertices[0], vertices[i + 1])) / 2.0;
        return fabs(sum);
    }

    friend istream &operator >>(istream &in, Polygon &poly)
    {
        size_t n;
        in >> n;
        poly.resize(n);
        for (size_t i = 0; i < n; ++i)
            in >> poly[i];
        return in;
    }

    friend ostream &operator <<(ostream &out, const Polygon poly)
    {
        out << poly.size();
        for (size_t i = 0; i < poly.size(); ++i)
            out << endl << poly.vertices[i];
        return out;
    }

    Point& operator[](int index)
    {
        return vertices[index];
    }
};

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

int main()
{
    return 0;
}
