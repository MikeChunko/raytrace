#include <iostream>

using namespace std;

// Vector in 3 dimensions
struct Vec {
    double x, y, z;

    Vec() : x(0), y(0), z(0) {};
    Vec(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {};

    Vec operator+(const Vec& v) const {
        return Vec(x + v.x, y + v.y, z + v.z);
    }

    Vec operator-(const Vec& v) const {
        return Vec(x - v.x, y - v.y, z - v.z);
    }

    Vec operator*(double c) const {
        return Vec(c * x, c * y, c * z);
    }

    Vec operator/(double c) const {
        return Vec(x / c, y / c, z / c);
    }
};

ostream& operator<<(ostream& os, Vec const& v) {
    os << '(' << v.x << ", " << v.y << ", " << v.z << ')';
    return os;
}

// Raytracing ray
struct Ray {
    Vec origin, direction;

    Ray() { origin = Vec(); direction = Vec();}
    Ray(const Vec& o, const Vec&d) : origin(o), direction(d) {};
};

ostream& operator<<(ostream& os, Ray const& r) {
    os << "[Origin: " << r.origin << ", Direction: " << r.direction << "]";
    return os;
}


int main() {
    Ray test = Ray(Vec(1, 2, 3), Vec(0, 0, 1));
    cout << test << endl;
    return 0;
}