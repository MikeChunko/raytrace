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
    Ray(const Vec& o, const Vec& d) : origin(o), direction(d) {};
};

ostream& operator<<(ostream& os, Ray const& r) {
    os << "[Origin: " << r.origin << ", Direction: " << r.direction << "]";
    return os;
}

// Clamp a color vector to 8bit colors
void color_clamp (Vec& v) {
    v.x = (v.x > 255) ? 255 : (v.x < 0) ? 0 : v.x;
    v.y = (v.y > 255) ? 255 : (v.y < 0) ? 0 : v.y;
    v.z = (v.z > 255) ? 255 : (v.z < 0) ? 0 : v.z;
}


int main() {
    const int HEIGHT = 20, WIDTH = 20;  // Image height and width

    // TODO: Scene creation

    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            Ray ray = Ray(Vec(x,y,0), Vec(0,0,1));  // Initial ray
            Vec color = Vec(0,0,0);  // Initialize to black (nothing)

            // TODO: Turn pseudocode into real code
            /*
            intersecting_obj = find_intersect(&objects, ray);
            if (intersecting_obj != null) {
                some ray tracing maginc
            }
            */

           cout << color;
           if (x < WIDTH - 1)
            cout << ", ";
        }
        cout << endl;
    }
}