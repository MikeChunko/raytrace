#include <iostream>
#include <fstream>
#include <cmath>

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

    // Return a normalized unit vector
    Vec normalize() const {
        double mean = sqrt(x*x + y*y + z*z);
        return Vec(x/mean, y/mean, z/mean);
    }
};

// Return the dot product of two vectors
double dot(const Vec& a, const Vec& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

ostream& operator<<(ostream& os, Vec const& v) {
    os << '(' << v.x << ", " << v.y << ", " << v.z << ')';
    return os;
}

// Raytracing ray
struct Ray {
    Vec origin, direction;

    Ray() { origin = Vec(); direction = Vec(); }
    Ray(const Vec& o, const Vec& d) : origin(o), direction(d) {};
};

ostream& operator<<(ostream& os, Ray const& r) {
    os << "[Origin: " << r.origin << ", Direction: " << r.direction << "]";
    return os;
}

// Represent a sphere
// Used in scene creation
struct Sphere {
    Vec center, color;
    double radius;

    Sphere() { center = Vec(); color = Vec(); radius = 0; }
    Sphere(const Vec& c, const Vec& col, double r) : center(c), color(col), radius(r) {}

    bool intersect(const Ray& ray, double &t) const {
        const Vec oc = ray.origin - center;
        const double b = 2 * dot(oc, ray.direction),
        			 c = dot(oc, oc) - radius*radius;
        double disc = b*b - 4*c;

        if (disc < 1e-3)  // Reasonably close to 0 should be accepted as being 0
            return false;

        disc = sqrt(disc);
        const double t0 = -b - disc, t1 = -b + disc;

        t = (t0 < t1) ? t0 : t1;
        return true;
    }

    Vec normal(const Vec& p) const {
        return (p - center) / radius;
    }
};

// Clamp a color vector to 8bit colors
void color_clamp(Vec& v) {
    v.x = (v.x > 255) ? 255 : (v.x < 0) ? 0 : v.x;
    v.y = (v.y > 255) ? 255 : (v.y < 0) ? 0 : v.y;
    v.z = (v.z > 255) ? 255 : (v.z < 0) ? 0 : v.z;
}


int main() {
    const int HEIGHT = 200, WIDTH = 200;  // Image height and width
    const Vec white(255,255,255), // Quick reference colors
              black(0,0,0),
              green(0,255,0);

    // Header for .ppm
    ofstream out = ofstream("result.ppm");
    out << "P3\n" << WIDTH << ' ' << HEIGHT << " 255\n";

    // Scene creation
    const Sphere sphere(Vec(.5*WIDTH, .5*HEIGHT, 50), green, .5*WIDTH),  // Sphere centered and with diameter half as wide as the screen
                 light(Vec(1.5*WIDTH,1.5*HEIGHT, 50), white, 1);  // Point light source at the same depth as the sphere

	double t;
	Vec color;

    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            color = black;  // Initialize to black (nothing)
            const Ray ray(Vec(x,y,0), Vec(0,0,1));  // Initial ray

            if (sphere.intersect(ray, t)) {  // Do some funky math (not sure on how this works)
				const Vec p = ray.origin + ray.direction*t,
						  l = light.center - p,
						  n = sphere.normal(p);
				const double dt = dot(l.normalize(), n.normalize());

                color = (sphere.color + light.color*dt) * .5;
			}

            // Output color at current pixel
            color_clamp(color);
            out << (int)color.x << ' ' << (int)color.y << ' ' << (int)color.z << "\n";
        }
        cout << endl;
    }
}