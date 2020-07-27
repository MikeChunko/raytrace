#include <fstream>
#include <vector>
#include <limits>
#include <cmath>

using namespace std;

// Vector in 3 dimensions
struct Vec {
    double x, y, z;

    Vec() : x(0), y(0), z(0) {};
    Vec(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {};

    // Vector operations
    Vec operator+(const Vec& v) const {
        return Vec(x + v.x, y + v.y, z + v.z);
    }
    Vec operator-(const Vec& v) const {
        return Vec(x - v.x, y - v.y, z - v.z);
    }

    // Scalar operations
    Vec operator*(double c) const {
        return Vec(c * x, c * y, c * z);
    }
    Vec operator/(double c) const {
        return Vec(x / c, y / c, z / c);
    }

    // Return a normalized unit vector (magnitude of 1)
    Vec normalize() const {
        double length = sqrt(x*x + y*y + z*z);
        return Vec(x/length, y/length, z/length);
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

    // Return true if ray intersects with the sphere
    // If true, t will hold the minimum intersecting t-value
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

    // Return the vector normal to the circle and p
    Vec normal(const Vec& p) const {
        return (p - center).normalize();
    }
};

// Return the initial ray center at the given (x, y) coordinates
// Use perspective to determine the direction
Ray get_initial_ray(int x, int y, int width, int height, int fov=30) {
    double invWidth = 1 / double(width), invHeight = 1 / double(height),
           aspectratio = width / double(height),
           angle = tan(M_PI * 0.5 * fov / 180.),
           xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio,
           yy = (2 * ((y + 0.5) * invHeight) - 1) * angle;
    return Ray(Vec(x,y,0), Vec(xx,yy,1).normalize());
}

// Clamp a color vector to 8bit colors
// (Nothing below 0 or above 255)
void color_clamp(Vec& v) {
    v.x = (v.x > 255) ? 255 : (v.x < 0) ? 0 : v.x;
    v.y = (v.y > 255) ? 255 : (v.y < 0) ? 0 : v.y;
    v.z = (v.z > 255) ? 255 : (v.z < 0) ? 0 : v.z;
}

int main() {
    const int WIDTH = 600, HEIGHT = 600;  // Image width and height
    const Vec white(255,244,229),  // Quick reference colors
              black(0,0,0),
              red(230,0,0),
              green(0,230,0),
              blue(0,0,230);

    // Header for .ppm
    ofstream out = ofstream("result.ppm");
    out << "P3\n" << WIDTH << ' ' << HEIGHT << " 255\n";

    // Scene creation
    // Assume objects,back() is the light source
    vector<Sphere> objects; // Assume every obect in the sceneis a sphere for now
    objects.push_back(Sphere(Vec(.5*WIDTH, .5*HEIGHT, .425*WIDTH), green, .4*WIDTH));  // Sphere centered and with diameter half as wide as the screen
    objects.push_back(Sphere(Vec(.55*WIDTH, .3*HEIGHT, .04*WIDTH), red, .2*WIDTH));
    objects.push_back(Sphere(Vec(.35*WIDTH, .7*HEIGHT, .04*WIDTH), blue, .2*WIDTH));
    objects.push_back(Sphere(Vec(1.5*WIDTH, 1.5*HEIGHT, .25*WIDTH), white, 1));  // Point light source at the same depth as the sphere

	Vec color;
    Ray ray;
	double t;
    int min_t;

    // Render each pixel
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            // Initialize variables
            color = black;
            ray = get_initial_ray(x,y,WIDTH,HEIGHT);
            min_t = numeric_limits<int>::max();
            Sphere min_obj;

            // Find closest interecting object
            for (auto obj = objects.begin(); obj < objects.end() - 1; obj++) {
                t = numeric_limits<int>::max();
                if (obj->intersect(ray, t) && t < min_t) {
                    min_t = t;
                    min_obj = *obj;
                }
            }

            // There is an object the ray intersects with
            if (min_t != INFINITY) {
                const Vec p = ray.origin + ray.direction*min_t,
                          n = min_obj.normal(p),
                          l = (objects.back().center - p).normalize();
                const double dt = dot(l, n);

                color = (min_obj.color + objects.back().color*dt) * .5;  // Set color
            }

            // Output color at current pixel
            color_clamp(color);
            out << (int)color.x << ' ' << (int)color.y << ' ' << (int)color.z << "\n";
        }
    }
}
