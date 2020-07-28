#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

// Vector in 3 dimensions
struct Vec {
    float x, y, z;

    Vec() : x(0), y(0), z(0) {};
    Vec(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {};

    // Vector operations
    Vec operator+(const Vec& v) const {
        return Vec(x + v.x, y + v.y, z + v.z);
    }
    Vec operator-(const Vec& v) const {
        return Vec(x - v.x, y - v.y, z - v.z);
    }

    bool operator==(const Vec& v) const {
        return x == v.x && y == v.y && z == v.z;
    }

    // Scalar operations
    Vec operator*(float c) const {
        return Vec(c * x, c * y, c * z);
    }
    Vec operator/(float c) const {
        return Vec(x / c, y / c, z / c);
    }

    // Return a normalized unit vector (magnitude of 1)
    Vec normalize() const {
        float length = sqrt(x*x + y*y + z*z);
        return Vec(x/length, y/length, z/length);
    }
};

// Return the dot product of two vectors
float dot(const Vec& a, const Vec& b) {
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
    float radius;

    Sphere() { center = Vec(); color = Vec(); radius = 0; }
    Sphere(const Vec& c, const Vec& col, float r) : center(c), color(col), radius(r) {}

    // Return true if ray intersects with the sphere
    // If true, t will hold the minimum intersecting t-value
    bool intersect(const Ray& ray, float &t) const {
        const Vec l = center - ray.origin;
        const float tca = dot(l, ray.direction);
        if (tca < 0)
            return false;

        const float d2 = dot(l, l) - tca*tca;
        if (d2 > radius*radius)
            return false;

        const float thc = sqrt(radius*radius - d2),
                    t0 = tca - thc,
                    t1 = tca + thc;

        t = (t0 < t1) ? t0 : t1;
        return true;
    }

    // Return the vector normal to the circle and p
    Vec normal(const Vec& p) const {
        return (p - center).normalize();
    }

    bool operator==(const Sphere& s) const {
        return center == s.center && color == s.color && radius == s.radius;
    }
};

// Return the initial ray center at the given (x, y) coordinates
// Use perspective to determine the direction
Ray get_initial_ray(int x, int y, int width, int height, int fov=30) {
    const float invWidth = 1 / float(width), invHeight = 1 / float(height),
                aspectratio = width / float(height),
                angle = tan(M_PI * 0.5 * fov / 180.),
                xdir = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio,
                ydir = (2 * ((y + 0.5) * invHeight) - 1) * angle;
    return Ray(Vec(x,y,0), Vec(xdir,ydir,1).normalize());
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
    const Vec white(230,230,230),  // Quick reference colors
              whitish(255,242,226),
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
    objects.push_back(Sphere(Vec(.5*WIDTH, .5*HEIGHT, 0.425*WIDTH), green, .35*WIDTH));  // Sphere centered and with diameter half as wide as the screen
    objects.push_back(Sphere(Vec(.55*WIDTH, .3*HEIGHT, .04*WIDTH), red, .15*WIDTH));
    objects.push_back(Sphere(Vec(.35*WIDTH, .7*HEIGHT, .075*WIDTH), blue, .12*WIDTH));
    objects.push_back(Sphere(Vec(.5*WIDTH, 10001*WIDTH, 0), white, 10000*WIDTH));  // Making a flat surface is too much effort
    objects.push_back(Sphere(Vec(1.*WIDTH, .4*HEIGHT, -.2), whitish, 1));  // Point light source at the same depth as the sphere

	Vec color;
    Ray ray;
	float t, min_t;

    // Render each pixel
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            // Initialize variables
            color = black;
            ray = get_initial_ray(x,y,WIDTH,HEIGHT);
            min_t = INFINITY;
            Sphere min_obj;

            // Find closest interecting object
            for (auto obj = objects.begin(); obj < objects.end() - 1; obj++) {
                t = INFINITY;
                if (obj->intersect(ray, t) && t < min_t) {
                    min_t = t;
                    min_obj = *obj;
                }
            }

            // There is an object the ray intersects with
            if (min_t != INFINITY) {
                const Vec p = ray.origin + ray.direction*min_t,  // Position of intersection
                          n = min_obj.normal(p),  // Normal at the intersection point
                          l = (objects.back().center - p).normalize();  // Direction to the light
                const float dt = dot(n, l);

                bool isShadow = false;
                for (auto obj = objects.begin(); obj < objects.end() - 1 && !isShadow; obj++) {
                    if (*obj == min_obj)
                        continue;

                    t = INFINITY;
                    if (obj->intersect(Ray(p + n * 1e-4, l), t) && t < min_t)
                        isShadow = true;
                }

                if (!isShadow)
                    color = (min_obj.color + (objects.back().color*dt)) * .5;  // Set color
                else
                    color = min_obj.color * .075;
                
            }

            // Output color at current pixel
            color_clamp(color);
            out << (int)color.x << ' ' << (int)color.y << ' ' << (int)color.z << "\n";
        }
    }

    out.close();
}
