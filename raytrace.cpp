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

    bool operator!=(const Vec& v) const {
        return x != v.x || y != v.y || z != v.z;
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
        const float length = sqrt(x*x + y*y + z*z);
        return Vec(x/length, y/length, z/length);
    }
};

// Return the dot product of two vectors
float inline dot(const Vec& a, const Vec& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

// Ray composed of an origin point and a direction
struct Ray {
    Vec origin, direction;

    Ray() { origin = Vec(); direction = Vec(); }
    Ray(const Vec& o, const Vec& d) : origin(o), direction(d) {};
};

// Represent a sphere
// Used in scene creation
struct Sphere {
    Vec center, color;
    float radius, radius_sqr;

    Sphere() { center = Vec(); color = Vec(); radius = 0; radius_sqr = 0; }
    Sphere(const Vec& c, const Vec& col, float r) : center(c), color(col), radius(r), radius_sqr(r*r) {}

    // Return true if ray intersects with the sphere
    // If true, t will hold the minimum intersecting t-value
    bool intersect(const Ray& ray, float& t) const {
        const Vec l = center - ray.origin;
        const float tca = dot(l, ray.direction);
        if (tca < 0)
            return false;

        const float d2 = dot(l, l) - tca*tca;
        if (d2 > radius_sqr)
            return false;

        const float thc = sqrt(radius_sqr - d2),
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
                angle = tan(M_PI * .5 * fov / 180.),
                xdir = (2 * ((x + .5) * invWidth) - 1) * angle * aspectratio,
                ydir = (2 * ((y + .5) * invHeight) - 1) * angle;
    return Ray(Vec(x,y,0), Vec(xdir,ydir,1).normalize());
}

// Return the closest interecting object with r
Sphere min_intersect(const vector<Sphere>& objs, const Ray& r, float& min_t) {
    Sphere min_obj;
    float t = INFINITY;

    for (auto obj = objs.begin(); obj < objs.end(); obj++) {
        if (obj->intersect(r, t) && t < min_t) {
            min_t = t;
            min_obj = *obj;
        }
    }

    return min_obj;
}

// Return true if there is an intersecting object between s and the light source
// I.e. if s should have a shadow cast on it
bool shadow(const Sphere& s, const vector<Sphere>& objs, const Ray& r, float min_t) {
    float t = INFINITY;
    for (auto obj = objs.begin(); obj < objs.end(); obj++) {
        if (*obj == s)
            continue;

        if (obj->intersect(r, t) && t < min_t)
            return true;
    }

    return false;
}

// Average two color vectors into the first
void inline color_average(Vec& v, const Vec w) {
    if (v != Vec(0,0,0)) {  // Don't bother averaging against black
        v.x = sqrt(v.x*v.x + w.x*w.x);
        v.y = sqrt(v.y*v.y + w.y*w.y);
        v.z = sqrt(v.z*v.z + w.z*w.z);
    } else
        v = w;
}

// Clamp a color vector to 8bit colors
// (Nothing below 0 or above 255)
void inline color_clamp(Vec& v) {
    v.x = (v.x > 255) ? 255 : (v.x < 0) ? 0 : v.x;
    v.y = (v.y > 255) ? 255 : (v.y < 0) ? 0 : v.y;
    v.z = (v.z > 255) ? 255 : (v.z < 0) ? 0 : v.z;
}

int main() {
    const int WIDTH = 3840, HEIGHT = 2160;  // Image width and height
    const Vec WHITE(230,230,230),  // Quick reference colors
              BLACK(0,0,0),
              RED(230,0,0),
              GREEN(0,230,0),
              BLUE(0,0,230),
              YELLOW(230,230,0),
              PINK(230,0,230);

    // Header for .ppm
    ofstream out = ofstream("result.ppm");
    out << "P3\n" << WIDTH << ' ' << HEIGHT << " 255\n";

    // Scene creation
    vector<Sphere> objects;  // Assume every obect in the scene is a sphere for now
    objects.push_back(Sphere(Vec(.5*WIDTH, .5*HEIGHT, .425*WIDTH), GREEN, .35*WIDTH));
    objects.push_back(Sphere(Vec(.55*WIDTH, .3*HEIGHT, .04*WIDTH), RED, .15*WIDTH));
    objects.push_back(Sphere(Vec(.35*WIDTH, .7*HEIGHT, .075*WIDTH), BLUE, .12*WIDTH));
    objects.push_back(Sphere(Vec(.5*WIDTH, 10001*WIDTH, 0), WHITE, 10000*WIDTH));  // Making a flat surface is too much effort

    vector<Sphere> lights;  // Point light sources
    lights.push_back(Sphere(Vec(1*WIDTH, .4*HEIGHT, -.2), PINK, 1));
    lights.push_back(Sphere(Vec(0, .4*HEIGHT, -.2), YELLOW, 1));

	Vec color, new_color;
    Ray ray;
	float min_t;

    // Render each pixel
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            // Initialize variables
            color = BLACK;
            ray = get_initial_ray(x,y,WIDTH,HEIGHT);
            min_t = INFINITY;
            Sphere min_obj = min_intersect(objects, ray, min_t);

            // There is an object the ray intersects with
            if (min_t != INFINITY) {
                // Calculate for every light source
                for (auto light = lights.begin(); light < lights.end(); light++) {
                    color_clamp(color);
                    const Vec p = ray.origin + ray.direction*min_t,  // Position of intersection
                              n = min_obj.normal(p),  // Normal at the intersection point
                              l = (light->center - p).normalize();  // Direction to the light
                    const float dt = dot(n, l);

                    // Check for and handle shadows
                    if (!shadow(min_obj, objects, Ray(p, l), min_t)) {
                        new_color = (min_obj.color + (light->color*dt)) * .5;
                        color_clamp(new_color);
                        color_average(color, new_color);
                    } else {
                        color_average(color, min_obj.color * .075);
                    }
                }
            }

            // Output color at current pixel
            color_clamp(color);
            out << (int)color.x << ' ' << (int)color.y << ' ' << (int)color.z << "\n";
        }
    }

    out.close();
}
