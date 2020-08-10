#include <fstream>
#include <vector>
#include <cmath>

#define MAX_DEPTH 5
#define FOV 30

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
    bool reflective;
    float ior;  // Index of refraction equal to 0 means non-refractive

    Sphere() { center = Vec(); color = Vec(); radius = 0; radius_sqr = 0; reflective = false; ior = 0.f; }
    Sphere(const Vec& c, const Vec& col, float r, bool ref, float ior) : center(c), color(col), radius(r), radius_sqr(r*r), reflective(ref), ior(ior) {}

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

// Sphere with light output
struct PointLight : Sphere {
    float brightness;

    PointLight(const Vec& c, const Vec& col, float r, float b) : Sphere(c, col, r, false, 0.f), brightness(b) {}
};

// Return the initial ray center at the given (x, y) coordinates
// Use perspective to determine the direction
Ray get_initial_ray(int x, int y, int width, int height) {
    const float invWidth = 1 / float(width), invHeight = 1 / float(height),
                aspectratio = width / float(height),
                angle = tan(M_PI * .5 * FOV / 180.),
                xdir = (2 * ((x + .5) * invWidth) - 1) * angle * aspectratio,
                ydir = (2 * ((y + .5) * invHeight) - 1) * angle;
    return Ray(Vec(x,y,0), Vec(xdir,ydir,1).normalize());
}

// Return the way used to trace reflection
Ray inline get_reflection_ray(const Vec dir, const Vec p, const Vec n) {
    Vec reflection_dir = (dir -  n * 2 * dot(dir, n)).normalize();
    return Ray(p, reflection_dir);
}

// Return the ray used to trace refraction
Ray inline get_refraction_ray(const Vec dir, const Vec p, const Vec n, float ior) {
    const float eta = 1/ior,
                cosi = -dot(n, dir),
                k = 1 - eta*eta * (1 - cosi*cosi);
    Vec refraction_dir = dir * eta + n * (eta * cosi - sqrt(k));

    return Ray(p, refraction_dir);
}

// Return the closest interecting object with r
Sphere min_intersect(const vector<Sphere>& objs, const Ray& r, float& min_t, const Sphere ignore) {
    Sphere min_obj;
    float t = INFINITY;

    for (auto obj = objs.begin(); obj < objs.end(); obj++) {
        if (obj->intersect(r, t) && t < min_t && !(*obj == ignore)) {
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

// Return transmission of reflection obtained by Fresnel's equation
float fresnel(const Vec& dir, const Vec& n, float ior) {
    float cosi = -dot(dir, n),
          etai = 1, etat = ior;

    if (cosi > 0)
        swap(etai, etat);

    const float sint = etai / etat * sqrtf(max(0.f, 1 - cosi*cosi));

    if (sint >= 1)
        return 1;
    else {
        const float cost = sqrtf(max(0.f, 1 - sint*sint));
        cosi = fabsf(cosi);
        const float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost)),
                    Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));

        return (Rs*Rs + Rp*Rp)/2;
    }
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

// Return the calculated color for the given ray using ray tracing
Vec raytrace(const Ray& ray, const vector<Sphere> objects, const vector<Sphere> lights, const int depth, const Sphere ignore){
     // Initialize variables
    Vec color(0,0,0);
    float min_t = INFINITY;
    Sphere min_obj = min_intersect(objects, ray, min_t, ignore);

    // There is an object the ray intersects with
    if (min_t != INFINITY) {
        // Calculate for every light source
        for (auto light = lights.begin(); light < lights.end(); light++) {
            color_clamp(color);
            const Vec p = ray.origin + ray.direction*min_t,  // Position of intersection
                      n = min_obj.normal(p),  // Normal at the intersection point
                      l = (light->center - p).normalize();  // Direction to the light
            const float dt = dot(n, l);

            // Handle reflection and refraction
            if ((min_obj.reflective || min_obj.ior != 0.f) && depth < MAX_DEPTH) {
                Vec reflection_color, refraction_color;

                if (min_obj.reflective) {
                    const Ray reflection_ray = get_reflection_ray(ray.direction, p, n);
                    reflection_color = raytrace(reflection_ray, objects, lights, depth + 1, min_obj);
                }

                if (min_obj.ior != 0.f) {
                    const Ray refraction_ray = get_refraction_ray(ray.direction, p, n, min_obj.ior);
                    refraction_color = raytrace(refraction_ray, objects, lights, depth + 1, min_obj);
                }

                // Combine reflection and refraction
                const float kr = fresnel(ray.direction, n, min_obj.ior);
                if (min_obj.reflective && min_obj.ior != 0.f) {
                    color = reflection_color * kr;
                    color_average(color, refraction_color * (1 - kr));
                } else if (min_obj.ior == 0.f)  // Only reflective
                    color = reflection_color * kr;
                else  // Only refractive
                    color = refraction_color * (1 - kr);
                color_average(color, (min_obj.color + (light->color*dt)) * .15);  // Mix in original color
                return color;
            }

            // Check for and handle shadows
            if (!shadow(min_obj, objects, Ray(p, l), min_t)) {
                Vec new_color = (min_obj.color + (light->color*dt)) * .5;
                color_clamp(new_color);
                color_average(color, new_color);
            } else {
                color_average(color, min_obj.color * .075);
            }
        }
    }

    return color;
}

int main() {
    const int WIDTH = 600, HEIGHT = 600,  // Image width and height
              SIZE = (WIDTH > HEIGHT) ? WIDTH : HEIGHT;
    const Vec WHITE(230,230,230),  // Quick reference colors
              BLACK(0,0,0),
              RED(230,0,0),
              GREEN(0,230,0),
              BLUE(0,0,230);

    // Header for .ppm
    ofstream out = ofstream("result.ppm");
    out << "P3\n" << WIDTH << ' ' << HEIGHT << " 255\n";

    // Scene creation
    vector<Sphere> objects;  // Assume every obect in the scene is a sphere for now
    objects.push_back(Sphere(Vec(.5*SIZE, .5*SIZE, .425*SIZE),  GREEN, .35*SIZE, true,  0.f));
    objects.push_back(Sphere(Vec(.55*SIZE, .3*SIZE, .04*SIZE),  RED,   .15*SIZE, false, 0.f));
    objects.push_back(Sphere(Vec(.35*SIZE, .7*SIZE, .075*SIZE), BLUE,  .12*SIZE, false, 0.f));
    objects.push_back(Sphere(Vec(.5*SIZE, 10001*SIZE, 0),    WHITE,    10000*WIDTH, false, 0.f));  // Making a flat surface is too much effort

    vector<Sphere> lights;  // Point light sources
    lights.push_back(PointLight(Vec(WIDTH, .4*HEIGHT, -.2), WHITE, 1, 1.0));

	Vec color;
    Ray ray;

    // Render each pixel
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            ray = get_initial_ray(x,y,WIDTH,HEIGHT);
            color = raytrace(ray,objects,lights,0,lights.back());

            // Output color at current pixel
            color_clamp(color);
            out << (int)color.x << ' ' << (int)color.y << ' ' << (int)color.z << "\n";
        }
    }

    out.close();
}
