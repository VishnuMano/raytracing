#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

class Vector3 {
public:
    float x, y, z;

    Vector3() : x(0), y(0), z(0) {}
    Vector3(float x, float y, float z) : x(x), y(y), z(z) {}
};

class Ray {
public:
    Vector3 origin;
    Vector3 direction;

    Ray(const Vector3& origin, const Vector3& direction)
        : origin(origin), direction(direction) {}
};

class Sphere {
public:
    Vector3 center;
    float radius;

    Sphere(const Vector3& center, float radius)
        : center(center), radius(radius) {}

    bool intersect(const Ray& ray, float& t) const {
        Vector3 oc = ray.origin - center;
        float a = dot(ray.direction, ray.direction);
        float b = 2.0f * dot(oc, ray.direction);
        float c = dot(oc, oc) - radius * radius;
        float discriminant = b * b - 4 * a * c;
        
        if (discriminant > 0) {
            float t1 = (-b - sqrt(discriminant)) / (2 * a);
            float t2 = (-b + sqrt(discriminant)) / (2 * a);
            
            t = (t1 < t2) ? t1 : t2;
            return true;
        }
        
        return false;
    }

private:
    float dot(const Vector3& a, const Vector3& b) const {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }
};

int main() {
    int image_width = 800;
    int image_height = 400;

    std::vector<Vector3> image(image_width * image_height);

    Vector3 camera_position(0, 0, 0);
    float viewport_height = 2.0;
    float viewport_width = 4.0;
    float focal_length = 1.0;
    
    Vector3 horizontal(viewport_width, 0, 0);
    Vector3 vertical(0, viewport_height, 0);
    Vector3 lower_left_corner = camera_position - horizontal / 2 - vertical / 2 - Vector3(0, 0, focal_length);

    Sphere sphere(Vector3(0, 0, -1), 0.5);

    for (int j = image_height - 1; j >= 0; --j) {
        for (int i = 0; i < image_width; ++i) {
            float u = float(i) / float(image_width - 1);
            float v = float(j) / float(image_height - 1);

            Ray ray(camera_position, lower_left_corner + u * horizontal + v * vertical - camera_position);
            float t;

            if (sphere.intersect(ray, t)) {
                Vector3 hit_point = ray.origin + t * ray.direction;
                Vector3 normal = (hit_point - sphere.center) / sphere.radius;
                Vector3 color = 0.5f * (normal + Vector3(1, 1, 1));
                image[j * image_width + i] = color;
            } else {
                
                Vector3 white(1, 1, 1);
                Vector3 blue(0.5, 0.7, 1.0);
                Vector3 color = (1.0 - v) * white + v * blue;
                image[j * image_width + i] = color;
            }
        }
    }

    std::ofstream image_file("output.ppm");
    image_file << "P3\n" << image_width << " " << image_height << "\n255\n";
    
    for (int i = 0; i < image.size(); ++i) {
        int ir = static_cast<int>(255.99 * image[i].x);
        int ig = static_cast<int>(255.99 * image[i].y);
        int ib = static_cast<int>(255.99 * image[i].z);
        image_file << ir << " " << ig << " " << ib << "\n";
    }

    image_file.close();

    return 0;
}