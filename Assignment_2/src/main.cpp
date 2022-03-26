// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

void raytrace_sphere()
{
    std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

    const std::string filename("sphere_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800);
    MatrixXd R = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd G = MatrixXd::Zero(800, 800);
    MatrixXd B = MatrixXd::Zero(800, 800);
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the
    // unit square (-1,1) in x and y
    Vector3d image_origin(-1, 1, 1);
    Vector3d x_displacement(2.0 / A.cols(), 0, 0);
    Vector3d y_displacement(0, -2.0 / A.rows(), 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < A.cols(); ++i)
    {
        for (unsigned j = 0; j < A.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            Vector3d ray_origin = pixel_center;
            Vector3d ray_direction = camera_view_direction;

            // Intersect with the sphere
            // NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
            Vector2d ray_on_xy(ray_origin(0), ray_origin(1));
            const double sphere_radius = 0.9;

            if (ray_on_xy.norm() < sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection(
                    ray_on_xy(0), ray_on_xy(1),
                    sqrt(sphere_radius * sphere_radius - ray_on_xy.squaredNorm()));

                // Compute normal at the intersection point
                Vector3d ray_normal = ray_intersection.normalized();

                // Simple diffuse model
                // C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;
                R(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;
                G(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;
                B(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                // C(i, j) = std::max(C(i, j), 0.);
                R(i, j) = std::max(R(i, j), 0.);
                G(i, j) = std::max(G(i, j), 0.);
                B(i, j) = std::max(B(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}
 
Vector3d ray_pgram_parameters(Vector3d u, Vector3d v, Vector3d d, Vector3d ray_o, Vector3d p_o)
{
    Matrix3d A;
    A << u, v, -d;

    Vector3d k = ray_o + (-1)*p_o;
    Vector3d uvt = A.partialPivLu().solve(k);

    return uvt;
    0<=uvt(0) && uvt(0)<=1 && 0<=uvt(1) && uvt(1)<=1 && uvt(2)>0;
}

void raytrace_parallelogram()
{
    std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

    const std::string filename("plane_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(0, 0.7, -10);
    const Vector3d pgram_v(1, 0.4, 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;
            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // Check if the ray intersects with the parallelogram
            
            Vector3d uvt = ray_pgram_parameters(pgram_u, pgram_v, ray_direction, ray_origin, pgram_origin);
            if (0 <= uvt(0) && uvt(0) <= 1 && 0 <= uvt(1) && uvt(1) <= 1 && uvt(2) > 0)
            {
                // The ray hit the parallelogram, compute the exact intersection point
                Vector3d ray_intersection= ray_origin + uvt(2)*ray_direction;
                ray_intersection = ray_intersection.normalized();

                // Compute normal at the intersection point
                Vector3d ray_normal = pgram_u.cross(pgram_v);
                ray_normal = ray_normal.normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_perspective()
{
    std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

    const std::string filename("plane_perspective.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Parameters of the parallelogram (position of the lower-left corner +
    // two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(0, 0.7, -10);
    const Vector3d pgram_v(1, 0.4, 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray (origin point and direction)
            const Vector3d ray_origin = camera_origin;
            const Vector3d ray_direction = pixel_center - camera_origin;

            // Check if the ray intersects with the parallelogram
            Vector3d uvt = ray_pgram_parameters(pgram_u, pgram_v, ray_direction, ray_origin, pgram_origin);
            if (0 <= uvt(0) && uvt(0) <= 1 && 0 <= uvt(1) && uvt(1) <= 1 && uvt(2) > 0)
            {
                // The ray hit the parallelogram, compute the exact intersection
                // point
                Vector3d ray_intersection= ray_origin + uvt(2)*ray_direction;
                ray_intersection = ray_intersection.normalized();

                // Compute normal at the intersection point
                Vector3d ray_normal = pgram_u.cross(pgram_v);
                ray_normal = ray_normal.normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_shading()
{
    std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

    const std::string filename("shading.png");
    MatrixXd R = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd G = MatrixXd::Zero(800, 800);
    MatrixXd B = MatrixXd::Zero(800, 800);
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the
    // unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / A.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / A.rows(), 0);

    //Sphere setup
    const Vector3d sphere_center(0, 0, 0);
    const double sphere_radius = 0.9;

    //material params
    const Vector3d diffuse_color(1, 0, 1);
    const double specular_exponent = 100;
    const Vector3d specular_color(0., 0, 1);

    // Single light source
    const Vector3d light_position(-1, 1, 1);
    double ambient = 0.1;

    for (unsigned i = 0; i < A.cols(); ++i)
    {
        for (unsigned j = 0; j < A.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            //
            double AA = ray_direction.dot(ray_direction);
            double BB = 2 * ray_origin.dot(ray_direction) + ray_direction.dot(sphere_center);
            double CC = (sphere_center-ray_origin).dot(sphere_center-ray_origin) - sphere_radius * sphere_radius;

            double delta = BB * BB - 4 * AA * CC;

            if(delta<0)
                continue;

            double t1 = (-BB + sqrt(delta)) / (2 * AA);
            double t2 = (-BB - sqrt(delta)) / (2 * AA);

            double t = std::min(t1,t2); // Find closest to camera
            if(t<0)
                t = std::max(t1,t2);

            if(t<0)
                continue;
            //
            // Intersect with the sphere
            // NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
            Vector2d ray_on_xy(ray_origin(0), ray_origin(1));
            const double sphere_radius = 0.9;

            if (ray_on_xy.norm() < sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection = ray_origin + t * ray_direction;
                    //(
                    //ray_on_xy(0), ray_on_xy(1),
                    //sqrt(sphere_radius * sphere_radius - ray_on_xy.squaredNorm()));

                // Compute normal at the intersection point ----------------------------
                Vector3d ray_normal = (ray_intersection - sphere_center).normalized();

                Vector3d Li = (light_position-ray_intersection).normalized();
                Vector3d Vi = (camera_origin - ray_intersection).normalized();
                Vector3d Hi = (Vi+Li).normalized();

                double ndotl = ray_normal.dot(Li);
                if(ndotl < 0)
                    ndotl = 0;
                double ndoth = ray_normal.dot(Hi);
                if(ndoth < 0)
                    ndoth = 0;

                // TODO: Add shading parameter here
                const double diffuse = (light_position - ray_intersection).normalized().dot(ray_normal);
                const double specular = (light_position - ray_intersection).normalized().dot(ray_normal);

                // Simple diffuse model
                // C(i, j) = ambient + diffuse(i, j) + specular(i, j);
                R(i, j) = ambient + diffuse*diffuse_color(0)*ndotl + specular_color(0)*specular*pow(ndoth,specular_exponent);
                G(i, j) = ambient + diffuse*diffuse_color(1)*ndotl + specular_color(1)*specular*pow(ndoth,specular_exponent);
                B(i, j) = ambient + diffuse*diffuse_color(2)*ndotl + specular_color(2)*specular*pow(ndoth,specular_exponent);

                // Clamp to zero
                R(i, j) = std::max(R(i, j), 0.);
                G(i, j) = std::max(G(i, j), 0.);
                B(i, j) = std::max(B(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}

int main()
{
    raytrace_sphere();
    raytrace_parallelogram();
    raytrace_perspective();
    raytrace_shading();

    return 0;
}
