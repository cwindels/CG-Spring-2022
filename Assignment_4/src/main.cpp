////////////////////////////////////////////////////////////////////////////////
// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <fstream>
#include <algorithm>
#include <numeric>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

////////////////////////////////////////////////////////////////////////////////
// Class to store tree
////////////////////////////////////////////////////////////////////////////////
class AABBTree
{
public:
    class Node
    {
    public:
        AlignedBox3d bbox;
        int parent;   // Index of the parent node (-1 for root)
        int left;     // Index of the left child (-1 for a leaf)
        int right;    // Index of the right child (-1 for a leaf)
        int triangle; // Index of the node triangle (-1 for internal nodes)
    };

    std::vector<Node> nodes;
    int root;

    AABBTree() = default;                           // Default empty constructor
    AABBTree(const MatrixXd &V, const MatrixXi &F); // Build a BVH from an existing mesh
};

////////////////////////////////////////////////////////////////////////////////
// Scene setup, global variables
////////////////////////////////////////////////////////////////////////////////
const std::string data_dir = DATA_DIR;
const std::string filename("raytrace.png");
const std::string mesh_filename(data_dir + "bunny.off");

//Camera settings
const double focal_length = 2;
const double field_of_view = 0.7854; //45 degrees
const bool is_perspective = true;
const Vector3d camera_position(0, 0, 2);

//Maximum number of recursive calls
const int max_bounce = 5;

// Triangle Mesh
MatrixXd vertices; // n x 3 matrix (n points)
MatrixXi facets;   // m x 3 matrix (m triangles)
AABBTree bvh;

// Objects
std::vector<Vector3d> sphere_centers;
std::vector<double> sphere_radii;
std::vector<Matrix3d> parallelograms;

//Material for the object, same material for all objects
const Vector4d obj_ambient_color(0.0, 0.5, 0.0, 0);
const Vector4d obj_diffuse_color(0.5, 0.5, 0.5, 0);
const Vector4d obj_specular_color(0.2, 0.2, 0.2, 0);
const double obj_specular_exponent = 256.0;
const Vector4d obj_reflection_color(0.7, 0.7, 0.7, 0);

// Precomputed (or otherwise) gradient vectors at each grid node
const int grid_size = 20;
std::vector<std::vector<Vector2d>> grid;

//Lights
std::vector<Vector3d> light_positions;
std::vector<Vector4d> light_colors;
//Ambient light
const Vector4d ambient_light(0.2, 0.2, 0.2, 0);

//Fills the different arrays
void setup_scene()
{
    //Loads file
    std::ifstream in(mesh_filename);
    std::string token;
    in >> token;
    int nv, nf, ne;
    in >> nv >> nf >> ne;
    vertices.resize(nv, 3);
    facets.resize(nf, 3);
    for (int i = 0; i < nv; ++i)
    {
        in >> vertices(i, 0) >> vertices(i, 1) >> vertices(i, 2);
    }
    for (int i = 0; i < nf; ++i)
    {
        int s;
        in >> s >> facets(i, 0) >> facets(i, 1) >> facets(i, 2);
        assert(s == 3);
    }

    //setup tree
    bvh = AABBTree(vertices, facets);

    grid.resize(grid_size + 1);
    for (int i = 0; i < grid_size + 1; ++i)
    {
        grid[i].resize(grid_size + 1);
        for (int j = 0; j < grid_size + 1; ++j)
            grid[i][j] = Vector2d::Random().normalized();
    }

    //Spheres
    sphere_centers.emplace_back(10, 0, 1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(7, 0.05, -1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(4, 0.1, 1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(1, 0.2, -1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(-2, 0.4, 1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(-5, 0.8, -1);
    sphere_radii.emplace_back(1);

    sphere_centers.emplace_back(-8, 1.6, 1);
    sphere_radii.emplace_back(1);

    //parallelograms
    parallelograms.emplace_back();
    parallelograms.back() << -100, 100, -100,
        -1.25, 0, -1.2,
        -100, -100, 100;

    //Lights
    light_positions.emplace_back(8, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(6, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(4, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(2, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(0, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(-2, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(-4, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0); 
}

////////////////////////////////////////////////////////////////////////////////
// BVH Code
////////////////////////////////////////////////////////////////////////////////

AlignedBox3d bbox_from_triangle(const Vector3d &a, const Vector3d &b, const Vector3d &c)
{
    AlignedBox3d box;
    box.extend(a);
    box.extend(b);
    box.extend(c);
    return box;
}

// returns index 
/*
int build_tree(int i, int j, Node n, std::vector<Node> nodes) // i and j are indices of triangle barycenters, we take a subset of triangles from i to j
{
    if (j - i = 1) // in a leaf
    {
        // append node, nodes.append(n);
        // set n.triangle = i
        // set n.parent, left and right
        // n.bbox = bbox_from_triangle(a,b,c);
        // return box_id
    }

    //find largest axis
    AlignedBox3d barybox;
    for(int k = i; k < j; k++) // looping over triangles from i to j
    {
        //barybox.extend(trianglebarycenterarray[k]);
    }
}
*/
AABBTree::AABBTree(const MatrixXd &V, const MatrixXi &F)
{
    // Compute the centroids of all the triangles in the input mesh
    MatrixXd centroids(F.rows(), V.cols());
    centroids.setZero();
    for (int i = 0; i < F.rows(); ++i)
    {
        for (int k = 0; k < F.cols(); ++k)
        {
            centroids.row(i) += V.row(F(i, k));
        }
        centroids.row(i) /= F.cols();
    }

    // TODO

    // Top-down approach.
    // Split each set of primitives into 2 sets of roughly equal size,
    // based on sorting the centroids along one direction or another.
    // returns index
    
}

////////////////////////////////////////////////////////////////////////////////
// Intersection code
////////////////////////////////////////////////////////////////////////////////

Vector3d ray_triangle_parameters(Vector3d u, Vector3d v, Vector3d d, Vector3d ray_o, Vector3d t_o)
{
    Matrix3d A;
    A << u, v, -d;

    Vector3d k = ray_o + (-1)*t_o;
    Vector3d uvt = A.partialPivLu().solve(k);

    return uvt;
}

double ray_triangle_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, const Vector3d &a, const Vector3d &b, const Vector3d &c, Vector3d &p, Vector3d &N)
{
    // Compute whether the ray intersects the given triangle.
    // If you have done the parallelogram case, this should be very similar to it.
    const Vector3d t_u = b - a;
    const Vector3d t_v = c - a;
    double t = -1;
    
    Vector3d uvt = ray_triangle_parameters(t_u, t_v, ray_direction, ray_origin, a);
    if (0 <= uvt(0) && 0 <= uvt(1) && uvt(2) > 0 && (uvt(0) + uvt(1)) <= 1) {
        // The ray hit the triangle, compute the exact intersection point
        Vector3d ray_intersection = ray_origin + uvt(2)*ray_direction;

        // Compute normal at the intersection point
        Vector3d ray_normal = t_u.cross(t_v);

        p = ray_intersection;
        N = ray_normal.normalized();
        t = uvt(2);
    }

    return t;
}

bool ray_box_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, const AlignedBox3d &box)
{
    // Compute whether the ray intersects the given box.
    // we are not testing with the real surface here anyway.
    double xmin = box.min().x();
    double xmax = box.max().x();
    double ymin = box.min().y();
    double ymax = box.max().y();

    double txmin;
    double txmax;
    double tymin;
    double tymax;

    if(ray_direction(0) >= 0)
    {
        txmin = (xmin - ray_origin(0))/ray_direction(0);
        txmax = (xmax - ray_origin(0))/ray_direction(0);
    } else
    {
        txmin = (xmax - ray_origin(0))/ray_direction(0);
        txmax = (xmin - ray_origin(0))/ray_direction(0);
    }
    if(ray_direction(1) >= 0)
    {
        tymin = (ymin - ray_origin(1))/ray_direction(1);
        tymax = (ymax - ray_origin(1))/ray_direction(1);
    } else
    {
        tymin = (ymin - ray_origin(1))/ray_direction(1);
        tymax = (ymax - ray_origin(1))/ray_direction(1);
    }

    if ((txmin > tymax) || (tymin > txmax))
        return false;
    else
        return true;
}

//Compute the intersection between a ray and a sphere, return -1 if no intersection
double ray_sphere_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, int index, Vector3d &p, Vector3d &N)
{
    const Vector3d sphere_center = sphere_centers[index];
    const double sphere_radius = sphere_radii[index];

    double t = -1;

    double AA = ray_direction.dot(ray_direction);
    double BB = 2 * ray_direction.dot(ray_origin - sphere_center);
    double CC = (ray_origin - sphere_center).dot(ray_origin - sphere_center) - (sphere_radius * sphere_radius);

    double delta = (BB * BB) - (4 * AA * CC);

    if(delta >= 0) {
        double t1 = (-BB + sqrt(delta)) / (2 * AA);
        double t2 = (-BB - sqrt(delta)) / (2 * AA);

        double tmpt = std::min(t1,t2); // Find closest to camera
        if(tmpt < 0)
            tmpt = std::max(t1,t2);
        if(tmpt >= 0)
            t = tmpt;
            p = ray_origin + t * ray_direction;
            N = (p - sphere_center).normalized();
    }

    return t;
}

Vector3d ray_pgram_parameters(Vector3d u, Vector3d v, Vector3d d, Vector3d ray_o, Vector3d p_o)
{
    Matrix3d A;
    A << u, v, -d;

    Vector3d k = ray_o + (-1)*p_o;
    Vector3d uvt = A.partialPivLu().solve(k);

    return uvt;
}

//Compute the intersection between a ray and a paralleogram, return -1 if no intersection
double ray_parallelogram_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, int index, Vector3d &p, Vector3d &N)
{
    const Vector3d pgram_origin = parallelograms[index].col(0);
    const Vector3d A = parallelograms[index].col(1);
    const Vector3d B = parallelograms[index].col(2);
    const Vector3d pgram_u = A - pgram_origin;
    const Vector3d pgram_v = B - pgram_origin;

    double t = -1;

    Vector3d uvt = ray_pgram_parameters(pgram_u, pgram_v, ray_direction, ray_origin, pgram_origin);
    if (0 <= uvt(0) && uvt(0) <= 1 && 0 <= uvt(1) && uvt(1) <= 1 && uvt(2) > 0) {
        // The ray hit the parallelogram, compute the exact intersection point
        Vector3d ray_intersection= ray_origin + uvt(2)*ray_direction;

        // Compute normal at the intersection point
        Vector3d ray_normal = pgram_u.cross(pgram_v);

        p = ray_intersection;
        N = -1 * ray_normal.normalized();
        t = uvt(2);
    }

    return t;
}

//Finds the closest intersecting object returns its index/ is this the nearest object
//In case of intersection it writes into p and N (intersection point and normals)
bool find_nearest_object(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &p, Vector3d &N)
{
    Vector3d tmp_p, tmp_N;
    bool object_hit = false;
    // TODO
    // Method (1): Traverse every triangle and return the closest hit.
    double closest_t = std::numeric_limits<double>::max(); //closest t is "+ infinity"

    for (int i = 0; i < facets.rows(); i++) {
        double aindex = facets.coeff(i,0);
        double bindex = facets.coeff(i,1);
        double cindex = facets.coeff(i,2);
       
        Vector3d a(vertices.coeff(aindex, 0), vertices.coeff(aindex, 1), vertices.coeff(aindex, 2));
        Vector3d b(vertices.coeff(bindex, 0), vertices.coeff(bindex, 1), vertices.coeff(bindex, 2));
        Vector3d c(vertices.coeff(cindex, 0), vertices.coeff(cindex, 1), vertices.coeff(cindex, 2));
       
        const double t = ray_triangle_intersection(ray_origin, ray_direction, a, b, c, tmp_p, tmp_N);
        if (t >= 0)
        {
            object_hit = true;
            //The point is before our current closest t
            if (t < closest_t)
            {
                closest_t = t;
                p = tmp_p;
                N = tmp_N;
            }
        }
    } 
    for (int i = 0; i < sphere_centers.size(); ++i)
    {
        //returns t and writes on tmp_p and tmp_N
        const double t = ray_sphere_intersection(ray_origin, ray_direction, i, tmp_p, tmp_N);
        //We have intersection
        if (t >= 0)
        {
            object_hit = true;
            //The point is before our current closest t
            if (t < closest_t)
            {
                closest_t = t;
                p = tmp_p;
                N = tmp_N;
            }
        }
    }

    for (int i = 0; i < parallelograms.size(); ++i)
    {
        //returns t and writes on tmp_p and tmp_N
        const double t = ray_parallelogram_intersection(ray_origin, ray_direction, i, tmp_p, tmp_N);
        //We have intersection
        if (t >= 0)
        {
            object_hit = true;
            //The point is before our current closest t
            if (t < closest_t)
            {
                closest_t = t;
                p = tmp_p;
                N = tmp_N;
            }
        }
    }

    // Method (2): Traverse the BVH tree and test the intersection with a
    // triangles at the leaf nodes that intersects the input ray.

    return object_hit;
}



////////////////////////////////////////////////////////////////////////////////
// Raytracer code
////////////////////////////////////////////////////////////////////////////////
//Checks if the light is visible
bool is_light_visible(const Vector3d &ray_origin, const Vector3d &ray_direction, const Vector3d &light_position)
{   
    Vector3d p, N;

    bool hit = find_nearest_object(ray_origin, ray_direction, p, N);

    return !hit;
} 

Vector4d shoot_ray(const Vector3d &ray_origin, const Vector3d &ray_direction, int max_bounce)
{
    //Intersection point and normal, these are output of find_nearest_object
    Vector3d p, N;

    const bool nearest_object = find_nearest_object(ray_origin, ray_direction, p, N);

    if (!nearest_object)
    {
        // Return a transparent color
        return Vector4d(0, 0, 0, 0);
    }

    // Ambient light contribution
    const Vector4d ambient_color = obj_ambient_color.array() * ambient_light.array();

    // Punctual lights contribution (direct lighting)
    Vector4d lights_color(0, 0, 0, 0);
    for (int i = 0; i < light_positions.size(); ++i)
    {
        const Vector3d &light_position = light_positions[i];
        const Vector4d &light_color = light_colors[i];

        const Vector3d Li = (light_position - p).normalized();
        const Vector3d Vi = (camera_position - p).normalized();
        const Vector3d Hi = (Li - ray_direction).normalized();

        // Shoot a shadow ray to determine if the light should affect the intersection point and call is_light_visible
        
        Vector3d epsilon = 0.000001 * Li;
        bool visible = is_light_visible(p + epsilon, Li, light_position);
        if(!visible)
            continue; 
        Vector4d diff_color = obj_diffuse_color;
        
        if (nearest_object == 4)
        {
            //Compute UV coodinates for the point on the sphere
            const double x = p(0) - sphere_centers[nearest_object][0];
            const double y = p(1) - sphere_centers[nearest_object][1];
            const double z = p(2) - sphere_centers[nearest_object][2];
            const double tu = acos(z / sphere_radii[nearest_object]) / 3.1415;
            const double tv = (3.1415 + atan2(y, x)) / (2 * 3.1415);
        }

        // Diffuse contribution
        const Vector4d diffuse = diff_color * std::max(Li.dot(N), 0.0);

        // Specular contribution
        const Vector4d specular = obj_specular_color * std::pow(std::max(N.dot(Hi), 0.0), obj_specular_exponent);

        // Attenuate lights according to the squared distance to the lights
        const Vector3d D = light_position - p;
        lights_color += (diffuse + specular).cwiseProduct(light_color) / D.squaredNorm();
    }

    Vector4d refl_color = obj_reflection_color;
    if (nearest_object)
    {
        refl_color = Vector4d(0.5, 0.5, 0.5, 0);
    }
    // Compute the color of the reflected ray and add its contribution to the current point color.
    // use refl_color
    Vector3d r = ray_direction - (N.dot(ray_direction)) * 1.99 * N;
    Vector4d reflection_color = refl_color.cwiseProduct(shoot_ray(p + 0.0001 * r, r, max_bounce--));

    // Rendering equation
    Vector4d C = ambient_color + lights_color + reflection_color;

    //Set alpha to 1
    C(3) = 1;

    return C;
}

////////////////////////////////////////////////////////////////////////////////

void raytrace_scene()
{
    std::cout << "Simple ray tracer." << std::endl;

    int w = 640;
    int h = 480;
    MatrixXd R = MatrixXd::Zero(w, h);
    MatrixXd G = MatrixXd::Zero(w, h);
    MatrixXd B = MatrixXd::Zero(w, h);
    MatrixXd A = MatrixXd::Zero(w, h); // Store the alpha mask

    // The camera always points in the direction -z
    // The sensor grid is at a distance 'focal_length' from the camera center,
    // and covers an viewing angle given by 'field_of_view'.
    double aspect_ratio = double(w) / double(h);
    //TODO
    double image_y = focal_length * tan(field_of_view/2);;
    double image_x = aspect_ratio * image_y;;

    // The pixel grid through which we shoot rays is at a distance 'focal_length'
    const Vector3d image_origin(-image_x, image_y, camera_position[2] - focal_length);
    const Vector3d x_displacement(2.0 / w * image_x, 0, 0);
    const Vector3d y_displacement(0, -2.0 / h * image_y, 0);

    for (unsigned i = 0; i < w; ++i)
    {
        for (unsigned j = 0; j < h; ++j)
        {
            const Vector3d pixel_center = image_origin + (i + 0.5) * x_displacement + (j + 0.5) * y_displacement;

            // Prepare the ray
            Vector3d ray_origin;
            Vector3d ray_direction;

            if (is_perspective)
            {
                // Perspective camera
                ray_origin = camera_position;
                ray_direction = (pixel_center - camera_position).normalized();
            }
            else
            {
                // Orthographic camera
                ray_origin = pixel_center;
                ray_direction = Vector3d(0, 0, -1);
            }

            const Vector4d C = shoot_ray(ray_origin, ray_direction, max_bounce);
            R(i, j) = C(0);
            G(i, j) = C(1);
            B(i, j) = C(2);
            A(i, j) = C(3);
        }
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    setup_scene();

    raytrace_scene();
    return 0;
}
