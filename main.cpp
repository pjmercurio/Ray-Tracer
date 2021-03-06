
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include "Eigen/Dense"
#include "CImg.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif


#include <time.h>
#include <math.h>

using namespace Eigen;
using namespace std;
using namespace cimg_library;


#define HEIGHT 1000
#define WIDTH 1000
#define RECURSION_DEPTH 2
#define EPSILON 0.000001    //Used when determining triangle-ray intersections
inline float sqr(float x) { return x*x; }


//****************************************************
// Global Variables
//****************************************************


// Initialize the CImg file with its dimensions and color
CImg<unsigned char> main_image(WIDTH, HEIGHT, 1, 3, 0);

char *imageName;

struct ray {
    Vector3f origin, distance;
    void normalize(){distance.normalize();}
};
struct camera {
    float ex, ey, ez, llx, lly, llz, lrx, lry, lrz, ulx, uly, ulz, urx, ury, urz;
};
struct sphere {
    float cx, cy, cz, r;
    int mat_num;
    Vector3f getNormal(Vector3f piece){ Vector3f normal = piece - Vector3f(cx, cy, cz); normal.normalize(); return normal;}
    bool equals(sphere other){ return cx == other.cx && cy == other.cx && cz == other.cz && r == other.r;}
};
struct triangle {
    float ax, ay, az, bx, by, bz, cx, cy, cz;
    bool computedNormal = false;
    Vector3f normal;
    int mat_num;
    Vector3f getV1(){ return Vector3f(ax, ay, az);}
    Vector3f getV2(){ return Vector3f(bx, by, bz);}
    Vector3f getV3(){ return Vector3f(cx, cy, cz);}
    Vector3f getNormal(){ if(computedNormal) return normal; else { normal = (getV3() - getV1()).cross(getV2() - getV1()); return normal;}}
    bool equals(triangle other){ return getV1() == other.getV1() && getV2() == other.getV2() && getV3() == other.getV3();}
    void print(){printf("A: %f, %f, %f   B: %f, %f, %f   C: %f, %f, %f", ax, ay, az, bx, by, bz, cx, cy, cz);}
};
struct ambient_light {
    float R, G, B;
};
struct point_light {
    float x, y, z, R, G, B;
    int falloff;
};
struct directional_light {
    float x, y, z, R, G, B;
};
struct material {
    float kar, kag, kab, kdr, kdg, kdb, ksr, ksg, ksb, ksp, krr, krg, krb;
};
struct transformation {
    char type;
    float x, y, z;
};
struct pixel {
    unsigned char R, G, B;
};

pixel calculateTrianglePhongShading(float x, float y, float z, Vector3f R, triangle t, int depth);
pixel calculatePhongShading(float x, float y, float z, Vector3f R, sphere s, int depth);


int camera_count = 1;
int ambient_light_count = 0;
int point_light_count = 0;
int directional_light_count = 0;
int triangle_count = 0;
int sphere_count = 0;
int material_count = 0;
int transformation_count = 0;

Vector3f current_translation(0, 0, 0);
Vector3f current_rotation(0, 0, 0);
Vector3f current_scale(0, 0, 0);
Vector3f ambient_term(0, 0, 0);
Vector3f diffuse_term(0, 0, 0);
Vector3f specular_term(0, 0, 0);


camera the_camera = {0,  0,  0,
    -1,  1, -1,
    1,  1, -1,
    1, -1, -1,
    -1, -1, -1};


ambient_light al_array[10];
point_light pl_array[10];
directional_light dl_array[10];
triangle triangle_array[10];
sphere sphere_array[10];
material material_array[10];
transformation transformation_array[10];
sphere piece_sphere;
triangle piece_tri;
pixel image_buffer[WIDTH][HEIGHT];
string image_name_str = "image";

//*************************************************************//
// Intersection methods for checking intersections with *rays* //
//*************************************************************//

float computeSphereIntersection(sphere s, ray r) {
    float a = r.distance.dot(r.distance), dist_x = r.origin(0) - s.cx, dist_y = r.origin(1) - s.cy, dist_z = r.origin(2) - s.cz;
    float b = 2.0 * r.distance(0) * dist_x + 2.0 * r.distance(1) * dist_y + 2.0 * r.distance(2) * dist_z;
    float c = dist_x * dist_x + dist_y * dist_y + dist_z * dist_z - s.r * s.r;
    float discriminant = b*b - 4*c*a;
    
    if (discriminant < 0)
        return -1;
    else
        return min((-b + sqrt(discriminant))/(2*a), (-b - sqrt(discriminant))/(2*a));
}

float computeTriangleIntersection(triangle tri, ray lightray)
{
    Vector3f V1 = tri.getV1(), V2 = tri.getV2(), V3 = tri.getV3(), origin = lightray.origin, distance = lightray.distance;
    distance.normalize();
    Vector3f e1, e2;  //Edge1, Edge2
    Vector3f P, Q, T;
    float det, u, v;
    float t;
    
    //Find vectors for two edges sharing V1
    e1 =  V2 - V1;
    e2 = V3 - V1;
    //Begin calculating determinant - also used to calculate u parameter
    P = distance.cross(e2);
    //if determinant is near zero, ray lies in plane of triangle
    det = e1.dot(P);
    
    //NOT CULLING
    if(det > -EPSILON && det < EPSILON)
        return -1;//no intersection
    
    //calculate distance from V1 to ray origin
    T = origin - V1;
    
    //Calculate u parameter and test bound
    u = T.dot(P) / det;
    //The intersection lies outside of the triangle
    if(u < 0.f || u > 1.f) return -1;
    
    //Prepare to test v parameter
    Q = T.cross(e1);
    
    //Calculate V parameter and test bound
    v = distance.dot(Q) / det;
    //The intersection lies outside of the triangle
    if(v < 0.f || u + v  > 1.f) return -1;  //no intersection
    
    t = e2.dot(Q)/ det;
    
    if(t > EPSILON)
        return t;
    
    return -1;
}

// Checks if an object is occluded from a point light's light ray to the object's "surface" at (x,y,z)
bool occluded_pl(float x, float y, float z, ray lightray){
    int num_intersected_objects = 0;
    for (int s_index = 0; s_index < sphere_count; s_index++) {
        float t_intersect = computeSphereIntersection(sphere_array[s_index], lightray);
        if (t_intersect >= 0 && t_intersect <= 1)
            num_intersected_objects++;
    }
    for (int t_index = 0; t_index < triangle_count; t_index++) {
        float t_intersect = computeTriangleIntersection(triangle_array[t_index], lightray);
        if (t_intersect >= 0 && t_intersect <= 1)
            num_intersected_objects++;
    }
    
    return num_intersected_objects > 1;
}

Vector3f reflected_color(ray reflection_ray, int reflection_depth_left, Vector3f piece, bool piece_is_tri){
    
    float closest_t = FLT_MAX;
    sphere closest_sphere;
    triangle closest_tri;
    bool closest_is_tri = false;
    
    for (int s_index = 0; s_index < sphere_count; s_index++) {
        float t_intersect = computeSphereIntersection(sphere_array[s_index], reflection_ray);
        if ((piece_is_tri || !piece_sphere.equals(sphere_array[s_index])) && t_intersect >= 0 && t_intersect < closest_t){
            closest_t = t_intersect;
            closest_sphere = sphere_array[s_index];
            closest_is_tri = false;
        }
    }
    for (int t_index = 0; t_index < triangle_count; t_index++) {
        float t_intersect = computeTriangleIntersection(triangle_array[t_index], reflection_ray);
        if ((!piece_is_tri || !triangle_array[t_index].equals(piece_tri)) && t_intersect >= 0 && t_intersect < closest_t){
            closest_t = t_intersect;
            closest_tri = triangle_array[t_index];
            closest_is_tri = true;
        }
    }
    //Found the closest object to reflect
    pixel reflected_pixel;
    Vector3f reflected_point(reflection_ray.origin(0) + reflection_ray.distance(0) * closest_t, reflection_ray.origin(1) + reflection_ray.distance(1) * closest_t, reflection_ray.origin(2) + reflection_ray.distance(2) * closest_t);
    if(!piece_is_tri && piece_sphere.equals(closest_sphere)){
        printf("OH NO SELF REFLECTION!!! SUCH PHYSICS");
        exit(-1);
    }
    if(closest_t == FLT_MAX)
        return Vector3f(0, 0, 0);
    else if(closest_is_tri)
        reflected_pixel = calculateTrianglePhongShading(reflected_point(0), reflected_point(1), reflected_point(2), Vector3f(0, 0, 0), closest_tri, reflection_depth_left);
    else
        reflected_pixel = calculatePhongShading(reflected_point(0), reflected_point(1), reflected_point(2), Vector3f(0, 0, 0), closest_sphere, reflection_depth_left);
    char id, r_id;
    if(piece_is_tri)
        id = 't';
    else
        id = 's';
    if(closest_is_tri)
        r_id = 't';
    else
        r_id = 's';
    //printf("location of reflected %c pixel %f, %f, %f - location of reflection %c pixel %f, %f, %f - color added ", id, reflected_point(0), reflected_point(1), reflected_point(2), r_id, piece(0), piece(1), piece(2));
    
    return Vector3f(reflected_pixel.R, reflected_pixel.G, reflected_pixel.B);
}

pixel calculatePhongShading(float x, float y, float z, Vector3f R, sphere s, int depth) {
    
    if(z > 0)
        exit(-1);
    Vector3f piece(x, y, z);
    Vector3f surfaceNormal = s.getNormal(piece);
    Vector3f view_ray_direction = piece;
    view_ray_direction.normalize();
    piece_sphere = s;
    
    for (int a = 0; a < ambient_light_count; a++){
        R += Vector3f(al_array[a].R * material_array[s.mat_num].kar, al_array[a].G * material_array[s.mat_num].kag, al_array[a].B * material_array[s.mat_num].kab);
    }
    
    if ((material_array[s.mat_num].krr > 0 || material_array[s.mat_num].krg > 0 || material_array[s.mat_num].krb > 0) && depth > 0) {
        Vector3f reflected_direction = 2 * (surfaceNormal.dot(-view_ray_direction) * surfaceNormal + view_ray_direction) - view_ray_direction;
        reflected_direction.normalize();
        ray reflected_ray;
        reflected_ray.origin = piece;
        reflected_ray.distance = reflected_direction;
        Vector3f resolved_reflected_color = reflected_color(reflected_ray, depth - 1, piece, false);
        if (resolved_reflected_color(0) > 0 || resolved_reflected_color(1) > 0 || resolved_reflected_color(2) > 0){
            Vector3f applied_reflection_color(resolved_reflected_color(0) * material_array[s.mat_num].krr, resolved_reflected_color(1) * material_array[s.mat_num].krg, resolved_reflected_color(2) * material_array[s.mat_num].krb);
            R(0) += applied_reflection_color(0);
            R(1) += applied_reflection_color(1);
            R(2) += applied_reflection_color(2);
        }
    }
    
    
    for (int p = 0; p < point_light_count; p++) {
        Vector3f lightloc(pl_array[p].x, pl_array[p].y, pl_array[p].z);
        Vector3f lightcolor;
        ray lightray;
        lightray.origin = lightloc;
        Vector3f distance = lightloc - piece;
        lightray.distance = piece - lightloc;
        
        if(occluded_pl(x, y, z, lightray))
            continue;
        
        float dist = sqrt(distance.dot(distance));
        float dist_squared = dist * dist;
        if(!pl_array[p].falloff)
            lightcolor = Vector3f(pl_array[p].R,pl_array[p].G,pl_array[p].B);
        else if(pl_array[p].falloff == 1)
            lightcolor = Vector3f(pl_array[p].R/dist,pl_array[p].G/dist,pl_array[p].B/dist);
        else
            lightcolor = Vector3f(pl_array[p].R/dist_squared,pl_array[p].G/dist_squared,pl_array[p].B/dist_squared);
        
        distance.normalize();
        
        float d = distance.dot(surfaceNormal);
        float color = fmax(d, 0);
        color *= d;
        
        Vector3f new_diffuse(material_array[s.mat_num].kdr*lightcolor(0),material_array[s.mat_num].kdg*lightcolor(1),material_array[s.mat_num].kdb*lightcolor(2));
        
        R += color * new_diffuse;
        
        Vector3f rd = (-1 * distance) + (2 * distance.dot(surfaceNormal) * surfaceNormal);
        rd.normalize();
        Vector3f v_vec(-x, -y, -z);
        v_vec.normalize();
        float rvec = rd.dot(v_vec);
        rvec = fmax(rvec, 0);
        Vector3f new_specular(material_array[s.mat_num].ksr*lightcolor(0),material_array[s.mat_num].ksg*lightcolor(1),material_array[s.mat_num].ksb*lightcolor(2));
        
        
        R += pow(rvec, material_array[s.mat_num].ksp) * new_specular;
        
    }
    
    for (int p = 0; p < directional_light_count; p++) {
        Vector3f lightdirection(dl_array[p].x, dl_array[p].y, dl_array[p].z);
        Vector3f lightcolor(dl_array[p].R,dl_array[p].G,dl_array[p].B);
        lightdirection.normalize();
        
        
        float d = lightdirection.dot(surfaceNormal);
        float color = fmax(d, 0);
        color *= d;
        
        Vector3f new_diffuse(material_array[s.mat_num].kdr*lightcolor(0),material_array[s.mat_num].kdg*lightcolor(1),material_array[s.mat_num].kdb*lightcolor(2));
        
        R += color * new_diffuse;
        
        Vector3f rd = (-1 * lightdirection) + (2 * lightdirection.dot(surfaceNormal) * surfaceNormal);
        rd.normalize();
        Vector3f v_vec(-x, -y, -z);
        v_vec.normalize();
        float rvec = rd.dot(v_vec);
        Vector3f new_specular(material_array[s.mat_num].ksr*lightcolor(0),material_array[s.mat_num].ksg*lightcolor(1),material_array[s.mat_num].ksb*lightcolor(2));
        
        R += pow(rvec, material_array[s.mat_num].ksp) * new_specular;
        
    }
    
    
    if (R(0) > 255) {
        R(0) = 255;
    }
    if (R(1) > 255) {
        R(1) = 255;
    }
    if (R(2) > 255) {
        R(2) = 255;
    }
    
    struct pixel R1 = {R(0), R(1), R(2)};
    
    return R1;
}

pixel calculateTrianglePhongShading(float x, float y, float z, Vector3f R, triangle t, int depth) {
    
    Vector3f piece(x, y, z);
    Vector3f surfaceNormal = t.getNormal();
    Vector3f view_ray_direction = piece;
    view_ray_direction.normalize();
    surfaceNormal.normalize();
    piece_tri = t;
    
    for (int a = 0; a < ambient_light_count; a++){
        R += Vector3f(al_array[a].R * material_array[t.mat_num].kar, al_array[a].G * material_array[t.mat_num].kag, al_array[a].B * material_array[t.mat_num].kab);
    }
    
    if ((material_array[t.mat_num].krr > 0 || material_array[t.mat_num].krg > 0 || material_array[t.mat_num].krb > 0) && depth > 0) {
        Vector3f reflected_direction = 2 * (surfaceNormal.dot(-view_ray_direction) * surfaceNormal + view_ray_direction) - view_ray_direction;
        reflected_direction.normalize();
        ray reflected_ray;
        reflected_ray.origin = piece;
        reflected_ray.distance = reflected_direction;
        Vector3f resolved_reflected_color = reflected_color(reflected_ray, depth-1, piece, true);
        if (resolved_reflected_color(0) > 0 || resolved_reflected_color(1) > 0 || resolved_reflected_color(2) > 0){
            Vector3f applied_reflection_color(resolved_reflected_color(0) * material_array[t.mat_num].krr, resolved_reflected_color(1) * material_array[t.mat_num].krg, resolved_reflected_color(2) * material_array[t.mat_num].krb);
            R(0) += applied_reflection_color(0);
            R(1) += applied_reflection_color(1);
            R(2) += applied_reflection_color(2);
        }
    }
    
    
    for (int p = 0; p < point_light_count; p++) {
        
        Vector3f lightloc(pl_array[p].x, pl_array[p].y, pl_array[p].z);
        Vector3f lightcolor;
        ray lightray;
        lightray.origin = lightloc;
        Vector3f distance = lightloc - piece;
        lightray.distance = piece - lightloc;
        float dist = sqrt(distance(0)*distance(0) + distance(1)*distance(1) + distance(2)*distance(2));
        float dist_squared = dist * dist;
        if(!pl_array[p].falloff)
            lightcolor = Vector3f(pl_array[p].R,pl_array[p].G,pl_array[p].B);
        else if(pl_array[p].falloff == 1)
            lightcolor = Vector3f(pl_array[p].R/dist,pl_array[p].G/dist,pl_array[p].B/dist);
        else
            lightcolor = Vector3f(pl_array[p].R/dist_squared,pl_array[p].G/dist_squared,pl_array[p].B/dist_squared);
        distance.normalize();
        
        if(occluded_pl(x, y, z, lightray))
            continue;
        
        float d = distance.dot(surfaceNormal);
        float color = fmax(d, 0);
        color *= d;
        
        Vector3f new_diffuse(material_array[t.mat_num].kdr*lightcolor(0),material_array[t.mat_num].kdg*lightcolor(1),material_array[t.mat_num].kdb*lightcolor(2));
        
        R += color * new_diffuse;
        
        Vector3f rd = (-1 * distance) + (2 * distance.dot(surfaceNormal) * surfaceNormal);
        rd.normalize();
        Vector3f v_vec(-x, -y, -z);
        v_vec.normalize();
        float rvec = rd.dot(v_vec);
        rvec = fmax(rvec, 0);
        Vector3f new_specular(material_array[t.mat_num].ksr*lightcolor(0),material_array[t.mat_num].ksg*lightcolor(1),material_array[t.mat_num].ksb*lightcolor(2));
        
        R += pow(rvec, material_array[t.mat_num].ksp) * new_specular;
        
    }
    
    for (int p = 0; p < directional_light_count; p++) {
        Vector3f lightdirection(dl_array[p].x, dl_array[p].y, dl_array[p].z);
        Vector3f lightcolor(dl_array[p].R,dl_array[p].G,dl_array[p].B);
        lightdirection.normalize();
        
        float d = lightdirection.dot(surfaceNormal);
        float color = fmax(d, 0);
        color *= d;
        
        Vector3f new_diffuse(material_array[t.mat_num].kdr*lightcolor(0),material_array[t.mat_num].kdg*lightcolor(1),material_array[t.mat_num].kdb*lightcolor(2));
        
        R += color * new_diffuse;
        
        Vector3f rd = (-1 * lightdirection) + (2 * lightdirection.dot(surfaceNormal) * surfaceNormal);
        rd.normalize();
        Vector3f v_vec(-x, -y, -z);
        v_vec.normalize();
        float rvec = rd.dot(v_vec);
        Vector3f new_specular(material_array[t.mat_num].ksr*lightcolor(0),material_array[t.mat_num].ksg*lightcolor(1),material_array[t.mat_num].ksb*lightcolor(2));
        
        R += pow(rvec, material_array[t.mat_num].ksp) * new_specular;
        
    }
    
    
    if (R(0) > 255) {
        R(0) = 255;
    }
    if (R(1) > 255) {
        R(1) = 255;
    }
    if (R(2) > 255) {
        R(2) = 255;
    }
    
    struct pixel R1 = {R(0), R(1), R(2)};
    
    return R1;
}

//****************************************************
// Helpers for lookAtObjects method
//****************************************************
float computeSphereIntersection(sphere s, float x, float y) {
    float a = ((x*x)*1.0 + (y*y)*1.0 + 1.0);
    float b = -2.0*x*s.cx - 2.0*y*s.cy + 2.0*s.cz;
    float c = s.cx*s.cx + s.cy*s.cy + s.cz*s.cz - s.r*s.r;
    
    float discriminant = b*b - 4*a*c;
    if (discriminant <= 0) {
        return FLT_MAX;
    }
    else {
        return min((-b + sqrt(discriminant))/(2*a), (-b - sqrt(discriminant))/(2*a));
    }
}

float computeTriangleIntersection(triangle t, float x, float y) {
    Matrix3f A;
    Vector3f b;
    
    Vector3f e1(t.bx - t.ax, t.by - t.ay, t.bz - t.az);
    Vector3f e2(t.cx - t.ax, t.cy - t.ay, t.cz - t.az);
    
    Vector3f cross_vector = e2.cross(e1);
    
    if (abs(cross_vector(0)) < EPSILON && abs(cross_vector(1)) < EPSILON && abs(cross_vector(2)) < EPSILON) {
        return FLT_MAX;
    }
    
    b << -t.ax, -t.ay, -t.az;
    A << -x, e1(0), e2(0),
    -y, e1(1), e2(1),
    1, e1(2), e2(2);
    
    Vector3f x_vector = A.colPivHouseholderQr().solve(b);
    
    if (x_vector(0) >= 0 && x_vector(1) >= 0 && x_vector(2) >= 0) {
        if (x_vector(1) + x_vector(2) <= 1)
            return x_vector(0);
    }
    return FLT_MAX;
}





//****************************************************
// Apply transformation to sphere
//****************************************************
void transformObjects() {
    if (transformation_count > 0) {
        for (int j = 0; j < sphere_count; j++) {
            Vector4f tempSphere(sphere_array[j].cx, sphere_array[j].cy, sphere_array[j].cz, 1);
            for (int i = transformation_count - 1; i >= 0; i--) {
                transformation curTransformation = transformation_array[i];
                if (curTransformation.type == 't') {
                    Matrix4f tempTranslator;
                    tempTranslator << 1, 0, 0, curTransformation.x,
                    0, 1, 0, curTransformation.y,
                    0, 0, 1, curTransformation.z,
                    0, 0, 0, 1;
                    Vector4f newSphere = tempTranslator*tempSphere;
                    sphere_array[j].cx = newSphere(0);
                    sphere_array[j].cy = newSphere(1);
                    sphere_array[j].cz = newSphere(2);
                }
                else if (curTransformation.type == 'r') {
                    // Do Rotation
                }
                else if (curTransformation.type == 's') {
                    Matrix4f tempScaler;
                    tempScaler << curTransformation.x, 0, 0, 0,
                    0, curTransformation.y, 0, 0,
                    0, 0, curTransformation.z, 0,
                    0, 0, 0, 1;
                    Vector4f newSphere = tempScaler*tempSphere;
                    sphere_array[j].cx = newSphere(0);
                    sphere_array[j].cy = newSphere(1);
                    sphere_array[j].cz = newSphere(2);
                }
                else {
                    printf("Invalid transformation type: %c. Skipping this one...\n",curTransformation.type);
                }
            }
        }
        
        for (int j = 0; j < triangle_count; j++) {
            Vector4f tempTriA(triangle_array[j].ax, triangle_array[j].ay, triangle_array[j].az, 1);
            Vector4f tempTriB(triangle_array[j].bx, triangle_array[j].by, triangle_array[j].bz, 1);
            Vector4f tempTriC(triangle_array[j].cx, triangle_array[j].cy, triangle_array[j].cz, 1);
            for (int i = transformation_count - 1; i >= 0; i--) {
                transformation curTransformation = transformation_array[i];
                if (curTransformation.type == 't') {
                    Matrix4f tempTranslator;
                    tempTranslator << 1, 0, 0, curTransformation.x,
                    0, 1, 0, curTransformation.y,
                    0, 0, 1, curTransformation.z,
                    0, 0, 0, 1;
                    Vector4f newTriA = tempTranslator*tempTriA;
                    Vector4f newTriB = tempTranslator*tempTriB;
                    Vector4f newTriC = tempTranslator*tempTriC;
                    triangle_array[j].ax = newTriA(0);
                    triangle_array[j].ay = newTriA(1);
                    triangle_array[j].az = newTriA(2);
                    triangle_array[j].bx = newTriB(0);
                    triangle_array[j].by = newTriB(1);
                    triangle_array[j].bz = newTriB(2);
                    triangle_array[j].cx = newTriC(0);
                    triangle_array[j].cy = newTriC(1);
                    triangle_array[j].cz = newTriC(2);
                }
                else if (curTransformation.type == 'r') {
                    Matrix4f tempRotatorX;
                    printf("COS(curTransformation.x): %f\n",cos(curTransformation.x));
                    tempRotatorX << 1, 0, 0, 0,
                    0, cos(curTransformation.x), -sin(curTransformation.x), 0,
                    0, sin(curTransformation.x), cos(curTransformation.x), 0,
                    0, 0, 0, 1;
                    Matrix4f tempRotatorY;
                    tempRotatorY << cos(curTransformation.y), 0, sin(curTransformation.y), 0,
                    0, 1, 0, 0,
                    -sin(curTransformation.y), 0, cos(curTransformation.y), 0,
                    0, 0, 0, 1;
                    Matrix4f tempRotatorZ;
                    tempRotatorZ << cos(curTransformation.z), -sin(curTransformation.z), 0, 0,
                    sin(curTransformation.z), cos(curTransformation.z), 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1;
                    if (curTransformation.x != 0) {
                        Vector4f newTriA = tempRotatorX*tempTriA;
                        Vector4f newTriB = tempRotatorX*tempTriB;
                        Vector4f newTriC = tempRotatorX*tempTriC;
                        triangle_array[j].ax = newTriA(0);
                        triangle_array[j].ay = newTriA(1);
                        triangle_array[j].az = newTriA(2);
                        triangle_array[j].bx = newTriB(0);
                        triangle_array[j].by = newTriB(1);
                        triangle_array[j].bz = newTriB(2);
                        triangle_array[j].cx = newTriC(0);
                        triangle_array[j].cy = newTriC(1);
                        triangle_array[j].cz = newTriC(2);
                    }
                    if (curTransformation.y != 0) {
                        Vector4f newTriA = tempRotatorY*tempTriA;
                        Vector4f newTriB = tempRotatorY*tempTriB;
                        Vector4f newTriC = tempRotatorY*tempTriC;
                        triangle_array[j].ax = newTriA(0);
                        triangle_array[j].ay = newTriA(1);
                        triangle_array[j].az = newTriA(2);
                        triangle_array[j].bx = newTriB(0);
                        triangle_array[j].by = newTriB(1);
                        triangle_array[j].bz = newTriB(2);
                        triangle_array[j].cx = newTriC(0);
                        triangle_array[j].cy = newTriC(1);
                        triangle_array[j].cz = newTriC(2);
                    }
                    if (curTransformation.z != 0) {
                        Vector4f newTriA = tempRotatorZ*tempTriA;
                        Vector4f newTriB = tempRotatorZ*tempTriB;
                        Vector4f newTriC = tempRotatorZ*tempTriC;
                        triangle_array[j].ax = newTriA(0);
                        triangle_array[j].ay = newTriA(1);
                        triangle_array[j].az = newTriA(2);
                        triangle_array[j].bx = newTriB(0);
                        triangle_array[j].by = newTriB(1);
                        triangle_array[j].bz = newTriB(2);
                        triangle_array[j].cx = newTriC(0);
                        triangle_array[j].cy = newTriC(1);
                        triangle_array[j].cz = newTriC(2);
                    }
                }
                else if (curTransformation.type == 's') {
                    Matrix4f tempScaler;
                    tempScaler << curTransformation.x, 0, 0, 0,
                    0, curTransformation.y, 0, 0,
                    0, 0, curTransformation.z, 0,
                    0, 0, 0, 1;
                    Vector4f newTriA = tempScaler*tempTriA;
                    Vector4f newTriB = tempScaler*tempTriB;
                    Vector4f newTriC = tempScaler*tempTriC;
                    triangle_array[j].ax = newTriA(0);
                    triangle_array[j].ay = newTriA(1);
                    triangle_array[j].az = newTriA(2);
                    triangle_array[j].bx = newTriB(0);
                    triangle_array[j].by = newTriB(1);
                    triangle_array[j].bz = newTriB(2);
                    triangle_array[j].cx = newTriC(0);
                    triangle_array[j].cy = newTriC(1);
                    triangle_array[j].cz = newTriC(2);
                }
                else {
                    printf("Invalid transformation type: %c. Skipping this one...\n",curTransformation.type);
                }
            }
        }
    }
}


//****************************************************
// Iterate thru all objects in scene and return object w/ closest intersection
//****************************************************
pixel lookAtObjects(float x, float y) {
    float t_intersect = FLT_MAX;
    bool isTriangle = false;
    triangle nearestTriangle;
    sphere nearestSphere;
    for (int i = 0; i < sphere_count; i++) {
        float temp = computeSphereIntersection(sphere_array[i], x, y);
        if (temp < t_intersect) {
            t_intersect = temp;
            nearestSphere = sphere_array[i];
        }
    }
    for (int p = 0; p < triangle_count; p++) {
        float temp = computeTriangleIntersection(triangle_array[p], x, y);
        if (temp < t_intersect) {
            t_intersect = temp;
            nearestTriangle = triangle_array[p];
            isTriangle = true;
        }
    }
    
    pixel temp = {0, 0, 0};
    // No intersection, paint black pixel
    if (t_intersect == FLT_MAX) {
        return temp;
    }
    // First intersection is a triangle, calculate for triangle
    else if (isTriangle == true) {
        if(t_intersect > EPSILON){
            Vector3f R(0.0, 0.0, 0.0);
            temp = calculateTrianglePhongShading(t_intersect*x, t_intersect*y, -t_intersect, R, nearestTriangle, RECURSION_DEPTH);
            return temp;
        }
    }
    // First intersection is a sphere, calculate for sphere
    else {
        if(t_intersect > EPSILON){
            Vector3f R(0.0, 0.0, 0.0);
            temp = calculatePhongShading(t_intersect*x, t_intersect*y, -t_intersect, R, nearestSphere, RECURSION_DEPTH);
            return temp;
        }
    }
    return temp;
}



//****************************************************
// Trace a ray through pixel at (x, y)
//****************************************************
pixel traceRay(int x, int y) {
    float x1 = -1.0 + ((2.0*x)/1000.0);
    float y1 = -1.0 + ((2.0*y)/1000.0);
    pixel temp = lookAtObjects(x1, y1);
    return(temp);
}


//****************************************************
// Fill the pixel buffer with traced rays
//****************************************************
void fillBuffer() {
    for (int i = 0; i < WIDTH; i++) {
        for (int j = 0; j < HEIGHT; j++){
            pixel temp = traceRay(i, j);
            image_buffer[i][HEIGHT-j-1].R = temp.R;
            image_buffer[i][HEIGHT-j-1].G = temp.G;
            image_buffer[i][HEIGHT-j-1].B = temp.B;
        }
    }
}



// Returns a random float between 0 and 1
float random_float()
{
    float the_number = (float)rand()/RAND_MAX;
    printf("rand: %f\n",the_number);
    return the_number;
}


struct point {
    float x,y,z;
};

vector<triangle> objParser(const char* szFileName) {
    ifstream myfile; //use all the method of the ifstream
    vector<string*> coord;
    vector<point> vertices;
    vector<triangle> triangles;
    
    myfile.open(szFileName);
    if(myfile.is_open())
    {
        //The file is open,use while loop to check file
        char buf[256]; //create buffer to read line
        while(!myfile.eof())//while we are not at the end of the file.
        {
            myfile.getline(buf,256);
            coord.push_back(new string(buf));
        }
        for(int i = 0; i <coord.size(); i++)
        {
            //check if it's a comment or not
            if(coord[i]->c_str()[0]=='#')
                continue;
            else if(coord[i]->c_str()[0]=='v' && coord[i]->c_str()[1]==' ')
            {
                char tmp;
                float tmp_x,tmp_y,tmp_z;
                
                sscanf(coord[i]->c_str(),"v %f %f %f",&tmp_x,&tmp_y,&tmp_z); //read in the 3
                point p = {tmp_x, tmp_y, tmp_z};
                vertices.push_back(p);
            }
            else if (coord[i]->c_str()[0]=='f' && coord[i]->c_str()[1]==' ') {
                char tmp;
                int tmp_x,tmp_y,tmp_z;
                sscanf(coord[i]->c_str(),"f %d %d %d",&tmp_x,&tmp_y,&tmp_z);
                triangle thisTriangle;
                thisTriangle.ax = vertices.at(tmp_x-1).x;
                thisTriangle.ay = vertices.at(tmp_x-1).y;
                thisTriangle.az = vertices.at(tmp_x-1).z;
                thisTriangle.bx = vertices.at(tmp_y-1).x;
                thisTriangle.by = vertices.at(tmp_y-1).y;
                thisTriangle.bz = vertices.at(tmp_y-1).z;
                thisTriangle.cx = vertices.at(tmp_z-1).x;
                thisTriangle.cy = vertices.at(tmp_z-1).y;
                thisTriangle.cz = vertices.at(tmp_z-1).z;
                thisTriangle.mat_num = material_count - 1;
                triangles.push_back(thisTriangle);
            }
        }
        //Delete coord to avoid memory leaks
        for(int i = 0; i < coord.size();i++)
            delete coord[i];
    }
    else{
        //Display error message, cause the file connot be opened
        cout << "Error occured while opening file" << endl;
    }
    //once complete close file.
    myfile.close();
    printf("%lu", triangles.size());
    return triangles;
}


//****************************************************
// Parse Input
//****************************************************
void parseInput(int argc, char *argv[]) {
    printf("parsing input");
    if (argc == 1) {
        printf("No arguments given, drawing a blank image...\n");
    }
    int i = 1;
    while (i < argc-1) {
        if (strcmp(argv[i], "cam") == 0) {
            if (camera_count == 1) {
                printf("ERROR: There may only be one camera in a given scene. Skipping this.\n");
            }
            else {
                camera temp;
                temp.ex = atof(argv[i+1]);
                temp.ey = atof(argv[i+2]);
                temp.ez = atof(argv[i+3]);
                temp.llx = atof(argv[i+4]);
                temp.lly = atof(argv[i+5]);
                temp.llz = atof(argv[i+6]);
                temp.lrx = atof(argv[i+7]);
                temp.lry = atof(argv[i+8]);
                temp.lrz = atof(argv[i+9]);
                temp.ulx = atof(argv[i+10]);
                temp.uly = atof(argv[i+11]);
                temp.ulz = atof(argv[i+12]);
                temp.urx = atof(argv[i+13]);
                temp.uly = atof(argv[i+14]);
                temp.urz = atof(argv[i+15]);
                the_camera = temp;
                camera_count += 1;
            }
            i += 16;
        }
        else if (strcmp(argv[i], "sph") == 0) {
            sphere temp;
            temp.cx = atof(argv[i+1]);
            temp.cy = atof(argv[i+2]);
            temp.cz = atof(argv[i+3]);
            temp.r = atof(argv[i+4]);
            temp.mat_num = material_count-1;
            sphere_array[sphere_count] = temp;
            sphere_count += 1;
            i += 5;
        }
        else if (strcmp(argv[i], "tri") == 0) {
            triangle temp;
            temp.ax = atof(argv[i+1]);
            temp.ay = atof(argv[i+2]);
            temp.az = atof(argv[i+3]);
            temp.bx = atof(argv[i+4]);
            temp.by = atof(argv[i+5]);
            temp.bz = atof(argv[i+6]);
            temp.cx = atof(argv[i+7]);
            temp.cy = atof(argv[i+8]);
            temp.cz = atof(argv[i+9]);
            temp.mat_num = material_count-1;
            triangle_array[triangle_count] = temp;
            triangle_count += 1;
            i += 10;
        }
        else if (strcmp(argv[i], "ltp") == 0) {
            point_light temp;
            temp.x = atof(argv[i+1]);
            temp.y = atof(argv[i+2]);
            temp.z = atof(argv[i+3]);
            temp.R = atof(argv[i+4]);
            temp.G = atof(argv[i+5]);
            temp.B = atof(argv[i+6]);
            temp.falloff = atof(argv[i+7]);
            pl_array[point_light_count] = temp;
            point_light_count += 1;
            i += 8;
        }
        else if (strcmp(argv[i], "ltd") == 0) {
            directional_light temp;
            temp.x = atof(argv[i+1]);
            temp.y = atof(argv[i+2]);
            temp.z = atof(argv[i+3]);
            temp.R = atof(argv[i+4]);
            temp.G = atof(argv[i+5]);
            temp.B = atof(argv[i+6]);
            dl_array[directional_light_count] = temp;
            directional_light_count += 1;
            i += 7;
        }
        else if (strcmp(argv[i], "lta") == 0) {
            ambient_light temp;
            temp.R = atof(argv[i+1]);
            temp.G = atof(argv[i+2]);
            temp.B = atof(argv[i+3]);
            al_array[ambient_light_count] = temp;
            printf("added ambient");
            ambient_light_count += 1;
            i += 4;
        }
        else if (strcmp(argv[i], "mat") == 0) {
            material temp;
            temp.kar = atof(argv[i+1]);
            temp.kag = atof(argv[i+2]);
            temp.kab = atof(argv[i+3]);
            temp.kdr = atof(argv[i+4]);
            temp.kdg = atof(argv[i+5]);
            temp.kdb = atof(argv[i+6]);
            temp.ksr = atof(argv[i+7]);
            temp.ksg = atof(argv[i+8]);
            temp.ksb = atof(argv[i+9]);
            temp.ksp = atof(argv[i+10]);
            temp.krr = atof(argv[i+11]);
            temp.krg = atof(argv[i+12]);
            temp.krb = atof(argv[i+13]);
            printf("added material");
            material_array[material_count] = temp;
            material_count += 1;
            i += 14;
        }
        else if (strcmp(argv[i], "obj") == 0) {
            if(triangle_count >= 10)
                printf("ERROR: Reached triangle count maximum. Apparently triangles are hard, so this is for your own good");
            else if(!material_count)
                printf("ERROR: No material to color these imported objects with");
            else{
                for(auto &tri : objParser(argv[i+1])){
                    if(triangle_count >= 10)
                        break;
                    triangle_array[triangle_count] = tri;
                    tri.print();
                    triangle_count++;
                }
            }
            i += 2;
        }
        else if (strcmp(argv[i], "xft") == 0) {
            transformation temp;
            temp.type = 't';
            temp.x = atof(argv[i+1]);
            temp.y = atof(argv[i+2]);
            temp.z = atof(argv[i+3]);
            transformation_array[transformation_count] = temp;
            transformation_count += 1;
            i += 4;
        }
        else if (strcmp(argv[i], "xfr") == 0) {
            transformation temp;
            temp.type = 'r';
            temp.x = atof(argv[i+1]);
            temp.y = atof(argv[i+2]);
            temp.z = atof(argv[i+3]);
            transformation_array[transformation_count] = temp;
            transformation_count += 1;
            i += 4;
        }
        else if (strcmp(argv[i], "xfs") == 0) {
            transformation temp;
            temp.type = 's';
            temp.x = atof(argv[i+1]);
            temp.y = atof(argv[i+2]);
            temp.z = atof(argv[i+3]);
            transformation_array[transformation_count] = temp;
            transformation_count += 1;
            i += 4;
        }
        else if (strcmp(argv[i], "xfz") == 0) {
            transformation_count = 0;
            i += 1;
        }
        else if (strcmp(argv[i], "filename") == 0) {
            int len = strlen(argv[i+1]);
            printf("filename len: %i\n",len);
            image_name_str = argv[i+1];
            strncpy(&image_name_str[0], argv[i+1], len);
            (image_name_str + ".bmp").c_str();
            i += 2;
        }
        else {
            printf("ERROR! Unknown argument '%s'\n",argv[i]);
            i += 1000;
        }   
    }
}

//****************************************************
// The Draw Function, Opens Window and Fills in Pixels from Pixel Buffer
//****************************************************
void drawImage() {
    clock_t t1,t2;
    t1 = clock();
    fillBuffer();
    for (int i = 0; i < WIDTH; i++) {
        for (int j = 0; j < HEIGHT; j++){
            main_image(i, j, 0) = image_buffer[i][j].R;
            main_image(i, j, 1) = image_buffer[i][j].G;
            main_image(i, j, 2) = image_buffer[i][j].B;
        }
    }
    CImgDisplay draw_disp(main_image, image_name_str.c_str());
    t2 = clock();
    printf("Rendering completed in %.2f sec\n", ((float)t2 - (float)t1) / CLOCKS_PER_SEC);
    main_image.save((image_name_str + ".bmp").c_str());
    printf("Saved image to %s\n", (image_name_str + ".bmp").c_str());
    fclose(stdout);
    while (!draw_disp.is_closed()) {
        draw_disp.wait();
    }
    
}

//****************************************************
// The Main Function, Calls other big Functions
//****************************************************
int main(int argc, char *argv[]) {
    srand(time(NULL));
    std::stringstream image_name;
    for(int argi = 1; argi < argc; argi++)
        image_name << argv[argi] << "_";
    image_name << "recur_depth_" << RECURSION_DEPTH;
    image_name_str = image_name.str();
    if(image_name_str.length() > 40)
        image_name_str = "default";
    parseInput(argc, argv);
    transformObjects();
    drawImage();
    return 0;
}
