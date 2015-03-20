
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <cfloat>
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


#define LENGTH 1000
#define WIDTH 1000
#define RECURSION_DEPTH 0
#define EPSILON 0.000001    //Used when determining triangle-ray intersections
#define PI 3.14159265  // Should be used from mathlib
inline float sqr(float x) { return x*x; }


//****************************************************
// Global Variables
//****************************************************


// Initialize the CImg file with its dimensions and color
CImg<unsigned char> main_image(WIDTH, LENGTH, 1, 3, 0);

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
};
struct triangle {
    float ax, ay, az, bx, by, bz, cx, cy, cz;
    int mat_num;
    Vector3f getV1(){ return Vector3f(ax, ay, az);}
    Vector3f getV2(){ return Vector3f(bx, by, bz);}
    Vector3f getV3(){ return Vector3f(cx, cy, cz);}
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

pixel image_buffer[WIDTH][LENGTH];

//*************************************************************//
// Intersection methods for checking intersections with *rays* //
//*************************************************************//

float computeSphereIntersection(sphere s, ray r) {
  float a = 1, dist_x = r.origin(0) - s.cx, dist_y = r.origin(1) - s.cy, dist_z = r.origin(2) - s.cz;
  float b = 2.0 * r.distance(0) * dist_x + 2.0 * r.distance(1) * dist_y + 2.0 * r.distance(2) * dist_z;
  float c = dist_x * dist_x + dist_y * dist_y + dist_z * dist_z - s.r * s.r;
  float discriminant = b*b - 4*c;
  
  if (discriminant <= 0)
    return -1;
  else
    return -min((-b + discriminant)/2*a, (-b - discriminant)/2*a);
}

float computeTriangleIntersection(triangle tri, ray lightray)
{
    Vector3f V1 = tri.getV1(), V2 = tri.getV2(), V3 = tri.getV3(), origin = lightray.origin, distance = lightray.distance;
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
    if(u < 0.f || u > 1.f) return 0;
    
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
// TODO: Add polygon checks
bool occluded_pl(float x, float y, float z, ray lightray){
    int num_intersected_objects = 0;
    for (int s_index = 0; s_index < sphere_count; s_index++)
        if (computeSphereIntersection(sphere_array[s_index], lightray) >= 0)
            num_intersected_objects++;
    for (int t_index = 0; t_index < triangle_count; t_index++)
        if (computeTriangleIntersection(triangle_array[t_index], lightray) >= 0)
            num_intersected_objects++;
    printf('%d', num_intersected_objects);
    
    return num_intersected_objects > 1;
}

pixel calculatePhongShading(float x, float y, float z, Vector3f R, sphere s, int depth) {
    
    Vector3f piece(x, y, z);
    piece.normalize();
    Vector3f surfaceNormal(x-s.cx, y-s.cy, z+s.cz);
    surfaceNormal.normalize();
    
    for (int p = 0; p < point_light_count; p++) {
        Vector3f lightloc(pl_array[p].x, pl_array[p].y, pl_array[p].z);
        Vector3f lightcolor;
        ray lightray;
        lightray.origin = lightloc;
        Vector3f distance = lightloc - Vector3f(x, y, z);
        lightray.distance = distance;
        float dist = sqrt(distance(0)*distance(0) + distance(1)*distance(1) + distance(2)*distance(2));
        float dist_squared = dist * dist;
        if(!pl_array[p].falloff)
            lightcolor = Vector3f(pl_array[p].R,pl_array[p].G,pl_array[p].B);
        else if(pl_array[p].falloff == 1)
            lightcolor = Vector3f(pl_array[p].R/dist,pl_array[p].G/dist,pl_array[p].B/dist);
        else
            lightcolor = Vector3f(pl_array[p].R/dist_squared,pl_array[p].G/dist_squared,pl_array[p].B/dist_squared);
        
        lightloc = lightloc - piece;
        lightloc.normalize();
        
        Vector3f new_ambient(material_array[s.mat_num].kar*lightcolor(0),material_array[s.mat_num].kag*lightcolor(1),material_array[s.mat_num].kab*lightcolor(2));
        R += new_ambient; // Add ambient term for each point light
        
        if(occluded_pl(x, y, z, lightray))
            continue;
        
        float d = lightloc.dot(surfaceNormal);
        float color = fmax(d, 0);
        color *= d;
        
        Vector3f new_diffuse(material_array[s.mat_num].kdr*lightcolor(0),material_array[s.mat_num].kdg*lightcolor(1),material_array[s.mat_num].kdb*lightcolor(2));
        
        R += color * new_diffuse;
        
        Vector3f rd = (-1 * lightloc) + (2 * lightloc.dot(surfaceNormal) * surfaceNormal);
        rd.normalize();
        Vector3f v_vec(-x, -y, z);
        v_vec.normalize();
        float rvec = rd.dot(v_vec);
        rvec = fmax(rvec, 0);
        Vector3f new_specular(material_array[s.mat_num].ksr*lightcolor(0),material_array[s.mat_num].ksg*lightcolor(1),material_array[s.mat_num].ksb*lightcolor(2));
        
        
        R += pow(rvec, material_array[s.mat_num].ksp) * new_specular;
        
        if (material_array[s.mat_num].krr > 0 && depth > 0) {
            Vector3f V(x, y, z);
            Vector3f R1 = -V + 2*surfaceNormal*(V.dot(surfaceNormal));
            Vector3f F_vec = R1+V;
            return calculatePhongShading(F_vec(0), F_vec(1), F_vec(2), R, s, depth-1);
        }
        
    }
    
    for (int p = 0; p < directional_light_count; p++) {
        Vector3f lightloc(dl_array[p].x, dl_array[p].y, dl_array[p].z);
        Vector3f lightcolor(dl_array[p].R,dl_array[p].G,dl_array[p].B);
        lightloc = -1*lightloc;
        lightloc.normalize();
        
        Vector3f new_ambient(material_array[s.mat_num].kar*lightcolor(0),material_array[s.mat_num].kag*lightcolor(1),material_array[s.mat_num].kab*lightcolor(2));
        R += new_ambient; // Add ambient term for each point light
        
        float d = lightloc.dot(surfaceNormal);
        float color = fmax(d, 0);
        color *= d;
        
        Vector3f new_diffuse(material_array[s.mat_num].kdr*lightcolor(0),material_array[s.mat_num].kdg*lightcolor(1),material_array[s.mat_num].kdb*lightcolor(2));
        
        R += color * new_diffuse;
        
        Vector3f rd = (-1 * lightloc) + (2 * lightloc.dot(surfaceNormal) * surfaceNormal);
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
        piece.normalize();
        Vector3f surfaceNormal(-x, -y, z);
        surfaceNormal.normalize();

        for (int p = 0; p < point_light_count; p++) {
          
            Vector3f lightloc(pl_array[p].x, pl_array[p].y, pl_array[p].z);
            Vector3f lightcolor;
            ray lightray;
            lightray.origin = lightloc;
            Vector3f distance = lightloc - Vector3f(x, y, z);
            lightray.distance = distance;
            float dist = sqrt(distance(0)*distance(0) + distance(1)*distance(1) + distance(2)*distance(2));
            float dist_squared = dist * dist;
            if(!pl_array[p].falloff)
                lightcolor = Vector3f(pl_array[p].R,pl_array[p].G,pl_array[p].B);
            else if(pl_array[p].falloff == 1)
                lightcolor = Vector3f(pl_array[p].R/dist,pl_array[p].G/dist,pl_array[p].B/dist);
            else
                lightcolor = Vector3f(pl_array[p].R/dist_squared,pl_array[p].G/dist_squared,pl_array[p].B/dist_squared);
            lightloc = lightloc - piece;
          lightloc.normalize();

          Vector3f new_ambient(material_array[t.mat_num].kar*lightcolor(0),material_array[t.mat_num].kag*lightcolor(1),material_array[t.mat_num].kab*lightcolor(2));
          R += new_ambient; // Add ambient term for each point light

          if(occluded_pl(x, y, z, lightray))
              continue;
          float d = lightloc.dot(surfaceNormal);
          float color = fmax(d, 0);
          color *= d;

          Vector3f new_diffuse(material_array[t.mat_num].kdr*lightcolor(0),material_array[t.mat_num].kdg*lightcolor(1),material_array[t.mat_num].kdb*lightcolor(2));

          R += color * new_diffuse;

          Vector3f rd = (-1 * lightloc) + (2 * lightloc.dot(surfaceNormal) * surfaceNormal);
          rd.normalize();
          Vector3f v_vec(-x, -y, z);
          v_vec.normalize();
          float rvec = rd.dot(v_vec);
          rvec = fmax(rvec, 0);
          Vector3f new_specular(material_array[t.mat_num].ksr*lightcolor(0),material_array[t.mat_num].ksg*lightcolor(1),material_array[t.mat_num].ksb*lightcolor(2));


          R += pow(rvec, material_array[t.mat_num].ksp) * new_specular;

          if (material_array[t.mat_num].krr > 0 && depth > 0) {
            Vector3f V(x, y, z);
            Vector3f R1 = -V + 2*surfaceNormal*(V.dot(surfaceNormal));
            Vector3f F_vec = R1+V;
            return calculateTrianglePhongShading(F_vec(0), F_vec(1), F_vec(2), R, t, depth-1);
          }

        }

        for (int p = 0; p < directional_light_count; p++) {
          Vector3f lightloc(dl_array[p].x, dl_array[p].y, dl_array[p].z);
          Vector3f lightcolor(dl_array[p].R,dl_array[p].G,dl_array[p].B);
          lightloc = -1*lightloc;
          lightloc.normalize();

          Vector3f new_ambient(material_array[t.mat_num].kar*lightcolor(0),material_array[t.mat_num].kag*lightcolor(1),material_array[t.mat_num].kab*lightcolor(2));
          R += new_ambient; // Add ambient term for each point light

          float d = lightloc.dot(surfaceNormal);
          float color = fmax(d, 0);
          color *= d;

          Vector3f new_diffuse(material_array[t.mat_num].kdr*lightcolor(0),material_array[t.mat_num].kdg*lightcolor(1),material_array[t.mat_num].kdb*lightcolor(2));

          R += color * new_diffuse;

          Vector3f rd = (-1 * lightloc) + (2 * lightloc.dot(surfaceNormal) * surfaceNormal);
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
  float b = -2.0*x*s.cx - 2.0*y*s.cy - 2.0*s.cz;
  float c = s.cx*s.cx + s.cy*s.cy + s.cz*s.cz - s.r*s.r;

  float discriminant = b*b - 4*a*c;
  //printf("discriminant: %f\n",discriminant);
  if (discriminant <= 0) {
    return FLT_MAX;
  }
  else {
    return -min((-b + discriminant)/2*a, (-b - discriminant)/2*a);
  }
}

float computeTriangleIntersection(triangle t, float x, float y) {
  Matrix3f A;
  Vector3f b;

  Vector3f v1(t.bx - t.ax, t.by - t.ay, t.bz - t.az);
  Vector3f v2(t.cx - t.ax, t.cy - t.ay, t.cz - t.az);

  Vector3f cross_vector = v2.cross(v1);

  if (cross_vector(0) == 0 && cross_vector(1) == 0 && cross_vector(2) == 0) {
    return FLT_MAX;
  }

  b << -t.ax, -t.ay, -t.az;
  A << x, t.bx - t.ax, t.cx - t.ax,
         y, t.by - t.ay, t.cy - t.ay,
         1, t.bz - t.az, t.cz - t.az;

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
  float minDist = FLT_MAX;
  bool isTriangle = false;
  triangle nearestTriangle;
  sphere nearestSphere;
  for (int i = 0; i < sphere_count; i++) {
    float temp = computeSphereIntersection(sphere_array[i], x, y);
    if (temp < minDist) {
      minDist = temp;
      nearestSphere = sphere_array[i];
    }
  }
  for (int p = 0; p < triangle_count; p++) {
    float temp = computeTriangleIntersection(triangle_array[p], x, y);
    if (temp < minDist) {
      minDist = temp;
      nearestTriangle = triangle_array[p];
      isTriangle = true;
    }
  }

  // No intersection, paint black pixel
  if (minDist == FLT_MAX) {
    pixel temp = {0, 0, 0};
    return temp;
  }
  // First intersection is a triangle, calculate for triangle
  else if (isTriangle == true) {
    Vector3f R(0.0, 0.0, 0.0);
    pixel temp = calculateTrianglePhongShading(-minDist*x, -minDist*y, minDist, R, nearestTriangle, RECURSION_DEPTH);
    //pixel temp = {255, 255, 255};
    return temp;
  }
  // First intersection is a sphere, calculate for sphere
  else {
    Vector3f R(0.0, 0.0, 0.0);
    pixel temp = calculatePhongShading(-minDist*x, -minDist*y, minDist, R, nearestSphere, RECURSION_DEPTH);
    return temp;
  }
}



//****************************************************
// Trace a ray through pixel at (x, y)
//****************************************************
pixel traceRay(int x, int y) {
  float x1 = 1.0 - ((2.0*x)/1000.0);
  float y1 = -1.0 + ((2.0*y)/1000.0);
  pixel temp = lookAtObjects(x1, y1);
  return(temp);
}


//****************************************************
// Fill the pixel buffer with traced rays
//****************************************************
void fillBuffer() {
  for (int i = 0; i < WIDTH; i++) {
    for (int j = 0; j < LENGTH; j++){
        pixel temp = traceRay(i, j);
        image_buffer[i][j].R = temp.R;
        image_buffer[i][j].G = temp.G;
        image_buffer[i][j].B = temp.B;
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


//****************************************************
// Parse Input
//****************************************************
void parseInput(int argc, char *argv[]) {
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
      material_array[material_count] = temp;
      material_count += 1;
      i += 14;
    }
    else if (strcmp(argv[i], "obj") == 0) {
      ambient_light temp;
      temp.R = atof(argv[i+1]);
      temp.G = atof(argv[i+2]);
      temp.B = atof(argv[i+3]);
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
      printf("LEN: %i\n",len);
      imageName = new char[len+4];
      strncpy(imageName, argv[i+1], len);
      imageName[len] = '.';
      imageName[len+1] = 'b';
      imageName[len+2] = 'm';
      imageName[len+3] = 'p';
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
  fillBuffer();
  for (int i = 0; i < WIDTH; i++) {
    for (int j = 0; j < LENGTH; j++){
      main_image(i, j, 0) = image_buffer[i][j].R;
      main_image(i, j, 1) = image_buffer[i][j].G;
      main_image(i, j, 2) = image_buffer[i][j].B;
    }
  }
  CImgDisplay draw_disp(main_image, "Main Image");
  while (!draw_disp.is_closed()) {
    draw_disp.wait();
  }
}

//****************************************************
// The Main Function, Calls other big Functions
//****************************************************
int main(int argc, char *argv[]) {
  srand(time(NULL));
  parseInput(argc, argv);
  transformObjects();
  drawImage();
  main_image.save("Image.bmp");
  printf("Window closed, program terminated.\n");
  return 0;
}
