
#include <vector>
#include <iostream>
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


#define PI 3.14159265  // Should be used from mathlib
inline float sqr(float x) { return x*x; }


//****************************************************
// Global Variables
//****************************************************


// Initialize the CImg file with its dimensions and color
CImg<unsigned char> main_image(1000, 700, 1, 3, 0);

struct camera {
  float ex, ey, ez, llx, lly, llz, lrx, lry, lrz, ulx, uly, ulz, urx, ury, urz;
};
struct sphere {
  float cx, cy, cz, r;
};
struct triangle {
  float ax, ay, az, bx, by, bz, cx, cy, cz;
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

ambient_light al_array[5];
point_light pl_array[5];
directional_light dl_array[5];
triangle triangle_array[5];
sphere sphere_array[5];
material material_array[5];

pixel image_buffer[1000][1000];



//****************************************************
// Helpers for lookAtObjects method
//****************************************************
float computeSphereIntersection(sphere s, float x, float y) {
  float a = ((x*x)*1.0 + (y*y)*1.0 + 1.0);
  float b = -2.0*x*s.cx - 2.0*y*s.cy - 2.0*s.cz;
  float c = s.cx*s.cx + s.cy*s.cy + s.cz*s.cz - s.r*s.r;

  float discriminant = b*b - 4*a*c;
  printf("discriminant: %f\n",discriminant);
  if (discriminant <= 0) {
    return FLT_MAX;
  }
  else {
    return min((-b + discriminant)/2*a, (-b - discriminant)/2*a);
  }
}

float computeTriangleIntersection(triangle t, int x, int y) {

}


//****************************************************
// Iterate thru all objects in scene and return object w/ closest intersection
// 
//****************************************************
pixel lookAtObjects(float x, float y) {
  float minDist = FLT_MAX;
  bool isTriangle = false;
  triangle nearestTriangle;
  sphere nearestSphere;
  for (int i = 0; i < sphere_count; i++) {
    float temp = computeSphereIntersection(sphere_array[i], x, y);
    //printf("Temp: %f\n",temp);
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
  //printf("Mindist: %f\n",minDist);
  if (minDist == FLT_MAX) {
    pixel temp = {0, 0, 0};
    return temp;
  }
  else if (isTriangle == true) {

  }
  pixel temp = {255, 255, 255};
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
  for (int i = 0; i < 1000; i++) {
    for (int j = 0; j < 1000; j++){
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
      ambient_light temp;
      temp.R = atof(argv[i+1]);
      temp.G = atof(argv[i+2]);
      temp.B = atof(argv[i+3]);
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
  for (int i = 0; i < 1000; i++) {
    for (int j = 0; j < 700; j++){
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
  drawImage();
  printf("Window closed, program terminated.\n");
  return 0;
}
