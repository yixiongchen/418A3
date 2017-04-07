/*
  CSC418 - RayTracer code - Winter 2017 - Assignment 3&4

  Written Dec. 9 2010 - Jan 20, 2011 by F. J. Estrada
  Freely distributable for adacemic purposes only.

  Uses Tom F. El-Maraghi's code for computing inverse
  matrices. You will need to compile together with
  svdDynamic.c

  You need to understand the code provided in
  this file, the corresponding header file, and the
  utils.c and utils.h files. Do not worry about
  svdDynamic.c, we need it only to compute
  inverse matrices.

  You only need to modify or add code in sections
  clearly marked "TO DO"
*/

#include "utils.h"
#include <math.h>
#include  <time.h>


// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
int MAX_DEPTH;

void buildScene(void)
{
 // Sets up all objects in the scene. This involves creating each object,
 // defining the transformations needed to shape and position it as
 // desired, specifying the reflectance properties (albedos and colours)
 // and setting up textures where needed.
 // Light sources must be defined, positioned, and their colour defined.
 // All objects must be inserted in the object_list. All light sources
 // must be inserted in the light_list.
 //
 // To create hierarchical objects:
 //   Copy the transform matrix from the parent node to the child, and
 //   apply any required transformations afterwards.
 //
 // NOTE: After setting up the transformations for each object, don't
 //       forget to set up the inverse transform matrix!

 struct object3D *o;
 struct pointLS *l;
 struct point3D p;

 ///////////////////////////////////////
 // TO DO: For Assignment 3 you have to use
 //        the simple scene provided
 //        here, but for Assignment 4 you
 //        *MUST* define your own scene.
 //        Part of your mark will depend
 //        on how nice a scene you
 //        create. Use the simple scene
 //        provided as a sample of how to
 //        define and position objects.
 ///////////////////////////////////////

 // Simple scene for Assignment 3:
 // Insert a couple of objects. A plane and two spheres
 // with some transformations.

 // Let's add a plane
 // Note the parameters: ra, rd, rs, rg, R, G, B, alpha, r_index, and shinyness)
 //o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2);
 o=newPlane(1,0,0,.05,.55,.8,.75,1,1,2);
 //o=newPlane(0,0.75,0.75,0,.55,.8,.75,1,1,2);	// Note the plane is highly-reflective (rs=rg=.75) so we
						// should see some reflections if all is done properly.
						// Colour is close to cyan, and currently the plane is
						// completely opaque (alpha=1). The refraction index is
						// meaningless since alpha=1
 Scale(o,6,6,1);				// Do a few transforms...
 RotateZ(o,PI/1.20);
 RotateX(o,PI/2.25);
 Translate(o,0,-3,10);
 invert(&o->T[0][0],&o->Tinv[0][0]);		// Very important! compute
						// and store the inverse
						// transform for this object!
 insertObject(o,&object_list);			// Insert into object list
 
 // Let's add a couple spheres
 //o=newSphere(.05,.95,.35,.35,1,.25,.25,1,1,6);
 o=newSphere(1,0,0,.35,1,.25,.25,1,1,6);
 Scale(o,.75,.5,1.5);
 RotateY(o,PI/2);
 Translate(o,-1.45,1.1,3.5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 //o=newSphere(.05,.95,.95,.75,.75,.95,.55,1,1,6);
 o=newSphere(1,0,0,.75,.75,.95,.55,1,1,6);
 Scale(o,.5,2.0,1.0);
 RotateZ(o,PI/1.5);
 Translate(o,1.75,1.25,5.0);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);


 //add area light source
 int sx = 1.0;
 int sy = 1.0;
 int lx = 1;
 int ly = 1;
 double px = 0;
 double py =15.5;
 double pz =-5.0;
 double r = 0.95;
 double g = 0.95;
 double b = 0.95;
 double n_p_x = 1;
 double n_p_y = 2;
 double n_p_z = 1;
 addAreaLight(sx, sy, n_p_x, n_p_y, n_p_z, px, py, pz, lx, ly, r, g, b, &object_list, &light_list);

 // End of simple scene for Assignment 3
 // Keep in mind that you can define new types of objects such as cylinders and parametric surfaces,
 // or, you can create code to handle arbitrary triangles and then define objects as surface meshes.
 //
 // Remember: A lot of the quality of your scene will depend on how much care you have put into defining
 //           the relflectance properties of your objects, and the number and type of light sources
 //           in the scene. 
}



//Assignment 4:
void buildOwnScene(void){
  struct object3D *o;
  struct pointLS *l;
  struct point3D p;
 
  //earth
  o=newSphere(.05,.95,.35,0.52,0,0,0,1,1,6);
  Scale(o,2,2,2);
  RotateZ(o,PI/0.5);
  RotateY(o,PI/0.5);
  RotateZ(o,PI/1);
  Translate(o,1.75, 1.25, 5.0);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  loadTexture(o, "./texture/earth.ppm");
  insertObject(o,&object_list);

  
  //moon
  o=newSphere(.05,.45,.85,.75,1,0.7,0.2,1,1,6);
  Scale(o,1,1,1);
  Translate(o,-2.5, 2, 4.5);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  loadTexture(o, "./texture/moon.ppm");
  insertObject(o,&object_list);


  //marz
  o=newSphere(.05,.45,.35,.35,0.99,0.21,0,1,1,6);
  Scale(o,1.1,1.1,1.1);
  Translate(o,-2.5,-0.5,4.5);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  loadTexture(o, "./texture/mars.ppm");
  insertObject(o,&object_list);
  
  
  //space plane
  o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2);       
  Scale(o,10,10,1);       
  RotateZ(o,PI/1.20);
  RotateX(o,PI/2.25);
  Translate(o,0,-3,10);
  invert(&o->T[0][0],&o->Tinv[0][0]);  
  loadTexture(o, "./texture/water.ppm");  
  insertObject(o,&object_list);      


  //cylinder
  o=newCylinder(.05,.45,.45,.01,0.99,0.53,0,1,1,6);
  Scale(o,0.3,2,0.3);
  RotateZ(o,PI/1.6);
  RotateY(o,PI/4);
  Translate(o,0.3,-1.9,4.5);
  invert(&o->T[0][0],&o->Tinv[0][0]);            
  insertObject(o,&object_list);      


 //Insert a single point light source.
 int sx = 1;
 int sy = 1;
 int lx = 5;
 int ly = 5;
 double px = 0;
 double py =15.5;
 double pz =-5.0;
 double r = 0.95;
 double g = 0.95;
 double b = 0.95;
 double n_p_x = 1;
 double n_p_y = 2;
 double n_p_z = 1;
 addAreaLight(sx, sy, n_p_x, n_p_y, n_p_z, px, py, pz, lx, ly, r, g, b, &object_list, &light_list);
}



void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double a, double b, struct colourRGB *col)
{
 // This function implements the shading model as described in lecture. It takes
 // - A pointer to the first object intersected by the ray (to get the colour properties)
 // - The coordinates of the intersection point (in world coordinates)
 // - The normal at the point
 // - The ray (needed to determine the reflection direction to use for the global component, as well as for
 //   the Phong specular component)
 // - The current recursion depth
 // - The (a,b) texture coordinates (meaningless unless texture is enabled)
 //
 // Returns:
 // - The colour for this ray (using the col pointer)
 //

 struct colourRGB tmp_col;	// Accumulator for colour components
 double R,G,B;			// Colour for the object in R G and B

 // This will hold the colour as we process all the components of
 // the Phong illumination model
 tmp_col.R=0;
 tmp_col.G=0;
 tmp_col.B=0;

 if (obj->texImg==NULL)		// Not textured, use object colour
 {

  R=obj->col.R;
  G=obj->col.G;
  B=obj->col.B;
 }
 else
 {
  // Get object colour from the texture given the texture coordinates (a,b), and the texturing function
  // for the object. Note that we will use textures also for Photon Mapping.
  
  obj->textureMap(obj->texImg, a, b, &R,&G,&B);
  
 }

 //////////////////////////////////////////////////////////////
 // TO DO: Implement this function. Refer to the notes for
 // details about the shading model.
 //////////////////////////////////////////////////////////////


 //light source list
 struct pointLS *source;
 // loop over light in the light list
 for(source=light_list; source != NULL; source=source->next){
  double a_tmp = 0;
  double b_tmp = 0;
  int num_lights = 0;
  struct object3D *hit_obj;
  struct point3D p_tmp;
  struct point3D n_tmp;

  //light source point position
  struct point3D p0 = source->p0;
  //light source color
  struct colourRGB lightcolor = source->col; 

  //create a shadow ray from intersection point to light source point
  struct point3D *shadow_d = newPoint(p0.px - p->px, p0.py - p->py, p0.pz - p->pz, 1);
  struct ray3D *shadowRay = newRay(p, shadow_d);

  //find intersection for shadow ray
  double lambda = -1;
  findFirstHit(shadowRay, &lambda, obj, &hit_obj, &p_tmp, &n_tmp, &a_tmp, &b_tmp);
  
  //if there is a intersection in shadow ray path
  if (lambda > 0 && lambda < 1){
    // use ambient for shadow
    tmp_col.R += 0;
    tmp_col.G += 0;
    tmp_col.B += 0;   
  }
  
  //there is no intersection in shadow ray path then compute color
  else{
    //get albedos a, d, s
    struct albedosPhong alb = obj->alb;
    double ka = alb.ra;
    double kd = alb.rd;
    double ks = alb.rs;
    //get light properties
    double la =lightcolor.R;
    double ld =lightcolor.G;
    double ls =lightcolor.B;

    //view point ray
    struct point3D * V_direct = newPoint(-ray->d.px, -ray->d.py, -ray->d.pz, 0.0);
    //normal ray
    struct point3D * N_direct = newPoint(n->px, n->py, n->pz, 0.0);
    //incident light ray
    struct point3D * L_direct = newPoint(p0.px-p->px, p0.py-p->py, p0.pz-p->pz, 0.0);
   
    normalize(L_direct);
    normalize(N_direct);
    normalize(V_direct);

    // calculate R
    double num = dot(L_direct, N_direct);
    double R_px = 2*num * N_direct->px - L_direct->px;
    double R_py = 2*num * N_direct->py - L_direct->py;
    double R_pz = 2*num * N_direct->pz - L_direct->pz;

    //reflected light ray
    struct point3D * R_direct = newPoint(R_px, R_py, R_pz, 0);
    //normalize
    normalize(R_direct); 
    
    //calculate phong
    double max_one = max(0, dot(N_direct, L_direct));

    //phong formula 
    double max_two = max(0, pow(dot(V_direct, R_direct), obj->shinyness));
    double phong = ka + kd*max_one + ks*max_two;

    //update color
    tmp_col.R += R * (phong * la);
    tmp_col.G += G * (phong * ld);
    tmp_col.B += B * (phong * ls);
   
    //free memory
    free(L_direct);
    free(N_direct);
    free(R_direct);
    free(V_direct);
  }
 }

//update color phong model
col->R = tmp_col.R;
col->G = tmp_col.G;
col->B = tmp_col.B;

//printf("in shade color is %f %f %f\n",col->R , col->R , col->R);
//recrusion with depth
 if (depth < MAX_DEPTH){

  struct point3D *new_direct =newPoint(0, 0, 0, 0);
  struct point3D* Ray_d = newPoint(-ray->d.px, -ray->d.py, -ray->d.pz, 0);
  // calculate reflective ray direction
  normalize(n);
  normalize(Ray_d);
  double num = dot(Ray_d, n);
  new_direct->px = 2*num * n->px - Ray_d->px;
  new_direct->py = 2*num * n->py - Ray_d->py;
  new_direct->pz = 2*num * n->pz - Ray_d->pz;
  new_direct->pw = 0;
  normalize(new_direct);
  
  /*
    Implementation: Glossy Reflection
  */
  struct point3D * u = cross(new_direct, n);
  normalize(u);
  struct point3D * v = cross(new_direct, u);
  normalize(v);
  double roughness = 0.05;
  double reflect_R, reflect_G, reflect_B;
  reflect_R = 0;
  reflect_G = 0;
  reflect_B = 0;
  // number of reflection ray
  int count = 0;
  srand(time(NULL));
  // generate direction for each sample ray
  while(count < 10){
    struct point3D* reflect_ray = newPoint(new_direct->px, new_direct->py, new_direct->pz, new_direct->pw);
    double rand_num1 = (double)rand() / (double)RAND_MAX;
    double rand_num2 = (double)rand() / (double)RAND_MAX;
    double theta =  2 * PI * (rand_num1 *roughness);
    double phi = 2 * PI *(rand_num2 *roughness);
    double x = sin(theta) * cos(phi);
    double y = sin(theta) * sin(phi);
    double z = cos(theta);
    // Convert sample to world coordinates using the orthonormal basis
    reflect_ray->px =  x*u->px + y*v->px + z*reflect_ray->px;
    reflect_ray->py =  x*u->py + y*v->py + z*reflect_ray->py;
    reflect_ray->pz =  x*u->pz + y*v->pz + z*reflect_ray->pz;
    reflect_ray->pw =  0;
    normalize(reflect_ray);

    //check relfective angle 
    double dot_v = dot(n, reflect_ray);
    if(dot_v >= 0 ){
      //create reflective ray
      struct ray3D *New_Ray = newRay(p, reflect_ray);
      //call ray trace
      struct colourRGB sub_color;
      sub_color.R = 0;
      sub_color.G = 0;
      sub_color.B = 0;
      // store color from relfective ray
      rayTrace(New_Ray, depth+1, &sub_color, obj);
      reflect_R += sub_color.R;
      reflect_G += sub_color.G;
      reflect_B += sub_color.B;
      count += 1;
      free(New_Ray);
    }
    free(reflect_ray);
  }
  // get average color from samples
  col->R += (reflect_R / count);
  col->G += (reflect_G / count);
  col->B += (reflect_B / count);
  
  //free 
  free(u);
  free(v);
  free(new_direct);
  free(Ray_d);
}
return;

}


void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Find the closest intersection between the ray and any objects in the scene.
 // It returns:
 //   - The lambda at the intersection (or < 0 if no intersection)
 //   - The pointer to the object at the intersection (so we can evaluate the colour in the shading function)
 //   - The location of the intersection point (in p)
 //   - The normal at the intersection point (in n)
 //
 // Os is the 'source' object for the ray we are processing, can be NULL, and is used to ensure we don't 
 // return a self-intersection due to numerical errors for recursive raytrace calls.
 //

 /////////////////////////////////////////////////////////////
 // TO DO: Implement this function. See the notes for
 // reference of what to do in here
 /////////////////////////////////////////////////////////////

 double lambda_closest;
 struct object3D* object;
 struct point3D p_tmp;
 struct point3D n_tmp;
 double a_tmp, b_tmp;
 
 //loop over all objects in scene
 for(object=object_list; object!=NULL; object=object->next){

  //ensure we don't return a self-intersection
  if(object != Os ){
    //find the intersection object and get the lambda
    object->intersect(object, ray, &lambda_closest, &p_tmp, &n_tmp, &a_tmp, &b_tmp);

    if((*lambda == -1 || *lambda > lambda_closest) && ( lambda_closest > 0)){
      // get the lamda
      *lambda = lambda_closest;
      //get the pointer to the object at the intersection
      *obj = object;
      //get intersection point
      *p = p_tmp;
      //get the normal at the intersection
      *n= n_tmp; 
      // get texture coordinate
      *a = a_tmp;
      *b = b_tmp;
    }
  }
 }
}



void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os)
{
 // Ray-Tracing function. It finds the closest intersection between
 // the ray and any scene objects, calls the shading function to
 // determine the colour at this intersection, and returns the
 // colour.
 //
 // Os is needed for recursive calls to ensure that findFirstHit will
 // not simply return a self-intersection due to numerical
 // errors. For the top level call, Os should be NULL. And thereafter
 // it will correspond to the object from which the recursive
 // ray originates.
 //

 double lambda;		// Lambda at intersection
 double a,b;		// Texture coordinates
 struct object3D *obj;	// Pointer to object at intersection
 struct point3D p;	// Intersection point
 struct point3D n;	// Normal at intersection
 struct colourRGB I;	// Colour returned by shading function

 if (depth>MAX_DEPTH)	// Max recursion depth reached. Return invalid colour.
 {
  return;
 }

 ///////////////////////////////////////////////////////
 // TO DO: Complete this function. Refer to the notes
 // if you are unsure what to do here.
 ///////////////////////////////////////////////////////

//find the first hitted object in object_list
lambda = -1;
a = 0;
b = 0;
findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);
//printf("a is %f  b is %f\n", a, b);
// there is not intersection
if(lambda <= 0){
  return;
}

if(obj == NULL){
  return;
}

//there is Intersection point exist
if (lambda > 0 ){
  //shading and phong
 
  rtShade(obj, &p, &n, ray, depth, a, b, &I); 

  //handle global phong: update current color + global
  if(Os != NULL){
    col->R += Os->alb.rg*I.R;
    col->G += Os->alb.rg*I.G;
    col->B += Os->alb.rg*I.B;
  }
  // if ray from camera
  else{
    col->R += I.R;
    col->G += I.G;
    col->B += I.B;
  }
}

}


int main(int argc, char *argv[])
{
 // Main function for the raytracer. Parses input parameters,
 // sets up the initial blank image, and calls the functions
 // that set up the scene and do the raytracing.
 struct image *im;	// Will hold the raytraced image
 struct view *cam;	// Camera and view for this scene
 int sx;		// Size of the raytraced image
 int antialiasing;	// Flag to determine whether antialiaing is enabled or disabled
 char output_name[1024];	// Name of the output file for the raytraced .ppm image
 struct point3D e;		// Camera view parameters 'e', 'g', and 'up'
 struct point3D g;
 struct point3D up;
 double du, dv;			// Increase along u and v directions for pixel coordinates
 struct point3D pc,d;		// Point structures to keep the coordinates of a pixel and
				// the direction or a ray
 struct ray3D *ray;		// Structure to keep the ray from e to a pixel
 struct colourRGB col;		// Return colour for raytraced pixels
 struct colourRGB background;   // Background colour
 int i,j;			// Counters for pixel coordinates
 unsigned char *rgbIm;

 if (argc<5)
 {
  fprintf(stderr,"RayTracer: Can not parse input parameters\n");
  fprintf(stderr,"USAGE: RayTracer size rec_depth antialias output_name\n");
  fprintf(stderr,"   size = Image size (both along x and y)\n");
  fprintf(stderr,"   rec_depth = Recursion depth\n");
  fprintf(stderr,"   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
  fprintf(stderr,"   output_name = Name of the output file, e.g. MyRender.ppm\n");
  exit(0);
 }
 sx=atoi(argv[1]);
 MAX_DEPTH=atoi(argv[2]);
 if (atoi(argv[3])==0) antialiasing=0; else antialiasing=1;
 strcpy(&output_name[0],argv[4]);

 fprintf(stderr,"Rendering image at %d x %d\n",sx,sx);
 fprintf(stderr,"Recursion depth = %d\n",MAX_DEPTH);
 if (!antialiasing) fprintf(stderr,"Antialising is off\n");
 else fprintf(stderr,"Antialising is on\n");
 fprintf(stderr,"Output file name: %s\n",output_name);

 object_list=NULL;
 light_list=NULL;

 // Allocate memory for the new image
 im=newImage(sx, sx);
 if (!im)
 {
  fprintf(stderr,"Unable to allocate memory for raytraced image\n");
  exit(0);
 }
 else rgbIm=(unsigned char *)im->rgbdata;

 ///////////////////////////////////////////////////
 // TO DO: You will need to implement several of the
 //        functions below. For Assignment 3, you can use
 //        the simple scene already provided. But
 //        for Assignment 4 you need to create your own
 //        *interesting* scene.
 ///////////////////////////////////////////////////
 //buildScene();
 buildOwnScene();		// Create a scene. This defines all the
			// objects in the world of the raytracer

 //////////////////////////////////////////
 // TO DO: For Assignment 3 you can use the setup
 //        already provided here. For Assignment 4
 //        you may want to move the camera
 //        and change the view parameters
 //        to suit your scene.
 //////////////////////////////////////////

 // Mind the homogeneous coordinate w of all vectors below. DO NOT
 // forget to set it to 1, or you'll get junk out of the
 // geometric transformations later on.

 // Camera center is at (0,0,-1)
 e.px=0;
 e.py=0;
 e.pz=-1;
 e.pw=1;

 // To define the gaze vector, we choose a point 'pc' in the scene that
 // the camera is looking at, and do the vector subtraction pc-e.
 // Here we set up the camera to be looking at the origin, so g=(0,0,0)-(0,0,-1)
 g.px=0;
 g.py=0;
 g.pz=1;
 g.pw=0;

 // Define the 'up' vector to be the Y axis
 up.px=0;
 up.py=1;
 up.pz=0;
 up.pw=0;

 // Set up view with given the above vectors, a 4x4 window,
 // and a focal length of -1 (why? where is the image plane?)
 // Note that the top-left corner of the window is at (-2, 2)
 // in camera coordinates.
 cam=setupView(&e, &g, &up, -3, -2, 2, 4);
 if (cam==NULL)
 {
  fprintf(stderr,"Unable to set up the view and camera parameters. Our of memory!\n");
  cleanup(object_list,light_list);
  deleteImage(im);
  exit(0);
 }

 // Set up background colour here
 background.R=0;
 background.G=0;
 background.B=0;

 // Do the raytracing
 //////////////////////////////////////////////////////
 // TO DO: You will need code here to do the raytracing
 //        for each pixel in the image. Refer to the
 //        lecture notes, in particular, to the
 //        raytracing pseudocode, for details on what
 //        to do here. Make sure you undersand the
 //        overall procedure of raytracing for a single
 //        pixel.
 //////////////////////////////////////////////////////
 du=cam->wsize/(sx-1);		// du and dv. In the notes in terms of wl and wr, wt and wb,
 dv=-cam->wsize/(sx-1);		// here we use wl, wt, and wsize. du=dv since the image is
				// and dv is negative since y increases downward in pixel
				// coordinates and upward in camera coordinates.

 fprintf(stderr,"View parameters:\n");
 fprintf(stderr,"Left=%f, Top=%f, Width=%f, f=%f\n",cam->wl,cam->wt,cam->wsize,cam->f);
 fprintf(stderr,"Camera to world conversion matrix (make sure it makes sense!):\n");
 printmatrix(cam->C2W);
 fprintf(stderr,"World to camera conversion matrix\n");
 printmatrix(cam->W2C);
 fprintf(stderr,"\n");

 fprintf(stderr,"Rendering row: ");

 srand((unsigned int)time(NULL));



// do ray trace for object
 #pragma omp parallel for
   for (j=0;j<sx;j++)		// For each of the pixels in the image
   {
    fprintf(stderr,"%d/%d, ",j,sx);
    for (i=0;i<sx;i++)
    {
      ///////////////////////////////////////////////////////////////////
      // TO DO - complete the code that should be in this loop to do the
      //         raytracing!
      ///////////////////////////////////////////////////////////////////

      //implementation of anti-aliasing
      if (antialiasing){
        //split each piexl into sub-pixel (4x4)
        //[i, i+1] and [j, j+1]
        struct point3D sub_pc,sub_d; 
        float ran_num_x, ran_num_y;
        struct colourRGB sub_color;
        sub_color.R = 0.0;
        sub_color.G = 0.0;
        sub_color.B = 0.0;
        int x;
        int y;
        int super_x = 3;
        int super_y = 3;
        // for each sub pixel
        for(x=0; x<super_x; x++){
          for (y=0; y<super_y; y++){
            float random_1 = (float)rand() / (float)RAND_MAX;
            float r_1 = random_1 * 0.25;
            ran_num_x= (x*0.25) + r_1;

            float random_2 = (float)rand() / (float)RAND_MAX;
            float r_2 = random_2 * 0.25;
            ran_num_y= (y*0.25) + r_2;

            sub_pc.px = cam->wl+(i+ran_num_x)*du;
            sub_pc.py = cam->wt+(j+ran_num_y)*dv;
            sub_pc.pz = cam->f;
            sub_pc.pw = 1;

            sub_d.px = sub_pc.px - cam->e.px;
            sub_d.py = sub_pc.py - cam->e.py;
            sub_d.pz = sub_pc.pz - cam->e.pz;
            sub_d.pw = 0;

            //Convert to world-space
            matVecMult(cam->C2W, &sub_pc);
            matVecMult(cam->C2W, &sub_d);
            col.R = 0.0;
            col.G = 0.0;
            col.B = 0.0;
            // create ray
            ray = newRay(&sub_pc, &sub_d);
            // ray trace
            rayTrace(ray, 0, &col, NULL);
            sub_color.R += col.R;
            sub_color.G += col.G;
            sub_color.B += col.B;
          }  
        }

        //Final Color of pixel is the average of all the samples.
        sub_color.R  =  sub_color.R /( super_x * super_y);
        sub_color.G  =  sub_color.G /( super_x * super_y);
        sub_color.B  =  sub_color.B /( super_x * super_y);


        if(sub_color.R > 1){
           rgbIm[(j*sx+i)*3] =255;
        }
        else{
           rgbIm[(j*sx+i)*3] =sub_color.R*255;
        }

        if(sub_color.G > 1){
           rgbIm[(j*sx+i)*3+1]=255;
        }
        else{
           rgbIm[(j*sx+i)*3+1]=sub_color.G*255;
        }

        if(sub_color.B > 1){
           rgbIm[(j*sx+i)*3+2] =255;
        }
        else{
           rgbIm[(j*sx+i)*3+2]=sub_color.B*255;
        }
      }
      
      //no anti-aliasing
      else{

        // Pixel in Camera Coordinates
        pc.px = cam->wl+(i+0.5)*du;
        pc.py = cam->wt+(j+0.5)*dv;
        pc.pz = cam->f;
        pc.pw = 1;

        d.px = pc.px - cam->e.px;
        d.py = pc.py - cam->e.py;
        d.pz = pc.pz - cam->e.pz;
        d.pw = 0;

        //Convert to world-space
        matVecMult(cam->C2W, &pc);
        matVecMult(cam->C2W, &d);
        col.R = 0.0;
        col.G = 0.0;
        col.B = 0.0;
         
        // create ray
        ray = newRay(&pc, &d);
        // ray trace
        rayTrace(ray, 0, &col, NULL);

        // create image
        if(col.R > 1){
           rgbIm[(j*sx+i)*3] =255;
        }
        else{
           rgbIm[(j*sx+i)*3] =col.R*255;
        }

        if(col.G > 1){
           rgbIm[(j*sx+i)*3+1]=255;
        }
        else{
           rgbIm[(j*sx+i)*3+1] =col.G*255;
        }

        if(col.B > 1){
           rgbIm[(j*sx+i)*3+2] =255;
        }
        else{
           rgbIm[(j*sx+i)*3+2] =col.B*255;
        }

     }
    } // end for i
   } // end for j

 fprintf(stderr,"\nDone!\n");

 // Output rendered image
 imageOutput(im,output_name);

 // Exit section. Clean up and return.
 cleanup(object_list,light_list);		// Object and light lists
 deleteImage(im);				// Rendered image
 free(cam);					// camera view
 exit(0);
}
