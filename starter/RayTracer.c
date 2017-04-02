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
 o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2);	// Note the plane is highly-reflective (rs=rg=.75) so we
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
 o=newSphere(.05,.95,.35,.35,1,.25,.25,1,1,6);
 Scale(o,.75,.5,1.5);
 RotateY(o,PI/2);
 Translate(o,-1.45,1.1,3.5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(.05,.95,.95,.75,.75,.95,.55,1,1,6);
 Scale(o,.5,2.0,1.0);
 RotateZ(o,PI/1.5);
 Translate(o,1.75,1.25,5.0);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 // Insert a single point light source.
 p.px=0;
 p.py=15.5;
 p.pz=-5.5;
 p.pw=1;
 l=newPLS(&p,.95,.95,.95);
 insertPLS(l,&light_list);

 // End of simple scene for Assignment 3
 // Keep in mind that you can define new types of objects such as cylinders and parametric surfaces,
 // or, you can create code to handle arbitrary triangles and then define objects as surface meshes.
 //
 // Remember: A lot of the quality of your scene will depend on how much care you have put into defining
 //           the relflectance properties of your objects, and the number and type of light sources
 //           in the scene.
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
  obj->textureMap(obj->texImg,a,b,&R,&G,&B);
 }

 //////////////////////////////////////////////////////////////
 // TO DO: Implement this function. Refer to the notes for
 // details about the shading model.
 //////////////////////////////////////////////////////////////
 double lambda = -1;
 double a_tmp;
 double b_tmp;
 struct object3D *hit_obj;
 struct point3D p_tmp;
 struct point3D n_tmp;
 //light source list
 struct pointLS *source;
 // loop over light in the scene
 for(source=light_list; source != NULL; source=source->next){
  //light source point
  struct point3D p0 = source->p0;
  struct colourRGB lightcolor = source->col; 

  //create a shadow ray from intersection to light
  struct point3D *shadow_d = newPoint(p0.px - p->px, p0.py - p->py, p0.pz - p->pz, 1);
  struct ray3D *shadowRay = newRay(p, shadow_d);

  //variables to indicate object intersected by shadow ray
   
  //find intersection for shadow
  findFirstHit(shadowRay, &lambda, obj, &hit_obj, &p_tmp, &n_tmp, &a_tmp, &b_tmp);
  
  
  //if there is a intersection in shadow ray path
  if (lambda > 0 && lambda < 1){
    col->R = 0;
    col->G = 0;
    col->B = 0;
    
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
    //normalize(R_direct); 
    //calculate phong
    double max_one = max(0, dot(N_direct, L_direct));
    //phong formula 
    double max_two = max(0, pow(dot(V_direct, R_direct), obj->shinyness));
    double phong = ka + kd*max_one + ks*max_two;

    //update color
    tmp_col.R += R * (phong * la);
    tmp_col.G += G * (phong * ld);
    tmp_col.B += B * (phong * ls);

    //printf("col->R is %f\n", col->R);
    free(L_direct);
    free(N_direct);
    free(R_direct);
    free(V_direct);
  }
 }

 // recrusion with depth
//  if (depth < MAX_DEPTH){
//   //reflectray  = reflectRay = create a reflected ray at the intersection;
//   struct ray3D *New_Ray = (struct ray3D*)malloc(sizeof(struct ray3D));
//   struct point3D *new_direct = (struct point3D *)malloc(sizeof(struct point3D));
//   struct point3D* Ray_d = newPoint(-ray->d.px, -ray->d.py, -ray->d.pz, 0);

//   // calculate new reflective ray direction
//   double num = dot(Ray_d, n);
//   new_direct->px = 2*num * n->px - Ray_d->px;
//   new_direct->py = 2*num * n->py - Ray_d->py;
//   new_direct->pz = 2*num * n->pz - Ray_d->pz;
//   new_direct->pw = 0;
//   normalize(new_direct);
//   New_Ray = newRay(p, new_direct);
//   //call ray trace
//   struct colourRGB recur_color;
//   rayTrace(New_Ray, depth+1, &recur_color, obj);
//   tmp_col.R += recur_color.R;
//   tmp_col.G += recur_color.G;
//   tmp_col.B += recur_color.B;
//   free(new_direct);
//   free(New_Ray);
//   free(Ray_d);
// }

//Be sure to update 'col' with the final colour computed here!
col->R = tmp_col.R;
col->G = tmp_col.G;
col->B = tmp_col.B;
printf("shade color is %f, %f, %f\n", col->R, col->G, col->B);


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
 //loop over all objects in scene
 for(object=object_list; object!=NULL; object=object->next){
  
  //ensure we don't return a self-intersection
  if(object == Os){
    continue;
  }
  
  //find the intersection object and get the lambda
  object->intersect(object, ray, &lambda_closest, &p_tmp, &n_tmp, a, b);
  //printf("the lambda computed from ray is %f\n",lambda_tmp);
  if((*lambda == -1 || *lambda > lambda_closest) && ( lambda_closest > 0)){
    // get the lamda
    *lambda = lambda_closest;
    //get the pointer to the object at the intersection
    *obj = object;
    //get intersection point
    *p = p_tmp;
    //geth the normal at the intersection
    *n= n_tmp;  
  }
 }
 return;
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
 //printf("b\n");
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
findFirstHit(ray, &lambda, NULL, &obj, &p, &n, &a, &b);
//printf("lambda in ray is %f\n", lambda);
// there is not intersection
if(lambda <= 0){
  return;
}

//intersection point exist
if (obj != Os && lambda > 0){
  //shading and phong
  //printf("color\n");
  rtShade(obj, &p, &n, ray, depth, a, b, &I); 
  //printf("color is %f, %f, %f\n", I.R, I.G, I.B);
  //assign color to image
  col->R = I.R;
  col->G = I.G;
  col->B = I.B;
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
 buildScene();		// Create a scene. This defines all the
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
 for (j=0;j<sx;j++)		// For each of the pixels in the image
 {
  fprintf(stderr,"%d/%d, ",j,sx);
  for (i=0;i<sx;i++)
  {
    ///////////////////////////////////////////////////////////////////
    // TO DO - complete the code that should be in this loop to do the
    //         raytracing!
    ///////////////////////////////////////////////////////////////////
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
    
    //printf("color is col %f , %f, %f\n", col.R, col.G, col.B);
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
