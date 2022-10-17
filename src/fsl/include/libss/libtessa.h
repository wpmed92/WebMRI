/* {{{ comments */

/*% cc -g sphere.c -o sphere -lm
 *
 * sphere - generate a triangle mesh approximating a sphere by
 *  recursive subdivision. First approximation is an platonic
 *  solid; each level of refinement increases the number of
 *  triangles by a factor of 4.
 *
 * Level 3 (128 triangles for an octahedron) is a good tradeoff if
 *  gouraud shading is used to render the database.
 *
 * Usage: sphere [level] [-p] [-c] [-f] [-t] [-i]
 *      level is an integer >= 1 setting the recursion level (default 1).
 *      -p causes generation of a PPHIGS format ASCII archive
 *          instead of the default generic output format.
 *      -c causes triangles to be generated with vertices in counterclockwise
 *          order as viewed from the outside in a RHS coordinate system.
 *          The default is clockwise order.
 *      -f generates triangle without per-vertex normals (PPHIGS only)
 *      -t starts with a tetrahedron instead of an octahedron
 *      -i starts with a icosahedron instead of an octahedron
 *      -v outputs 3D volume
 *      -x output xtopol format
 *      -r resort data into points list
 *      -n normal output
 *
 *  The subroutines print_object() and print_triangle() should
 *  be changed to generate whatever the desired database format is.
 *
 * Jon Leech (leech@cs.unc.edu) 3/24/89
 * icosahedral code added by Jim Buddenhagen (jb1556@daditz.sbc.com) 5/93
 * various alternative output options including list of points and their links added by Steve Smith (steve@fmrib.ox.ac.uk) 5/99
 *
 */

/* }}} */
/* {{{ struc setups */

typedef struct {
    double  x, y, z;
} point;

typedef struct {
    point     pt[3];    /* Vertices of triangle */
    double    area;     /* Unused; might be used for adaptive subdivision */
} triangle;

typedef struct {
    int       npoly;    /* # of triangles in object */
    triangle *poly;     /* Triangles */
} object;

#define MAX_CONNECTIONS 10

typedef struct { 
  double x,y,z,xnew,ynew,znew,nx,ny,nz,xorig,yorig,zorig; /* co-ords, updated co-ords, estimated normal, original co-ords */
  int    n[MAX_CONNECTIONS]; /* index numbers of connected points, terminated with -1 */
} points_struc;

/* }}} */
/* {{{ vertex setups */

extern int PPHIGSflag;
extern int Flatflag;

/* Six equidistant points lying on the unit sphere */
#define XPLUS {  1,  0,  0 }    /*  X */
#define XMIN  { -1,  0,  0 }    /* -X */
#define YPLUS {  0,  1,  0 }    /*  Y */
#define YMIN  {  0, -1,  0 }    /* -Y */
#define ZPLUS {  0,  0,  1 }    /*  Z */
#define ZMIN  {  0,  0, -1 }    /* -Z */

/* Vertices of a unit octahedron */
extern triangle octahedron[];

/* A unit octahedron */
extern object oct;

/* Vertices of a tetrahedron */
#define sqrt_3 0.5773502692
#define PPP {  sqrt_3,  sqrt_3,  sqrt_3 }   /* +X, +Y, +Z */
#define MMP { -sqrt_3, -sqrt_3,  sqrt_3 }   /* -X, -Y, +Z */
#define MPM { -sqrt_3,  sqrt_3, -sqrt_3 }   /* -X, +Y, -Z */
#define PMM {  sqrt_3, -sqrt_3, -sqrt_3 }   /* +X, -Y, -Z */

/* Structure describing a tetrahedron */
extern triangle tetrahedron[];

extern object tet;

/* Twelve vertices of icosahedron on unit sphere */
#define tau 0.8506508084      /* t=(1+sqrt(5))/2, tau=t/sqrt(1+t^2)  */
#define one 0.5257311121      /* one=1/sqrt(1+t^2) , unit sphere     */
#define ZA {  tau,  one,    0 }
#define ZB { -tau,  one,    0 }
#define ZC { -tau, -one,    0 }
#define ZD {  tau, -one,    0 }
#define YA {  one,   0 ,  tau }
#define YB {  one,   0 , -tau }
#define YC { -one,   0 , -tau }
#define YD { -one,   0 ,  tau }
#define XA {   0 ,  tau,  one }
#define XB {   0 , -tau,  one }
#define XC {   0 , -tau, -one }
#define XD {   0 ,  tau, -one }

/* Structure for unit icosahedron */
extern triangle icosahedron[];

/* A unit icosahedron */
extern object ico;

/* }}} */
/* {{{ proc defs */

void tessa(int, object**);
point *normalize(point*);
point *midpoint(point*, point*);
void flip_object(object*);
void print_object(object*, int);
void print_triangle(triangle*);
void pphigs_header(int);
void pphigs_trailer();
int points_list(object*, points_struc**);
void draw_triangle(image_struct*, double, double, double, double, double, double, double, double, double, FDT);
void draw_surface(image_struct*, FDT, FDT, points_struc*, int);
void fill_surface(image_struct*, points_struc*, int);
void output_volume(points_struc*, int);
void xtopol_output(points_struc*, int, char*);

/* }}} */
