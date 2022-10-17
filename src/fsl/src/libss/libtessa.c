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

#include "libss.h"
#include "libavw.h"
#include "libtessa.h"

/* }}} */
/* {{{ globals */

int PPHIGSflag = 0; /* Don't generate PPHIGS format output */
int Flatflag = 0;   /* Don't generate per-vertex normals */

triangle octahedron[] = {
    { { XPLUS, ZPLUS, YPLUS }, 0.0 },
    { { YPLUS, ZPLUS, XMIN  }, 0.0 },
    { { XMIN , ZPLUS, YMIN  }, 0.0 },
    { { YMIN , ZPLUS, XPLUS }, 0.0 },
    { { XPLUS, YPLUS, ZMIN  }, 0.0 },
    { { YPLUS, XMIN , ZMIN  }, 0.0 },
    { { XMIN , YMIN , ZMIN  }, 0.0 },
    { { YMIN , XPLUS, ZMIN  }, 0.0 }
};

object oct = {
    sizeof(octahedron) / sizeof(octahedron[0]),
    &octahedron[0]
};

triangle tetrahedron[] = {
    { { PPP, MMP, MPM }, 0.0 },
    { { PPP, PMM, MMP }, 0.0 },
    { { MPM, MMP, PMM }, 0.0 },
    { { PMM, PPP, MPM }, 0.0 }
};

object tet = {
    sizeof(tetrahedron) / sizeof(tetrahedron[0]),
    &tetrahedron[0]
};

triangle icosahedron[] = {
    { { YA, XA, YD }, 0.0 },
    { { YA, YD, XB }, 0.0 },
    { { YB, YC, XD }, 0.0 },
    { { YB, XC, YC }, 0.0 },
    { { ZA, YA, ZD }, 0.0 },
    { { ZA, ZD, YB }, 0.0 },
    { { ZC, YD, ZB }, 0.0 },
    { { ZC, ZB, YC }, 0.0 },
    { { XA, ZA, XD }, 0.0 },
    { { XA, XD, ZB }, 0.0 },
    { { XB, XC, ZD }, 0.0 },
    { { XB, ZC, XC }, 0.0 },
    { { XA, YA, ZA }, 0.0 },
    { { XD, ZA, YB }, 0.0 },
    { { YA, XB, ZD }, 0.0 },
    { { YB, ZD, XC }, 0.0 },
    { { YD, XA, ZB }, 0.0 },
    { { YC, ZB, XD }, 0.0 },
    { { YD, ZC, XB }, 0.0 },
    { { YC, XC, ZC }, 0.0 }
};

object ico = {
    sizeof(icosahedron) / sizeof(icosahedron[0]),
    &icosahedron[0]
};

/* }}} */
/* {{{ tessa */

void tessa(int maxlevel, object **old)
{
  int level, i;
  object *new;

    for (level = 1; level < maxlevel; level++) {
        /* Allocate a new object */
        new = (object *)malloc(sizeof(object));
        if (new == NULL) {
            fprintf(stderr, "Out of memory on subdivision level %d\n",level);
            exit(1);
        }
        new->npoly = (*old)->npoly * 4;

        /* Allocate 4* the number of points in the current approximation */
        new->poly  = (triangle *)malloc(new->npoly * sizeof(triangle));
        if (new->poly == NULL) {
            fprintf(stderr, "Out of memory on subdivision level %d\n",level);
            exit(1);
        }

        /* Subdivide each triangle in the old approximation and normalize
         *  the new points thus generated to lie on the surface of the unit
         *  sphere.
         * Each input triangle with vertices labelled [0,1,2] as shown
         *  below will be turned into four new triangles:
         *
         *                      Make new points
         *                          a = (0+2)/2
         *                          b = (0+1)/2
         *                          c = (1+2)/2
         *        1
         *       /\             Normalize a, b, c
         *      /  \
         *    b/____\ c         Construct new triangles
         *    /\    /\              [0,b,a]
         *   /  \  /  \             [b,1,c]
         *  /____\/____\            [a,b,c]
         * 0      a     2           [a,c,2]
         */
        for (i = 0; i < (*old)->npoly; i++) {
            triangle
                 *oldt = &(*old)->poly[i],
                 *newt = &new->poly[i*4];
            point a, b, c;

            a = *normalize(midpoint(&oldt->pt[0], &oldt->pt[2]));
            b = *normalize(midpoint(&oldt->pt[0], &oldt->pt[1]));
            c = *normalize(midpoint(&oldt->pt[1], &oldt->pt[2]));

            newt->pt[0] = oldt->pt[0];
            newt->pt[1] = b;
            newt->pt[2] = a;
            newt++;

            newt->pt[0] = b;
            newt->pt[1] = oldt->pt[1];
            newt->pt[2] = c;
            newt++;

            newt->pt[0] = a;
            newt->pt[1] = b;
            newt->pt[2] = c;
            newt++;

            newt->pt[0] = a;
            newt->pt[1] = c;
            newt->pt[2] = oldt->pt[2];
        }

	if (level > 1) {
	  free((*old)->poly);
	  free((*old));
	}

        /* Continue subdividing new triangles */
        (*old) = new;
    }
}

/* }}} */
/* {{{ normalize */

point *normalize(point *p)
{
    static point r;
    double mag;

    r = *p;
    mag = r.x * r.x + r.y * r.y + r.z * r.z;
    if (mag != 0.0) {
        mag = 1.0 / sqrt(mag);
        r.x *= mag;
        r.y *= mag;
        r.z *= mag;
    }

    return &r;
}

/* }}} */
/* {{{ midpoint */

/* Return the midpoint on the line between two points */
point *midpoint(point *a, point *b)
{
    static point r;

    r.x = (a->x + b->x) * 0.5;
    r.y = (a->y + b->y) * 0.5;
    r.z = (a->z + b->z) * 0.5;

    return &r;
}

/* }}} */
/* {{{ flip_object */

/* Reverse order of points in each triangle */
void flip_object(object *obj)
{
    int i;
    for (i = 0; i < obj->npoly; i++) {
        point tmp;
                       tmp = obj->poly[i].pt[0];
        obj->poly[i].pt[0] = obj->poly[i].pt[2];
        obj->poly[i].pt[2] = tmp;
    }
}

/* }}} */
/* {{{ print_object */

void print_object(object *obj, int level)
{
    int i;

    if (PPHIGSflag)
        pphigs_header(level);

    /* Spit out coordinates for each triangle */
    for (i = 0; i < obj->npoly; i++)
        print_triangle(&obj->poly[i]);

    if (PPHIGSflag)
        pphigs_trailer();
}

/* }}} */
/* {{{ print_triangle */

void print_triangle(triangle *t)
{
    int i;

    if (PPHIGSflag) {
        printf("\tpolygon 3 {\n");
        for (i = 0; i < 3; i++)
            if (Flatflag) {
                printf("\t\t%g\t%g\t%g ;\n",
                    t->pt[i].x, t->pt[i].y, t->pt[i].z);    /* Point */
            } else {
                printf("\t\t%g\t%g\t%g %g\t%g\t%g ;\n",
                    t->pt[i].x, t->pt[i].y, t->pt[i].z,     /* Point */
                    t->pt[i].x, t->pt[i].y, t->pt[i].z);    /* Normal */
            }
        printf("\t};\n");
    } else {
        /* Modify this to generate your favorite output format
         * Triangle vertices are in t->pt[0..2].{x,y,z}
         * A generic format is provided.
         */
        printf("triangle\n");
        for (i = 0; i < 3; i++)
            printf("\t%g %g %g\n", t->pt[i].x, t->pt[i].y, t->pt[i].z);
    }
}

/* }}} */
/* {{{ phigs */

void pphigs_header(int level)
{
    int dx, dy, dz;

    if (Flatflag)
        printf("structure flat%d posted {\n", level);
    else
        printf("structure sphere%d posted {\n", level);
    printf("\tcolor polygon {\n");
    printf("\t\t200 100  50   0     50 100 200   0\n");
    printf("\t};\n");

    switch (level) {
        case 1:
            dx = -2000; dy =  2000; dz = 0;
            break;
        case 2:
            dx =  2000; dy =  2000; dz = 0;
            break;
        case 3:
            dx = -2000; dy = -2000; dz = 0;
            break;
        case 4:
            dx =  2000; dy = -2000; dz = 0;
            break;
        case 5:
            dx =     0; dy =     0; dz = 0;
            break;
        default:
            dx = dy = dz = 0;
            break;
    }

    printf("\tmatrix Pre scale 1000 1000 1000;\n");
    printf("\tmatrix Pre translate %d %d %d ;\n", dx, dy, dz);
}

void pphigs_trailer()
{
    printf("};\n");
}

/* }}} */
/* {{{ points_list */

#define MISS 0.0001          /* to test if points are identical, this is max error */

int points_list(object *obj, points_struc **points)
{
  int i,j,k,pc=0,pointentry[3];

  *points = (points_struc *)malloc(sizeof(points_struc)*(obj->npoly+50)); /* Yes, I _could_ predict the number of points exactly */

  for (i = 0; i < obj->npoly; i++)
    {
      triangle *lp = &obj->poly[i];

      /* {{{ phase 1 - loop round triangle vertices finding/creating positions in points list */

      for (j = 0; j < 3; j++)
	{
	  for(k=0; k<pc; k++) /* test to see if this vertex is already in the points list */
	      if ( (fabs(lp->pt[j].x-(*points)[k].x)<MISS) &&
		   (fabs(lp->pt[j].y-(*points)[k].y)<MISS) &&
		   (fabs(lp->pt[j].z-(*points)[k].z)<MISS) )
		break;

	  if (k==pc) /* true if the vertex wasn't already in the points list - thus setup a new entry */
	    {
	      (*points)[k].n[0]=-1; /* terminate new list of connections */
	      (*points)[k].x=lp->pt[j].x;
	      (*points)[k].y=lp->pt[j].y;
	      (*points)[k].z=lp->pt[j].z;
	      pc++;
	    }

	  pointentry[j]=k; /* remember the index number for this vertex for phase 2 */

 	  /* printf("%d %d\n",i,k); */
	}

/* }}} */
      /* {{{ phase 2 - setup the point connections */

      for (j = 0; j < 3; j++)
	{
	  int l,m;

	  for(k=0; (*points)[pointentry[j]].n[k]>=0; k++); /* find end of connections list */

	  for(m=1; m<3; m++) /* (j+m)%3 points to the two other vertices */
	    {
	      for(l=0; l<k; l++) /* has connection to vertex (j+m)%3 already been made? */
		if ( (fabs(lp->pt[(j+m)%3].x-(*points)[(*points)[pointentry[j]].n[l]].x)<MISS) &&
		     (fabs(lp->pt[(j+m)%3].y-(*points)[(*points)[pointentry[j]].n[l]].y)<MISS) &&
		     (fabs(lp->pt[(j+m)%3].z-(*points)[(*points)[pointentry[j]].n[l]].z)<MISS) )
		  l=1000;
	      if (l<1000) /* make a new connection */
		{
		  (*points)[pointentry[j]].n[k]=pointentry[(j+m)%3];
		  k++;
		}
	    }

	  (*points)[pointentry[j]].n[k]=-1; /* re-terminate connections list */
	}

/* }}} */
    }

  /* {{{ phase 3 - reorder connections */

for(i=0; i<pc; i++)
{
  int l,m;

  for(k=0; (*points)[i].n[k]>=0; k++); /* k = number of connections */

  for(l=0; l<k-1; l++)
    {
      double adx, ady, adz, bdx, bdy, bdz, tmpf, nx, ny, nz, sint, cost, theta, min_theta=10, min_m=0;
      
      adx = (*points)[(*points)[i].n[l]].x - (*points)[i].x;
      ady = (*points)[(*points)[i].n[l]].y - (*points)[i].y;
      adz = (*points)[(*points)[i].n[l]].z - (*points)[i].z;
      tmpf = sqrt(adx*adx+ady*ady+adz*adz);
      adx/=tmpf; ady/=tmpf; adz/=tmpf;

      /* {{{ set m to next corrent connection */

      for(m=l+1; m<k; m++)
	{
	  bdx = (*points)[(*points)[i].n[m]].x - (*points)[i].x;
	  bdy = (*points)[(*points)[i].n[m]].y - (*points)[i].y;
	  bdz = (*points)[(*points)[i].n[m]].z - (*points)[i].z;
	  tmpf = sqrt(bdx*bdx+bdy*bdy+bdz*bdz);
	  bdx/=tmpf; bdy/=tmpf; bdz/=tmpf;
	  
	  nx = ady*bdz - adz*bdy;
	  ny = adz*bdx - adx*bdz;
	  nz = adx*bdy - ady*bdx;
	  
	  sint = sqrt ( (nx*nx) + (ny*ny) + (nz*nz) );
	  if ( (nx*(*points)[i].x) + (ny*(*points)[i].y) + (nz*(*points)[i].z) < 0 )
	    sint = -sint;
	  
	  cost = (adx*bdx)+(ady*bdy)+(adz*bdz);
	  
	  theta = atan2(sint,cost);
	  if (theta<0)
	    theta += 2.0 * M_PI;
	  
	  if (theta<min_theta)
	    {
	      min_theta=theta;
	      min_m=m;
	    }
	}

m=min_m;

/* }}} */

      if (m!=l+1) /* swap connections */
	{
	  int tmp;

	  tmp = (*points)[i].n[l+1];
	  (*points)[i].n[l+1] = (*points)[i].n[m];
	  (*points)[i].n[m] = tmp;
	}
    }
}

/* }}} */
  /* {{{ phase 4 - make orig co-ords */

for(i=0; i<pc; i++)
{
  (*points)[i].xorig=(*points)[i].x;
  (*points)[i].yorig=(*points)[i].y;
  (*points)[i].zorig=(*points)[i].z;
}

/* }}} */

  /* {{{ COMMENT full debug printout */

#ifdef FoldingComment

for(i=0; i<pc; i++)
{
  printf("%d %f %f %f",i,(*points)[i].x,(*points)[i].y,(*points)[i].z);
  for(j=0; j<MAX_CONNECTIONS; j++)
    printf(" %d",(*points)[i].n[j]);
  printf("\n");
}

#endif

/* }}} */

  return(pc);
}

/* }}} */
/* {{{ draw_triangle */

/*
_x is vector;  _^x is unit vector
vertices: _A, _B, _C
sides: _b = _B - _A , _c = _C - _A
general position = _A + beta * _^b + gamma * _^c
beta range: 0 -> |_b|
gamma range: 0 -> |_c| * (|_b|-beta) / |_b|
*/

/* intputs in voxels not mm */

void draw_triangle(image_struct *im, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, FDT colour)
{
  FDT *image=im->i;
  double bx=x2-x1, by=y2-y1, bz=z2-z1,
    cx=x3-x1, cy=y3-y1, cz=z3-z1,
    modb=sqrt(bx*bx+by*by+bz*bz),
    modc=sqrt(cx*cx+cy*cy+cz*cz),
    beta, gamma;
  int X,Y,Z, x_size=im->x, y_size=im->y, z_size=im->z;

  bx/=modb; by/=modb; bz/=modb;
  cx/=modc; cy/=modc; cz/=modc;

  for(beta=0; beta<=modb; beta+=0.5)  /* increments of 0.5 - seems to be necessary to make sure we have no gaps */
    for(gamma=0; gamma<=modc*(modb-beta)/modb; gamma+=0.5)
      {
	X=(int)(x1+beta*bx+gamma*cx);
	Y=(int)(y1+beta*by+gamma*cy);
	Z=(int)(z1+beta*bz+gamma*cz);
	if ( (X>=0) && (X<x_size) && (Y>=0) && (Y<y_size) && (Z>=0) && (Z<z_size) )
	  image[(Z*y_size+Y)*x_size+X]=colour;
      }
}

/* }}} */
/* {{{ draw_surface */

void draw_surface(image_struct *im, FDT colour, FDT vcolour, points_struc *points, int pc)
{
  FDT *image=im->i;
  int i, x_size=im->x, y_size=im->y, z_size=im->z;

  /* {{{ draw surface */

  for(i=0; i<pc; i++)
    {
      int k, l;
      
      for(k=0; points[i].n[k]>-1; k++); /* find number of connections */
      
      for(l=0; l<k; l++)
	if ( (i<points[i].n[l]) && (i<points[i].n[(l+1)%k]) )
	  draw_triangle(im,
			points[i].x/im->xv,points[i].y/im->yv,points[i].z/im->zv,
			points[points[i].n[l]].x/im->xv,points[points[i].n[l]].y/im->yv,points[points[i].n[l]].z/im->zv,
			points[points[i].n[(l+1)%k]].x/im->xv,points[points[i].n[(l+1)%k]].y/im->yv,points[points[i].n[(l+1)%k]].z/im->zv,
			colour);
    }

/* }}} */

  if (vcolour!=colour)
    /* {{{ draw vertices */

{
  for(i=0; i<pc; i++)
    {
      int X=(points[i].x/im->xv), Y=(points[i].y/im->yv), Z=(points[i].z/im->zv);

      if ( (X>=0) && (X<x_size) && (Y>=0) && (Y<y_size) && (Z>=0) && (Z<z_size) )
	image[(Z*y_size+Y)*x_size+X]=vcolour;
    }
}

/* }}} */
}

/* }}} */
/* {{{ COMMENT fill_surface ORIG */

#ifdef FoldingComment

void fill_surface(image_struct *im, points_struc *points, int pc)
{
  FDT *image=im->i;
  int i, x0, y0, z0, x, y, z, x_size=im->x, y_size=im->y, z_size=im->z;

  memset(image,0,sizeof(FDT)*x_size*y_size*z_size);

  draw_surface(im,(FDT)1,(FDT)1,points,pc);

  /* {{{ COMMENT find seed points x0, y0, z0 using surface normals */

#ifdef FoldingComment

for(i=0; i<pc; i++) /* try all vertices until one is inside images _and_ gives a good seed point */
{
  /* {{{ find approximate local surface normal */

  double tmpf, nx, ny, nz,
    adx = points[points[i].n[0]].x - points[i].x,
    ady = points[points[i].n[0]].y - points[i].y,
    adz = points[points[i].n[0]].z - points[i].z,
    bdx = points[points[i].n[1]].x - points[i].x,
    bdy = points[points[i].n[1]].y - points[i].y,
    bdz = points[points[i].n[1]].z - points[i].z;
  
  nx = ady*bdz - adz*bdy;
  ny = adz*bdx - adx*bdz;
  nz = adx*bdy - ady*bdx;
  tmpf = sqrt(nx*nx+ny*ny+nz*nz);
  nx/=(double)tmpf; ny/=(double)tmpf; nz/=(double)tmpf;

/* }}} */

  x0=(points[i].x/im->xv); y0=(points[i].y/im->yv); z0=(points[i].z/im->zv); /* start at a surface point */
  if ( (x0>=0) && (x0<x_size) && (y0>=0) && (y0<y_size) && (z0>=0) && (z0<z_size) )
    {
      for(tmpf=1.0; tmpf>0; tmpf+=0.1) /* slowly move inwards */
	{
	  x0=((points[i].x-nx*tmpf)/im->xv);
	  y0=((points[i].y-ny*tmpf)/im->yv);
	  z0=((points[i].z-nz*tmpf)/im->zv);
	  if ( (x0>=0) && (x0<x_size) && (y0>=0) && (y0<y_size) && (z0>=0) && (z0<z_size) &&
	       (image[(z0*y_size+y0)*x_size+x0]==0) )  /* found a valid seed point (probably!) */
	    {
	      tmpf=-1;
	      i=pc;
	    }
	}
    }
}

printf("Seed point=(%d,%d,%d)\n",x0,y0,z0);

#endif

/* }}} */
  /* {{{ find seed points x0, y0, z0 using CofG */

{
  double sx, sy, sz, sn;

  sx=sy=sz=sn=0;

  for(i=0; i<pc; i++)
    {
      sx += points[i].x;
      sy += points[i].y;
      sz += points[i].z;
      sn++;
    }

  x0 = sx/(sn*im->xv);
  y0 = sy/(sn*im->yv);
  z0 = sz/(sn*im->zv);

  /*printf("Seed point=(%d,%d,%d)\n",x0,y0,z0);*/
}

/* }}} */
  /* {{{ do the filling */

#define I(x,y,z) image[((z)*y_size+(y))*x_size+(x)]

I(x0,y0,z0)=2;

  for(z=z0; z<z_size; z++)
    for(y=0; y<y_size; y++)
      for(x=0; x<x_size; x++)
	{
	  if(I(x,y,z)==2)
	    {
	      int decx=0, decy=0, decz=0;

	      I(x,y,z)=1;

	      /* {{{ test 6-connected neighbourhood for included points */

{
  FDT *p=&I(x,y,z);
  int A=x_size*y_size;

  if( (x>0) && (*(p-1)==0) ) { *(p-1)=2; decx=1; }
  if( (x<x_size-1) && (*(p+1)==0) ) *(p+1)=2;
  if( (y>0) && (*(p-x_size)==0) ) { *(p-x_size)=2; decy=1; }
  if( (y<y_size-1) && (*(p+x_size)==0) ) *(p+x_size)=2;
  if( (z>0) && (*(p-A)==0) ) { *(p-A)=2; decz=1; }
  if( (z<z_size-1) && (*(p+A)==0) ) *(p+A)=2;
}

/* }}} */
	      /* {{{ apply decrements, for speed, if allowed */

if (decz)
{
     z--;
     x--;
}
else
{
  if (decy)
    {
      y--;
      x--;
    }
  else
    if (decx) x-=2;
}

/* }}} */
	    }
	}

/* }}} */
}

#endif

/* }}} */
/* {{{ fill_surface */

void fill_surface(image_struct *im, points_struc *points, int pc)
{
  FDT *image=im->i;
  int i, x0, y0, z0, x, y, z, x_size=im->x, y_size=im->y, z_size=im->z;

  memset(image,0,sizeof(FDT)*x_size*y_size*z_size);

  draw_surface(im,(FDT)1,(FDT)1,points,pc);

  /* {{{ find seed points x0, y0, z0 using CofG */

{
  double sx, sy, sz, sn;

  sx=sy=sz=sn=0;

  for(i=0; i<pc; i++)
    {
      sx += points[i].x;
      sy += points[i].y;
      sz += points[i].z;
      sn++;
    }

  x0 = sx/(sn*im->xv);
  y0 = sy/(sn*im->yv);
  z0 = sz/(sn*im->zv);

  /*printf("Seed point=(%d,%d,%d)\n",x0,y0,z0);*/
}

/* }}} */
  /* {{{ do the filling */

#define I(x,y,z) image[((z)*y_size+(y))*x_size+(x)]

I(x0,y0,z0)=2;

  for(z=z0; z<z_size; z++)
    for(y=0; y<y_size; y++)
      for(x=0; x<x_size; x++)
	{
	  if(I(x,y,z)==2)
	    {
	      int decx=0, decy=0, decz=0;

	      I(x,y,z)=1;

	      /* {{{ test 6-connected neighbourhood for included points */

{
  FDT *p=&I(x,y,z);
  int A=x_size*y_size;

  if( (x>0) && (*(p-1)==0) ) { *(p-1)=2; decx=1; }
  if( (x<x_size-1) && (*(p+1)==0) ) *(p+1)=2;
  if( (y>0) && (*(p-x_size)==0) ) { *(p-x_size)=2; decy=1; }
  if( (y<y_size-1) && (*(p+x_size)==0) ) *(p+x_size)=2;
  if( (z>0) && (*(p-A)==0) ) { *(p-A)=2; decz=1; }
  if( (z<z_size-1) && (*(p+A)==0) ) *(p+A)=2;
}

/* }}} */
	      /* {{{ apply decrements, for speed, if allowed */

if (decz)
{
     z--;
     x--;
}
else
{
  if (decy)
    {
      y--;
      x--;
    }
  else
    if (decx) x-=2;
}

/* }}} */
	    }
	}

/* }}} */
}

/* }}} */
/* {{{ output_volume */

#define LENGTH 256

/* create image with vertices shown */

void output_volume(points_struc *points, int pc)
{
  int i;
  image_struct im;

  im.x=LENGTH; im.y=LENGTH; im.z=LENGTH; im.t=1;
  init_image_struct(&im);

  for(i=0; i<pc; i++)
    {
      points[i].x = points[i].x*(0.45*LENGTH) + (0.5*LENGTH);
      points[i].y = points[i].y*(0.45*LENGTH) + (0.5*LENGTH);
      points[i].z = points[i].z*(0.45*LENGTH) + (0.5*LENGTH);
    }

  fill_surface(&im,points,pc);
  draw_surface(&im,(FDT)128,(FDT)255,points,pc);

  avw_write("grot",im);
}

/* }}} */
/* {{{ xtopol_output */

void xtopol_output(points_struc *points, int pc, char *fileroot)
{
  FILE *fd1, *fd2;
  char *thestring,filename[1000];
  int i,j,max_count=pc;

  sprintf(filename,"%s.coo",fileroot);
  fd1=fopen(filename,"wb");
  sprintf(filename,"%s.dat",fileroot);
  fd2=fopen(filename,"wb");

  thestring=(char*)malloc(max_count+1);
  thestring[max_count]=0;

  for (i = 0; i < pc; i++)
    {
      /*      fprintf(fd1,"%d %f %f %f\n",i,points[i].x,points[i].y,points[i].z);*/
      fprintf(fd1,". %f %f %f\n",points[i].x,points[i].y,points[i].z);
      memset(thestring,'0',max_count);
      for(j=0; points[i].n[j]>-1; j++)
	{
	  /* only make the link one way */
	  if (points[i].n[j]<i)
	    thestring[points[i].n[j]]='1';
	}
      fprintf(fd2,"%s\n",thestring);
    }
}

/* }}} */
