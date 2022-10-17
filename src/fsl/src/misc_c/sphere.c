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
/* {{{ includes */

#include "libss/libss.h"
#include "libss/libavw.h"
#include "libss/libtessa.h"

/* }}} */
/* {{{ main */

int main(ac, av)
     int ac;
     char *av[];
{
    object *old = &oct;         /* Default is octahedron */
    int     ccwflag = 0,        /* Reverse vertex order if true */
            i, pc,
            normal = 0,
            xtopol = 0,
            volume = 0,
            maxlevel = 1;       /* Maximum subdivision level */
    points_struc *points;

    /* {{{ Parse arguments */

    for (i = 1; i < ac; i++) {
        if (!strcmp(av[i], "-p"))
            PPHIGSflag = 1;
        else if (!strcmp(av[i], "-c"))
            ccwflag = 1;
        else if (!strcmp(av[i], "-v"))
            volume = 1;
        else if (!strcmp(av[i], "-n"))
            normal = 1;
        else if (!strcmp(av[i], "-x"))
            xtopol = 1;
        else if (!strcmp(av[i], "-f"))
            Flatflag = 1;
        else if (!strcmp(av[i], "-t"))
            old = &tet;
        else if (!strcmp(av[i], "-i"))
            old = &ico;
        else if (isdigit((int)av[i][0])) {
            if ((maxlevel = atoi(av[i])) < 1) {
                fprintf(stderr, "%s: # of levels must be >= 1\n", av[0]);
                exit(1);
            }
        } else {
            fprintf(stderr, "Usage: %s [-p] [-c] [-f] [-t] [-i] [-v] [-x] [-n] [n]\n", av[0]);
            exit(1);
        }
}

/* }}} */

    if (ccwflag)
        flip_object(old);

    tessa(maxlevel,&old);

    pc=points_list(old,&points);

    if (xtopol)
      xtopol_output(points,pc,"grot");

    if (volume)
      output_volume(points,pc);

    if (normal)
      print_object(old, maxlevel);

    return(0);
}

/* }}} */
