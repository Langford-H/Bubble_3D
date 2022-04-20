#include "grid/octree.h"
#include "navier-stokes/centered.h"

#include "navier-stokes/perfs.h"

#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"

#include "tension.h"


#if dimension == 3
  #include "lambda2.h"
#endif

#include "view.h"
#include "maxruntime.h"

/**
The density ratio is 1000 and the dynamic viscosity ratio 100. */

#define RHOR 1000.
#define MUR 100.

# define Ga 100.
# define Bo 0.12
# define MAXTIME 110


#define WIDTH 20.
#define R0 0.5

int LEVEL = 10;

u.t[bottom]=dirichlet(0);
u.n[top]=0;

/**
The main function can take two optional parameters: the maximum level
of adaptive refinement (as well as an optional maximum runtime). */


void bubble_position_setting_3D (scalar f)//向本函数输入两个空的scalar数据，函数将其填充
{
    vertex scalar phi[];//先定义一个level-set函数
    foreach_vertex() {//共32个气泡，对每一个单元格都针对32个气泡的levelset函数进行一次迭代，以确保其边角上的函数完整
    phi[] = HUGE;
        for (double xp = -5. ; xp <= 5.; xp += 10.)//此处i、j为气泡的圆心，每一次只更迭一个气泡
            for (double yp = -8. ; yp <= -4.  ; yp += 2.)
            	for (double zp = -5. ; zp <= 5.; zp += 10.)
               	phi[] = intersection (phi[], (sq(x - xp) + sq(y - yp) + sq(z - zp) - sq(R0)));//将每一次结果都重复输入进phi[]中

    }//着重需要注意intersection()函数的插入机制比较奇怪，由于并没有单独列出经过调试猜测应该是在每一次遍历时取小值，是故必须将levelset函数先设置为气泡外部为正，气泡内部为负值，然后反向，否则会出现全部都取为负数的情况
    refine( phi[] < 0 && level < LEVEL);
    fractions (phi, f);//利用上文中得到的phi[]将其带入fractions函数，并将两个scalar型函数填满
}

int main (int argc, char * argv[]) {
  maxruntime (&argc, argv);
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  
  /**
  We set the domain geometry and initial refinement. */
  
  size (WIDTH);
  origin (-10., -10., -10.);
  init_grid (128);

  /**
  We set the physical parameters: densities, viscosities and surface
  tension. */
  
  rho1 = 1.;
  rho2 = 1./RHOR;
  mu1 = 1./Ga;
  mu2 = 1./(MUR*Ga);
  f.sigma = 1./Bo;
  
  periodic(left);
  periodic(front);

  TOLERANCE = 1e-4;
  run();
}

/**
For the initial conditions, we first try to restore the simulation
from a previous "restart", if this fails we refine the mesh locally to
the maximum level, in a sphere of diameter 1.5 around the bubble. We
then initialise the volume fraction for a bubble initially at (0,Z0,0)
of diameter unity. */

event init (t = 0) {
  if (!restore (file = "restart")) {
	bubble_position_setting_3D(f);
	//porous(f);
	dump();
  }
}

/**
We add the acceleration of gravity (unity) in the downward (-y)
direction. */

event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= 1.;
}

/**
We adapt the mesh by controlling the error on the volume fraction and
velocity field. */

event adapt (i++) {
  double uemax = 1e-2;
  adapt_wavelet ({f,u}, (double[]){0.01,uemax,uemax,uemax}, LEVEL, 5);
}

/**
## Outputs

Every ten timesteps, we output the time, volume, position, and
velocity of the bubble. */

event logfile (i = 0 ; t <= MAXTIME; i += 10) {
  double xb = 0., yb = 0., zb = 0., sb = 0.;
  double vbx = 0., vby = 0., vbz = 0.;
  foreach(reduction(+:xb) reduction(+:yb) reduction(+:zb)
	  reduction(+:vbx) reduction(+:vby) reduction(+:vbz)
	  reduction(+:sb)) {
    double dv = f[]*dv();
    xb += x*dv;
    yb += y*dv;
    zb += z*dv;
    vbx += u.x[]*dv;
    vby += u.y[]*dv;
    vbz += u.z[]*dv;
    sb += dv;
  }
  fprintf (stderr,
	   "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", 
	   t, sb,
	   xb/sb, yb/sb, zb/sb,
	   vbx/sb, vby/sb, vbz/sb);
  fflush (stderr);
}

/*
event snapshot (t = 1; t <= MAXTIME; t++)
{
  scalar l2[], omegay[];
  lambda2 (u, l2);
  foreach()
    omegay[] = (u.z[1] - u.z[-1] - u.x[0,0,1] + u.x[0,0,-1])/(2.*Delta);
  boundary ({omegay});
  
  char name[80];
  sprintf (name, "dump-%03d", (int) t);
  dump (file = name);
}
*/
