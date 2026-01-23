#include "fractions.h"
#define BGHOSTS 2
#define EMBED 1

scalar cs[];
face vector fs[];

double (* metric_embed_factor) (Point, coord) = NULL;

#if TREE
# include "embed-tree.h"
#endif

#undef SEPS
#define SEPS 1e-30

#define cs_avg(a,i,j,k)							\
  ((a[i,j,k]*(1.5 + cs[i,j,k]) + a[i-1,j,k]*(1.5 + cs[i-1,j,k]))/	\
   (cs[i,j,k] + cs[i-1,j,k] + 3.))

#if dimension == 2
#define face_condition(fs, cs)						\
  (fs.x[i,j] > 0.5 && fs.y[i,j + (j < 0)] && fs.y[i-1,j + (j < 0)] &&	\
   cs[i,j] && cs[i-1,j])

foreach_dimension()
static inline double embed_face_gradient_x (Point point, scalar a, int i)
{
  int j = sign(fs.x[i,1] - fs.x[i,-1]);
  assert (cs[i] && cs[i-1]);
  if (face_condition (fs, cs))
    return ((1. + fs.x[i])*(a[i] - a[i-1]) +
	    (1. - fs.x[i])*(a[i,j] - a[i-1,j]))/(2.*Delta);
  return (a[i] - a[i-1])/Delta;
}

foreach_dimension()
static inline double embed_face_value_x (Point point, scalar a, int i)
{
  int j = sign(fs.x[i,1] - fs.x[i,-1]);
  return face_condition (fs, cs) ?
    ((1. + fs.x[i])*cs_avg(a,i,0,0) + (1. - fs.x[i])*cs_avg(a,i,j,0))/2. :
    cs_avg(a,i,0,0);
}

#else
foreach_dimension()
static inline coord embed_face_barycentre_z (Point point, int i)
{

  coord n1 = {0};
  double nn = 0.;
  scalar f = fs.z;
  foreach_dimension(2) {
    n1.x = (f[-1,-1,i] + 2.*f[-1,0,i] + f[-1,1,i] -
	    f[+1,-1,i] - 2.*f[+1,0,i] - f[+1,1,i]);
    nn += fabs(n1.x);
  }
  if (!nn)
    return (coord){0.,0.,0.};
  foreach_dimension(2)
    n1.x /= nn;

  coord n, p1, p;
  ((double *)&n)[0] = n1.x, ((double *)&n)[1] = n1.y;
  double alpha = line_alpha (f[0,0,i], n);
  line_center (n, alpha, f[0,0,i], &p1);
  p.x = ((double *)&p1)[0], p.y = ((double *)&p1)[1], p.z = 0.;
  return p;
}

#define face_condition(fs, cs)						\
  (fs.x[i,j,k] > 0.5 && (fs.x[i,j,0] > 0.5 || fs.x[i,0,k] > 0.5) &&	\
   fs.y[i,j + (j < 0),0] && fs.y[i-1,j + (j < 0),0] &&			\
   fs.y[i,j + (j < 0),k] && fs.y[i-1,j + (j < 0),k] &&			\
   fs.z[i,0,k + (k < 0)] && fs.z[i-1,0,k + (k < 0)] &&			\
   fs.z[i,j,k + (k < 0)] && fs.z[i-1,j,k + (k < 0)] &&			\
   cs[i-1,j,0] && cs[i-1,0,k] && cs[i-1,j,k] &&				\
   cs[i,j,0] && cs[i,0,k] && cs[i,j,k])

foreach_dimension()
static inline double embed_face_gradient_x (Point point, scalar a, int i)
{
  assert (cs[i] && cs[i-1]);
  coord p = embed_face_barycentre_x (point, i);

  int j = sign(p.y), k = sign(p.z);
  if (face_condition(fs, cs)) {
    p.y = fabs(p.y), p.z = fabs(p.z);
    return (((a[i,0,0] - a[i-1,0,0])*(1. - p.y) +
	     (a[i,j,0] - a[i-1,j,0])*p.y)*(1. - p.z) +
	    ((a[i,0,k] - a[i-1,0,k])*(1. - p.y) +
	     (a[i,j,k] - a[i-1,j,k])*p.y)*p.z)/Delta;
  }
  return (a[i] - a[i-1])/Delta;
}

foreach_dimension()
static inline double embed_face_value_x (Point point, scalar a, int i)
{
  coord p = embed_face_barycentre_x (point, i);

  int j = sign(p.y), k = sign(p.z);
  if (face_condition(fs, cs)) {
    p.y = fabs(p.y), p.z = fabs(p.z);
    return ((cs_avg(a,i,0,0)*(1. - p.y) + cs_avg(a,i,j,0)*p.y)*(1. - p.z) +
	    (cs_avg(a,i,0,k)*(1. - p.y) + cs_avg(a,i,j,k)*p.y)*p.z);
  }
  return cs_avg(a,i,0,0);
}
#endif

attribute {
  bool third;
}

#undef face_gradient_x
#define face_gradient_x(a,i)					\
  (a.third && fs.x[i] < 1. && fs.x[i] > 0. ?			\
   embed_face_gradient_x (point, a, i) :			\
   (a[i] - a[i-1])/Delta)

#undef face_gradient_y
#define face_gradient_y(a,i)					\
  (a.third && fs.y[0,i] < 1. && fs.y[0,i] > 0. ?		\
   embed_face_gradient_y (point, a, i) :			\
   (a[0,i] - a[0,i-1])/Delta)

#undef face_gradient_z
#define face_gradient_z(a,i)					\
  (a.third && fs.z[0,0,i] < 1. && fs.z[0,0,i] > 0. ?		\
   embed_face_gradient_z (point, a, i) :			\
   (a[0,0,i] - a[0,0,i-1])/Delta)

#undef face_value
#define face_value(a,i)							\
  (a.third && fs.x[i] < 1. && fs.x[i] > 0. ?				\
   embed_face_value_x (point, a, i) :					\
   cs_avg(a,i,0,0))

#undef center_gradient
#define center_gradient(a) (fs.x[] && fs.x[1] ? (a[1] - a[-1])/(2.*Delta) : \
			    fs.x[1] ? (a[1] - a[])/Delta :		    \
			    fs.x[]  ? (a[] - a[-1])/Delta : 0.)

static inline
double embed_geometry (Point point, coord * p, coord * n)
{
  *n = facet_normal (point, cs, fs);
  double alpha = plane_alpha (cs[], *n);
  double area = plane_area_center (*n, alpha, p);
  normalize (n);
  return area;
}

static inline
double embed_area_center (Point point, double * x1, double * y1, double * z1)
{
  double area = 0.;
  if (cs[] > 0. && cs[] < 1.) {
    coord n, p;
    area = embed_geometry (point, &p, &n);
    *x1 += p.x*Delta, *y1 += p.y*Delta, *z1 += p.z*Delta;
  }
  return area;
}

#define embed_pos() embed_area_center (point, &x, &y, &z)

double embed_interpolate (Point point, scalar s, coord p)
{
  assert (dimension == 2);
  int i = sign(p.x), j = sign(p.y);
  if (cs[i] && cs[0,j] && cs[i,j])

    return ((s[]*(1. - fabs(p.x)) + s[i]*fabs(p.x))*(1. - fabs(p.y)) +
	    (s[0,j]*(1. - fabs(p.x)) + s[i,j]*fabs(p.x))*fabs(p.y));
  else {

    double val = s[];
    foreach_dimension() {
      int i = sign(p.x);
      if (cs[i])
	val += fabs(p.x)*(s[i] - s[]);
      else if (cs[-i])
	val += fabs(p.x)*(s[] - s[-i]);
    }
    return val;
  }
}

struct Cleanup {
  scalar c;
  face vector s;
  double smin;
  bool opposite;
};

trace
int fractions_cleanup (scalar c, face vector s,
		       double smin = 0., bool opposite = false)
{

  int changed = 1, schanged = 0, i;
  for (i = 0; i < 100 && changed; i++) {

    foreach_face()
      if (s.x[] && ((!c[] || !c[-1]) || s.x[] < smin))
	s.x[] = 0.;

    changed = 0;
    foreach(reduction(+:changed))
      if (c[] > 0. && c[] < 1.) {
	int n = 0;
	foreach_dimension() {
	  for (int i = 0; i <= 1; i++)
	    if (s.x[i] > 0.)
	      n++;

	  if (opposite && s.x[] == 0. && s.x[1] == 0.)
	    c[] = 0., changed++;
	}

	if (n < dimension)
	  c[] = 0., changed++;
      }

    schanged += changed;
  }
  if (changed)
    fprintf (stderr, "src/embed.h:%d: warning: fractions_cleanup() did not converge after "
	     "%d iterations\n", LINENO, i);
  return schanged;
}

#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))

foreach_dimension()
static inline double dirichlet_gradient_x (Point point, scalar s, scalar cs,
					   coord n, coord p, double bc,
					   double * coef)
{
  foreach_dimension()
    n.x = - n.x;
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !fs.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (fs.x[i + (i < 0),j] && fs.y[i,j] && fs.y[i,j+1] &&
	  cs[i,j-1] && cs[i,j] && cs[i,j+1])
	v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fs.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!fs.y[i,j,k+m] || !fs.y[i,j+1,k+m] ||
	    !fs.z[i,j+m,k] || !fs.z[i,j+m,k+1] ||
	    !cs[i,j+m,k-1] || !cs[i,j+m,k] || !cs[i,j+m,k+1])
	  defined = false;
      if (defined)

	v[l] =
	  quadratic (z,
		     quadratic (y1,
				(s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
		     quadratic (y1,
				(s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
		     quadratic (y1,
				(s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif
      else
	break;
    }
  if (v[0] == nodata) {

    d[0] = max(1e-3, fabs(p.x/n.x));
    *coef = - 1./(d[0]*Delta);
    return bc/(d[0]*Delta);
  }

  *coef = 0.;
  if (v[1] != nodata)
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  return (bc - v[0])/(d[0]*Delta);
}

double dirichlet_gradient (Point point, scalar s, scalar cs,
			   coord n, coord p, double bc, double * coef)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return dirichlet_gradient_x (point, s, cs, n, p, bc, coef);
#else
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return dirichlet_gradient_x (point, s, cs, n, p, bc, coef);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return dirichlet_gradient_y (point, s, cs, n, p, bc, coef);
  return dirichlet_gradient_z (point, s, cs, n, p, bc, coef);
#endif
  return nodata;
}

bid embed;

static inline
coord embed_gradient (Point point, vector u, coord p, coord n)
{
  coord dudn;
  foreach_dimension() {
    bool dirichlet = false;
    double vb = u.x.boundary[embed] (point, point, u.x, &dirichlet);
    if (dirichlet) {
      double val;
      dudn.x = dirichlet_gradient (point, u.x, cs, n, p, vb, &val);
    }
    else
      dudn.x = vb;
    if (dudn.x == nodata)
      dudn.x = 0.;
  }
  return dudn;
}

trace
void embed_force (scalar p, vector u, face vector mu, coord * Fp, coord * Fmu)
{
  coord Fps = {0}, Fmus = {0};
  foreach (reduction(+:Fps) reduction(+:Fmus), nowarning)
    if (cs[] > 0. && cs[] < 1.) {

      coord n, b;
      double area = embed_geometry (point, &b, &n);
      area *= pow (Delta, dimension - 1);
      double Fn = area*embed_interpolate (point, p, b);
      foreach_dimension()
	Fps.x += Fn*n.x;

      if (constant(mu.x) != 0.) {
	double mua = 0., fa = 0.;
	foreach_dimension() {
	  mua += mu.x[] + mu.x[1];
	  fa  += fm.x[] + fm.x[1];
	}
	mua /= fa;

	assert (dimension == 2);
	coord dudn = embed_gradient (point, u, b, n);
	foreach_dimension()
	  Fmus.x -= area*mua*(dudn.x*(sq(n.x) + 1.) + dudn.y*n.x*n.y);
      }
    }

  *Fp = Fps; *Fmu = Fmus;
}

#if dimension == 2
double embed_vorticity (Point point, vector u, coord p, coord n)
{

  coord dudn = embed_gradient (point, u, p, n);

  return dudn.y*n.x - dudn.x*n.y;
}
#endif

double embed_flux (Point point, scalar s, face vector mu, double * val)
{

  *val = 0.;
  if (cs[] >= 1. || cs[] <= 0.)
    return 0.;

  bool dirichlet = false;
  double grad = s.boundary[embed] (point, point, s, &dirichlet);
  if (!grad && !dirichlet)
    return 0.;

  coord n = facet_normal (point, cs, fs), p;
  double alpha = plane_alpha (cs[], n);
  double area = plane_area_center (n, alpha, &p);
  if (metric_embed_factor)
    area *= metric_embed_factor (point, p);

  double coef = 0.;
  if (dirichlet) {
    normalize (&n);
    grad = dirichlet_gradient (point, s, cs, n, p, grad, &coef);
  }

  double mua = 0., fa = 0.;
  foreach_dimension() {
    mua += mu.x[] + mu.x[1];
    fa  += fm.x[] + fm.x[1];
  }
  *val = - mua/(fa + SEPS)*grad*area/Delta;
  return - mua/(fa + SEPS)*coef*area/Delta;
}

macro2
double dirichlet (double expr, Point point = point,
		  scalar s = _s, bool * data = data)
{
  return data ? embed_area_center (point, &x, &y, &z),
    *((bool *)data) = true, expr : 2.*expr - s[];
}

macro2
double dirichlet_homogeneous (double expr, Point point = point,
			      scalar s = _s, bool * data = data)
{
  return data ? *((bool *)data) = true, 0 : - s[];
}

macro2
double neumann (double expr, Point point = point,
		scalar s = _s, bool * data = data)
{
  return data ? embed_area_center (point, &x, &y, &z),
    *((bool *)data) = false, expr : Delta*expr + s[];
}

macro2
double neumann_homogeneous (double expr, Point point = point,
			    scalar s = _s, bool * data = data)
{
  return data ? *((bool *)data) = false, 0 : s[];
}

#if MULTIGRID
static inline double bilinear_embed (Point point, scalar s)
{
  if (!coarse(cs) || !coarse(cs,child.x))
    return coarse(s);
  #if dimension >= 2
  if (!coarse(cs,0,child.y) || !coarse(cs,child.x,child.y))
    return coarse(s);
  #endif
  #if dimension >= 3
  if (!coarse(cs,0,0,child.z) || !coarse(cs,child.x,0,child.z) ||
      !coarse(cs,0,child.y,child.z) ||
      !coarse(cs,child.x,child.y,child.z))
    return coarse(s);
  #endif
  return bilinear (point, s);
}

#define bilinear(point, s) bilinear_embed(point, s)
#endif

trace
void update_tracer (scalar f, face vector uf, face vector flux, double dt)
{

  scalar e[];
  foreach() {

    if (cs[] <= 0.)
      e[] = 0.;

    else if (cs[] >= 1.) {
      foreach_dimension()
	f[] += dt*(flux.x[] - flux.x[1])/Delta;
      e[] = 0.;
    }

    else {
      double umax = 0.;
      for (int i = 0; i <= 1; i++)
	foreach_dimension()
	  if (fabs(uf.x[i]) > umax)
	    umax = fabs(uf.x[i]);
      double dtmax = Delta*cm[]/(umax + SEPS);

      double F = 0.;
      foreach_dimension()
	F += flux.x[] - flux.x[1];
      F /= Delta*cm[];

      if (dt <= dtmax) {
	f[] += dt*F;
	e[] = 0.;
      }

      else {
	f[] += dtmax*F;
	double scs = 0.;
	foreach_neighbor(1)
	  scs += sq(cm[]);
	e[] = (dt - dtmax)*F*cm[]/scs;
      }
    }
  }

  foreach() {
    double se = 0.;
    foreach_neighbor(1)
      se += e[];
    f[] += cs[]*se;
  }
}

event metric (i = 0)
{
  if (is_constant (fm.x)) {
    foreach_dimension()
      assert (constant (fm.x) == 1.);
    fm = fs;
  }
  foreach_face()
    fs.x[] = 1.;
  if (is_constant (cm)) {
    assert (constant (cm) == 1.);
    cm = cs;
  }
  foreach()
    cs[] = 1.;

#if TREE
  cs.refine = embed_fraction_refine;

  cs.prolongation = fraction_refine;
  foreach_dimension()
    fs.x.prolongation = embed_face_fraction_refine_x;

#endif
  restriction ({cs, fs});
}

event defaults (i = 0 ) {
  display ("draw_vof (c = 'cs', s = 'fs', filled = -1, "
	   "fc = {0.5,0.5,0.5}, order = 2);");
}
