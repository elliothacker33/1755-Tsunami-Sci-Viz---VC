/**
# von Karman vortex street with ParaView output

Instrumented Basilisk case for flow around a cylinder at Re=160.
Outputs scalar fields to VTK legacy format and a PVD time series
for ParaView.
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "vtk.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

scalar f[];
scalar * tracers = {f};
scalar omega[];

face vector muv[];

double Reynolds = 160.;
int maxlevel = 9;
int base_n = 512;
int output_n = 512;
double output_dt = 0.1;
double t_end = 15.;
int output_movies = 0;

char out_dir[256] = "data/M";
char out_prefix[64] = "karman";
char results_dir[256] = "results";

static char pvd_path[512];
static int output_index = 0;

double D = 0.125, U0 = 1.;

static int env_int (const char *name, int def)
{
  const char *val = getenv (name);
  if (!val || !*val)
    return def;
  return atoi (val);
}

static double env_double (const char *name, double def)
{
  const char *val = getenv (name);
  if (!val || !*val)
    return def;
  return atof (val);
}

static void env_string (const char *name, const char *def,
                        char *out, int out_len)
{
  const char *val = getenv (name);
  if (!val || !*val)
    val = def;
  strncpy (out, val, out_len - 1);
  out[out_len - 1] = '\0';
}

static void ensure_dir (const char *path)
{
  char tmp[512];
  int len = (int) strlen (path);

  if (len == 0 || len >= sizeof (tmp))
    return;

  snprintf (tmp, sizeof (tmp), "%s", path);
  for (char *p = tmp + 1; *p; p++) {
    if (*p == '/') {
      *p = '\0';
      if (mkdir (tmp, 0775) && errno != EEXIST)
        return;
      *p = '/';
    }
  }
  mkdir (tmp, 0775);
}

static void read_env (void)
{
  Reynolds = env_double ("REYNOLDS", Reynolds);
  maxlevel = env_int ("MAXLEVEL", maxlevel);
  base_n = env_int ("BASE_N", base_n);
  output_n = env_int ("OUTPUT_N", base_n);
  output_dt = env_double ("OUTPUT_DT", output_dt);
  t_end = env_double ("T_END", t_end);
  output_movies = env_int ("OUTPUT_MOVIES", output_movies);

  env_string ("OUT_DIR", "data/M", out_dir, sizeof (out_dir));
  env_string ("OUT_PREFIX", "karman", out_prefix, sizeof (out_prefix));
  env_string ("RESULTS_DIR", "results", results_dir, sizeof (results_dir));
}

int main()
{
  read_env();

  L0 = 8. [1];
  origin (-0.5, -L0/2.);
  N = base_n;
  mu = muv;

  display_control (Reynolds, 10, 1000);
  display_control (maxlevel, 6, 12);

  run();
}

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*D*U0/Reynolds;
}

u.n[left]  = dirichlet(U0);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);

event init (t = 0)
{
  solid (cs, fs, intersection (intersection (0.5 - y, 0.5 + y),
                               sqrt(sq(x) + sq(y)) - D/2.));

  foreach()
    u.x[] = cs[] ? U0 : 0.;

  ensure_dir (out_dir);
  ensure_dir (results_dir);
  if (pid() == 0) {
    snprintf (pvd_path, sizeof (pvd_path), "%s/%s.pvd", out_dir, out_prefix);
    FILE * fp = fopen (pvd_path, "w");
    if (fp) {
      fputs ("<?xml version=\"1.0\"?>\n"
             "<VTKFile type=\"Collection\" version=\"0.1\" "
             "byte_order=\"LittleEndian\">\n"
             "  <Collection>\n", fp);
      fclose (fp);
    }
  }
}

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

event output (t += output_dt; t <= t_end)
{
  vorticity (u, omega);

  char vtk_name[512];
  snprintf (vtk_name, sizeof (vtk_name), "%s/%s_%06d.vtk",
            out_dir, out_prefix, output_index);
  FILE * fp = fopen (vtk_name, "w");
  if (fp) {
    scalar * list = {p, f, omega, cs, u.x, u.y};
    output_vtk (list, output_n, fp, true);
    fclose (fp);
  }

  if (pid() == 0) {
    FILE * pvd = fopen (pvd_path, "a");
    if (pvd) {
      fprintf (pvd,
               "    <DataSet timestep=\"%g\" file=\"%s_%06d.vtk\"/>\n",
               t, out_prefix, output_index);
      fclose (pvd);
    }
  }

  output_index++;
}

event movies (i += 4; t <= t_end)
{
  if (!output_movies)
    return;

  vorticity (u, omega);

  char vort_path[512];
  char tracer_path[512];
  snprintf (vort_path, sizeof (vort_path), "%s/%s_vort.mp4",
            results_dir, out_prefix);
  snprintf (tracer_path, sizeof (tracer_path), "%s/%s_f.mp4",
            results_dir, out_prefix);

  scalar m[];
  foreach()
    m[] = cs[] - 0.5;

  output_ppm (omega, file = vort_path,
              box = {{-0.5,-0.5},{7.5,0.5}},
              min = -10, max = 10, linear = true, mask = m);
  output_ppm (f, file = tracer_path,
              box = {{-0.5,-0.5},{7.5,0.5}},
              linear = false, min = 0, max = 1, mask = m);
}

event adapt (i++)
{
  adapt_wavelet ({cs,u,f}, (double[]){1e-2,3e-2,3e-2,3e-2},
                 maxlevel, 4);
}

event finalize (t = t_end)
{
  if (pid() == 0) {
    FILE * pvd = fopen (pvd_path, "a");
    if (pvd) {
      fputs ("  </Collection>\n</VTKFile>\n", pvd);
      fclose (pvd);
    }
  }

  return 1;
}
