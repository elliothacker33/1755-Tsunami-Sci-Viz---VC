extern scalar * tracers;
extern face vector uf;
extern double dt;

#if TREE
event defaults (i = 0) {
  for (scalar s in tracers) {
#if EMBED
    s.refine = s.prolongation = refine_embed_linear;
#else
    s.refine  = refine_linear;
#endif
    s.restriction = restriction_volume_average;
    s.dirty = true;
  }
}
#endif

#include "bcg.h"

event tracer_advection (i++,last) {
  advection (tracers, uf, dt);
}

event vof (i++, last);

event tracer_diffusion (i++,last);
