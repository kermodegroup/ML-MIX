#include "lammpsplugin.h"

#include "version.h"

#include <cstring>

#include "pair_hybrid_overlay_mlml.h"

using namespace LAMMPS_NS;

static Pair *mlmlcreator(LAMMPS *lmp)
{
  return new PairHybridOverlayMLML(lmp);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  plugin.version = LAMMPS_VERSION;
  plugin.style = "pair";
  plugin.name = "hybrid/overlay/mlml";
  plugin.info = "MLML hybrid overlay pair style v0.1";
  plugin.author = "Fraser Birks (fraser.birks@warwick.ac.uk)";
  plugin.creator.v2 = (lammpsplugin_factory2 *) &mlmlcreator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
