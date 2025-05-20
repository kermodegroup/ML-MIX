/* ----------------------------------------------------------------------
   This file is part of a plugin for LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2025) Fraser Birks
   This file was written by Fraser Birks and is distributed under the 
   GNU General Public License.

   This file is not part of the original LAMMPS distribution but is designed 
   to be used with LAMMPS. It is provided under the same GPLv2 license as LAMMPS 
   to ensure compatibility.

   See the LICENSE file for details.
------------------------------------------------------------------------- */

#include "lammpsplugin.h"

#include "version.h"

#include <cstring>

#include "fix_mlml.h"
#include "fix_mlml_kokkos.h"

using namespace LAMMPS_NS;

static Fix *mlmlcreator(LAMMPS *lmp, int argc, char **argv)
{
  return new FixMLML(lmp, argc, argv);
}

static Fix *mlmlkkcreator(LAMMPS *lmp, int argc, char **argv)
{
  return new FixMLMLKokkos<LMPDeviceType>(lmp, argc, argv);
}

static Fix *mlmlkkhostcreator(LAMMPS *lmp, int argc, char **argv)
{
  return new FixMLMLKokkos<LMPHostType>(lmp, argc, argv);
}

static Fix *mlmlkkdevicecreator(LAMMPS *lmp, int argc, char **argv)
{
  return new FixMLMLKokkos<LMPDeviceType>(lmp, argc, argv);
}


extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  // register base plugin
  plugin.version = LAMMPS_VERSION;
  plugin.style = "fix";
  plugin.name = "mlml";
  plugin.info = "MLML fix style v0.2.0";
  plugin.author = "Fraser Birks (fraser.birks@warwick.ac.uk)";
  plugin.creator.v2 = (lammpsplugin_factory2 *) &mlmlcreator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);

  // register learn/kk pair style
  plugin.name = "mlml/kk";
  plugin.info = "mlml/kk fix style";
  plugin.author = "Fraser Birks (fraser.birks@warwick.ac.uk)";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &mlmlkkcreator;
  (*register_plugin)(&plugin, lmp);

  // also register learn/kk/host pair style. only need to update changed fields
  plugin.name = "mlml/kk/host";
  plugin.info = "mlml/kk fix style host";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &mlmlkkhostcreator;
  (*register_plugin)(&plugin, lmp);

  // also register learn/kk/device pair style. only need to update changed fields
  plugin.name = "mlml/kk/device";
  plugin.info = "mlml/kk fix style device";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &mlmlkkdevicecreator;
  (*register_plugin)(&plugin, lmp);
  
}
