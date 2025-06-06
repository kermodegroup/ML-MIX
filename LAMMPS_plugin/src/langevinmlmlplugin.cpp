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

#include "fix_langevin_mlml.h"

using namespace LAMMPS_NS;

static Fix *mlmlcreator(LAMMPS *lmp, int argc, char **argv)
{
  return new FixLangevinMLML(lmp, argc, argv);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  plugin.version = LAMMPS_VERSION;
  plugin.style = "fix";
  plugin.name = "langevin/mlml";
  plugin.info = "Langevin MLML fix style v0.3.0";
  plugin.author = "Fraser Birks (fraser.birks@warwick.ac.uk)";
  plugin.creator.v2 = (lammpsplugin_factory2 *) &mlmlcreator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
