/*******************************************************************************
alignLib.h - Align-it

Copyright 2021 by OliverBScott and the Align-it contributors

This file is part of Align-it.

        Align-it is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published
        by the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        Align-it is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with Align-it.  If not, see <http://www.gnu.org/licenses/>.

***********************************************************************/

#ifndef __SILICOSIT_ALIGNIT_ALIGNLIB_H__
#define __SILICOSIT_ALIGNIT_ALIGNLIB_H__

// General
#include <tuple>

// Align-it
#include <alignment.h>
#include <aromFuncCalc.h>
#include <exitFuncCalc.h>
#include <chargeFuncCalc.h>
#include <functionMapping.h>
#include <hAccFuncCalc.h>
#include <hDonFuncCalc.h>
#include <hybridCalc.h>
#include <lipoFuncCalc.h>
#include <pharMerger.h>
#include <pharmacophore.h>
#include <result.h>
#include <siMath.h>
#include <solutionInfo.h>

// Toolkit
#include <GraphMol/ROMol.h>
using Molecule = RDKit::ROMol;

namespace alignit {

Pharmacophore calcPharmacophore(Molecule &mol, bool calcArom = true,
                                bool calcHDon = true, bool calcHAcc = true,
                                bool calcLipo = true, bool calcCharge = true,
                                bool calcHybrid = true, bool calcExits = true);

void mergePharmacophore(Pharmacophore &p);

Result alignPharmacophores(Pharmacophore &ref, Pharmacophore &db,
                           double epsilon, bool useNormals, bool useExclusion,
                           Molecule *dbMol);

std::tuple<Pharmacophore, Result>
alignMols(Molecule &refMol, Molecule &refDb, bool calcArom = true,
          bool calcHDon = true, bool calcHAcc = true, bool calcLipo = true,
          bool calcCharge = true, bool calcHybrid = true, bool merge = false,
          double epsilon = 0.5, bool useNormals = true,
          bool useExclusion = false, bool calcExits = true);

} // namespace alignit

#endif //__SILICOSIT_ALIGNIT_ALIGNLIB_H__
