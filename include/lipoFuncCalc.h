/*******************************************************************************
lipoFuncCalc.h - Align-it

Copyright 2012-2013 by Silicos-it, a division of Imacosi BVBA

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

#ifndef __SILICOSIT_ALIGNIT_LIPOFUNCCALC_H__
#define __SILICOSIT_ALIGNIT_LIPOFUNCCALC_H__

// General
#include <algorithm>
#include <list>
#include <set>

// Toolkit
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RingInfo.h>
using Molecule = RDKit::ROMol;

// Align-it
#include <defines.h>
#include <pharmacophore.h>

void lipoFuncCalc(Molecule *, Pharmacophore *);
void _lipoLabelAtoms(Molecule *);
void _lipoGroupAtoms(Molecule *, Pharmacophore *);

double _lipoCalcAccSurf(RDKit::Atom *, const RDKit::Conformer &);
std::list<RDKit::Atom *> _lipoGetNeighbors(RDKit::Atom *,
                                           const RDKit::Conformer &);
void _lipoLabelNeighbors(RDKit::Atom *, double);

#endif //__SILICOSIT_ALIGNIT_LIPOFUNCCALC_H__
