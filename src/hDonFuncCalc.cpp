/*******************************************************************************
hDonFuncCalc.cpp - Align-it

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

#include <hDonFuncCalc.h>

void hDonFuncCalc(RDKit::ROMol *mol, Pharmacophore *pharmacophore) {
    // Create for every hydrogen donor a pharmacophore point
    const auto &conf = mol->getConformer();
    for (auto atom : mol->atoms()) {
        if (atom->getAtomicNum() == 7 || atom->getAtomicNum() == 8) {
            if (atom->getFormalCharge() >= 0 &&
                (atom->getTotalNumHs(true) != 0)) {
                PharmacophorePoint p;
                p.func = HDON;
                p.point.x = conf.getAtomPos(atom->getIdx()).x;
                p.point.y = conf.getAtomPos(atom->getIdx()).y;
                p.point.z = conf.getAtomPos(atom->getIdx()).z;
                p.hasNormal = true;
                p.alpha = funcSigma[HDON];
                p.normal = _hDonCalcNormal(atom, conf);
                pharmacophore->push_back(p);
            }
        }
    }
}

Coordinate _hDonCalcNormal(RDKit::Atom *a, const RDKit::Conformer &conf) {
    Coordinate normal;
    const auto &p = conf.getAtomPos(a->getIdx());
    for (const auto &nbri :
         boost::make_iterator_range(a->getOwningMol().getAtomNeighbors(a))) {
        const auto aa = a->getOwningMol()[nbri];
        if (aa->getAtomicNum() == 1)
            continue;
        const auto &pp = conf.getAtomPos(aa->getIdx());
        normal.x += (pp.x - p.x);
        normal.y += (pp.y - p.y);
        normal.z += (pp.z - p.z);
    }
    double length(
        sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z));
    normal.x /= length;
    normal.y /= length;
    normal.z /= length;
    normal.x = -normal.x;
    normal.y = -normal.y;
    normal.z = -normal.z;
    normal.x += p.x;
    normal.y += p.y;
    normal.z += p.z;
    return normal;
}
