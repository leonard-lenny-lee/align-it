/*******************************************************************************
chargeFuncCalc.cpp - Align-it

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

#include <chargeFuncCalc.h>

void chargeFuncCalc(Molecule *m, Pharmacophore *pharmacophore) {
    // Create for every non-zero formal charge a pharmacophore point
    int charge;
    const auto &conf = m->getConformer();
    for (auto atom : m->atoms()) {
        charge = atom->getFormalCharge();
        if (charge < 0) {
            PharmacophorePoint p;
            const auto &point = conf.getAtomPos(atom->getIdx());
            p.func = NEGC;
            p.point.x = point.x;
            p.point.y = point.y;
            p.point.z = point.z;
            p.alpha = funcSigma[NEGC];
            p.hasNormal = false;
            pharmacophore->push_back(p);
        } else if (charge > 0) {
            PharmacophorePoint p;
            const auto &point = conf.getAtomPos(atom->getIdx());
            p.func = POSC;
            p.point.x = point.x;
            p.point.y = point.y;
            p.point.z = point.z;
            p.alpha = funcSigma[POSC];
            p.hasNormal = false;
            pharmacophore->push_back(p);
        }
    }
}
