/*******************************************************************************
hAccFuncCalc.cpp - Align-it

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

#include <hAccFuncCalc.h>

#include <GraphMol/PeriodicTable.h>

void hAccFuncCalc(RDKit::ROMol *mol, Pharmacophore *pharmacophore) {
    // Create for every hydrogen acceptor a pharmacophore point
    const auto &conf = mol->getConformer();
    for (auto atom : mol->atoms()) {
        if (atom->getAtomicNum() == 7 || atom->getAtomicNum() == 8) {
            if (atom->getFormalCharge() <= 0) {
                if (_hAccDelocalized(atom)) {
                    continue;
                }
                PharmacophorePoint p;
                p.func = HACC;
                p.point.x = conf.getAtomPos(atom->getIdx()).x;
                p.point.y = conf.getAtomPos(atom->getIdx()).y;
                p.point.z = conf.getAtomPos(atom->getIdx()).z;
                p.hasNormal = true;
                p.alpha = funcSigma[HACC];
                p.normal = _hAccCalcNormal(atom, conf);
                pharmacophore->push_back(p);
            }
        }
    }
}

double _hAccCalcAccSurf(RDKit::Atom *atom, const RDKit::Conformer &conf) {
    double radius(H_BOND_DIST);

    //---(1)-- create sphere with uniformly distributed points
    std::vector<Coordinate> sphere;
    std::vector<Coordinate>::iterator itS;

    const double arclength(1.0 / sqrt(sqrt(3.0) * DENSITY));
    double dphi(arclength / radius);
    int nlayer(ROUND(PI / dphi) + 1);

    double phi(0.0);
    for (int i(0); i < nlayer; ++i) {
        double rsinphi(radius * sin(phi));
        double z(radius * cos(phi));
        double dtheta((rsinphi == 0) ? PI * 2 : arclength / rsinphi);
        int tmpNbrPoints(ROUND(PI * 2 / dtheta));
        if (tmpNbrPoints <= 0) {
            tmpNbrPoints = 1;
        }
        dtheta = PI * 2.0 / tmpNbrPoints;
        double theta((i % 2) ? 0 : PI);
        for (int j(0); j < tmpNbrPoints; ++j) {
            Coordinate coord;
            coord.x = rsinphi * cos(theta) + conf.getAtomPos(atom->getIdx()).x;
            coord.y = rsinphi * sin(theta) + conf.getAtomPos(atom->getIdx()).y;
            coord.z = z + conf.getAtomPos(atom->getIdx()).z;
            sphere.push_back(coord);
            theta += dtheta;
            if (theta > PI * 2) {
                theta -= PI * 2;
            }
        }
        phi += dphi;
    }
    //---(2)-- define neighbors of atom
    std::list<RDKit::Atom *> aList(_hAccGetNeighbors(atom, conf));
    std::list<RDKit::Atom *>::iterator itA;

    //---(3) -- check for every sphere-point if it is accessible
    int nbrAccSurfPoints(0);
    double r;
    for (itS = sphere.begin(); itS != sphere.end(); ++itS) {
        bool isAccessible(true);
        for (itA = aList.begin(); itA != aList.end(); ++itA) {
            RDKit::Atom *n(*itA);
            const auto &p = conf.getAtomPos(n->getIdx());
            double distSq(((itS->x - p.x) * (itS->x - p.x)) +
                          ((itS->y - p.y) * (itS->y - p.y)) +
                          ((itS->z - p.z) * (itS->z - p.z)));
            r = RDKit::PeriodicTable::getTable()->getRvdw(n->getAtomicNum());
            double sumSq((r + H_RADIUS) * (r + H_RADIUS));
            if (distSq <= sumSq) {
                isAccessible = false;
                break;
            }
        }
        if (isAccessible) {
            ++nbrAccSurfPoints;
        }
    }
    return (nbrAccSurfPoints / (double)sphere.size());
}

std::list<RDKit::Atom *> _hAccGetNeighbors(RDKit::Atom *a,
                                           const RDKit::Conformer &conf) {
    std::list<RDKit::Atom *> aList;
    double r;
    for (auto atom : a->getOwningMol().atoms()) {
        if (atom == a)
            continue;
        const auto &p = conf.getAtomPos(a->getIdx());
        const auto &pp = conf.getAtomPos(atom->getIdx());
        r = RDKit::PeriodicTable::getTable()->getRvdw(a->getAtomicNum());
        double delta(H_BOND_DIST + H_RADIUS + r);
        double maxDistSq(delta * delta);
        double distSq((p.x - pp.x) * (p.x - pp.x) +
                      (p.y - pp.y) * (p.y - pp.y) +
                      (p.z - pp.z) * (p.z - pp.z));
        if (distSq <= maxDistSq) {
            aList.push_back(atom);
        }
    }
    return aList;
}

bool _hAccDelocalized(RDKit::Atom *a) {
    if (a->getAtomicNum() != 7) {
        return false;
    }
    if (a->getIsAromatic() && a->getTotalDegree() == 3) {
        return true;
    }
    for (const auto &nbri :
         boost::make_iterator_range(a->getOwningMol().getAtomNeighbors(a))) {
        const auto aa = a->getOwningMol()[nbri];
        if (aa->getIsAromatic() && a->getTotalDegree() == 3) {
            return true;
        }
        if (aa->getAtomicNum() == 6) {
            for (const auto &nbrj : boost::make_iterator_range(
                     aa->getOwningMol().getAtomBonds(aa))) {
                const auto bnd = aa->getOwningMol()[nbrj];
                const auto aaa = bnd->getOtherAtom(aa);
                if (aaa == a)
                    continue;
                if (bnd->getBondTypeAsDouble() == 2.0) {
                    if (aaa->getAtomicNum() == 8)
                        return true;
                    if (aaa->getAtomicNum() == 7)
                        return true;
                    if (aaa->getAtomicNum() == 16)
                        return true;
                }
            }
        } else if (aa->getAtomicNum() == 16) {
            for (const auto &nbrj : boost::make_iterator_range(
                     aa->getOwningMol().getAtomBonds(aa))) {
                const auto bnd = aa->getOwningMol()[nbrj];
                const auto aaa = bnd->getOtherAtom(aa);
                if (aaa == a)
                    continue;
                if ((bnd->getBondTypeAsDouble() == 2.0) &&
                    (aaa->getAtomicNum() == 8)) {
                    return true;
                }
            }
        }
    }
    return false;
}

Coordinate _hAccCalcNormal(RDKit::Atom *a, const RDKit::Conformer &conf) {
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
