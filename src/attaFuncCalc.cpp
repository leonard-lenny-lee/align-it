#include <attaFuncCalc.h>

void attaFuncCalc(RDKit::ROMol *m, Pharmacophore *pharmacophore) {
    // Create a point for every atom with a R group connection
    const auto &conf = m->getConformer();
    for (auto a : m->atoms()) {
        if (a->getAtomicNum() <= 1)
            continue;
        if (!_hasRGroup(a, m)) {
            continue;
        }
        PharmacophorePoint p;
        const auto &point = conf.getAtomPos(a->getIdx());
        p.func = FuncGroup::ATTA;
        p.point.x = point.x;
        p.point.y = point.y;
        p.point.z = point.z;
        p.alpha = funcSigma[FuncGroup::ATTA];
        p.hasNormal = true;
        p.normal = _attaCalcNormal(a, conf);
        pharmacophore->push_back(p);
    }
}

unsigned int countAttaFunc(const Pharmacophore &p) {
    int count(0);
    for (const auto &point : p) {
        if (point.func == FuncGroup::ATTA) {
            count++;
        }
    }
    return count;
}

bool _hasRGroup(RDKit::Atom *a, RDKit::ROMol *m) {
    for (const auto &nbri :
         boost::make_iterator_range(m->getAtomNeighbors(a))) {
        const auto &aa = (*m)[nbri];
        if (aa->getAtomicNum() == 0) {
            return true;
        }
    }
    return false;
}

Coordinate _attaCalcNormal(RDKit::Atom *a, const RDKit::Conformer &conf) {
    Coordinate normal;
    const auto &p = conf.getAtomPos(a->getIdx());
    for (const auto &nbri :
         boost::make_iterator_range(a->getOwningMol().getAtomNeighbors(a))) {
        const auto aa = a->getOwningMol()[nbri];
        if (aa->getAtomicNum() == 0) {
            const auto &pp = conf.getAtomPos(aa->getIdx());
            normal.x += (p.x - pp.x);
            normal.y += (p.y - pp.y);
            normal.z += (p.z - pp.z);
            break;
        }
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
