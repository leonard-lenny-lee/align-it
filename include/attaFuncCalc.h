#pragma once

// RDKit
#include <GraphMol/Atom.h>
#include <GraphMol/ROMol.h>

// Align-it
#include <pharmacophore.h>

void attaFuncCalc(RDKit::ROMol *, Pharmacophore *);

unsigned int countAttaFunc(const Pharmacophore &p);

bool _hasRGroup(RDKit::Atom *a, RDKit::ROMol *m);

Coordinate _attaCalcNormal(RDKit::Atom *a, const RDKit::Conformer &conf);
