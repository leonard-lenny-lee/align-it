#pragma once

// RDKit
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>

// Align-it
#include <pharmacophore.h>

void attaFuncCalc(RDKit::ROMol *, Pharmacophore *);
bool _hasRGroup(RDKit::Atom *a, RDKit::ROMol *m);
Coordinate _attaCalcNormal(RDKit::Atom *a, const RDKit::Conformer &conf);
