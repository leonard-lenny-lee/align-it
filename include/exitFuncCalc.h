#pragma once

// RDKit
#include <GraphMol/Atom.h>
#include <GraphMol/ROMol.h>

// Align-it
#include <pharmacophore.h>

void exitFuncCalc(RDKit::ROMol *, Pharmacophore *);

unsigned int countExit(const Pharmacophore &p);

Coordinate _exitCalcNormal(RDKit::Atom *a, const RDKit::Conformer &conf);
