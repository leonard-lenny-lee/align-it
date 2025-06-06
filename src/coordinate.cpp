/*******************************************************************************
coordinate.cpp - Align-it

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

#include <coordinate.h>

Coordinate::Coordinate() : x(0.0), y(0.0), z(0.0) {}

Coordinate::Coordinate(double x, double y, double z) : x(x), y(y), z(z) {}

std::ostream &operator<<(std::ostream &os, const Coordinate &A) {
    os << "(" << A.x << "," << A.y << "," << A.z << ")";
    return os;
};
