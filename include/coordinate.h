/*******************************************************************************
coordinate.h - Align-it

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

#ifndef __SILICOSIT_ALIGNIT_COORDINATE_H__
#define __SILICOSIT_ALIGNIT_COORDINATE_H__

// General
#include <iostream>

class Coordinate {
  public:
    double x;
    double y;
    double z;

    Coordinate(void);
    Coordinate(double, double, double);
};

std::ostream &operator<<(std::ostream &, const Coordinate &);

#endif //__SILICOSIT_ALIGNIT_COORDINATE_H__
