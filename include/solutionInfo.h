/*******************************************************************************
solutionInfo.h - Align-it

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

#ifndef __SILICOSIT_ALIGNIT_SOLUTIONINFO_H__
#define __SILICOSIT_ALIGNIT_SOLUTIONINFO_H__

// Align-it
#include <coordinate.h>
#include <siMath.h>

class SolutionInfo {
  public:
    SiMath::Vector rotor;
    double volume;
    unsigned int iterations;
    Coordinate center1;
    Coordinate center2;
    SiMath::Matrix rotation1;
    SiMath::Matrix rotation2;

    SolutionInfo(void);
    ~SolutionInfo(void);
};

#endif //__SILICOSIT_ALIGNIT_SOLUTIONINFO_H__
