//# Copyright (C) 2000,2001
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id: Array.h 21545 2015-01-22 19:36:35Z gervandiepen $

#ifndef SCIMATH_CHAUVENETCRITERIONSTATISTICS_TCC
#define SCIMATH_CHAUVENETCRITERIONSTATISTICS_TCC

#include <casacore/scimath/Mathematics/ChauvenetCriterionStatistics.h>

#include <casacore/scimath/Mathematics/ZScoreCalculator.h>

namespace casacore {

CASA_STATD
ChauvenetCriterionStatistics<CASA_STATP>::ChauvenetCriterionStatistics(
    Double zscore, Int maxIterations
)
  : IterativeRangeStatistics<CASA_STATP>(maxIterations),
    _zscore(zscore) {}

CASA_STATD
ChauvenetCriterionStatistics<CASA_STATP>::~ChauvenetCriterionStatistics() {}

CASA_STATD
ChauvenetCriterionStatistics<CASA_STATP>&
ChauvenetCriterionStatistics<CASA_STATP>::operator=(
    const ChauvenetCriterionStatistics<CASA_STATP>& other
) {
    if (this == &other) {
        return *this;
    }
    IterativeRangeStatistics<CASA_STATP>::operator=(other);
    _zscore = other._zscore;
    return *this;
}

CASA_STATD
CountedPtr<std::pair<AccumType, AccumType> > ChauvenetCriterionStatistics<CASA_STATP>::_setNewRange(
    const StatsData<AccumType>& sd
) {
    Double zScore = _zscore >= 0 ? _zscore : ZScoreCalculator::getMaxZScore((uInt64)sd.npts);
    return new std::pair<AccumType, AccumType>(
        sd.mean - zScore*sd.stddev, sd.mean + zScore*sd.stddev
    );
}

}

#endif
