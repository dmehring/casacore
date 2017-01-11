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

#ifndef SCIMATH_CHAUVENETCRITERIONSTATISTICS_H
#define SCIMATH_CHAUVENETCRITERIONSTATISTICS_H

#include <casacore/casa/aips.h>

#include <casacore/scimath/Mathematics/IterativeRangeStatistics.h>

namespace casacore {

// Class to calculate statistics using the so-called Chauvenet criterion. This method
// iteratively calculates statistics by discarding outliers on the basis of Chauvenet's
// criterion, until the specified maximum number of iterations is reached, or the final
// iteration results in no additional points being discarded.
// Alternatively, one can specify a z score which indicates the number of standard deviations
// beyond which to discard points. In this case, no iterating is done.

template <class AccumType, class DataIterator, class MaskIterator=const Bool*, class WeightsIterator=DataIterator>
class ChauvenetCriterionStatistics
    : public IterativeRangeStatistics<CASA_STATP> {
public:

    // If <src>zscore</src> is not negative, use that value to discard outliers beyond
    // zscore standard deviations from the mean, and compute statistics based on the
    // remaining data. If <src>zscore</src> is negative, use Chauvenet's Criterion to
    // determine which outliers to discard. <src>maxIterations</src> is the maximum
    // number of iterations to use before stopping. If negative, continue iterating until the
    // set zscore or Chauvenet's criterion is met (ie that there are no remaining outliers).
    ChauvenetCriterionStatistics(Double zscore=-1, Int maxIterations=0);

    virtual ~ChauvenetCriterionStatistics();

    // copy semantics
    ChauvenetCriterionStatistics<CASA_STATP>& operator=(
        const ChauvenetCriterionStatistics<CASA_STATP>& other
    );

    // get the algorithm that this object uses for computing stats
    virtual StatisticsData::ALGORITHM algorithm() const {
        return StatisticsData::CHAUVENETCRITERION;
    };

protected:

    CountedPtr<std::pair<AccumType, AccumType> > _setNewRange(const StatsData<AccumType>& sd);

private:

    Double _zscore;

};

}

#ifndef CASACORE_NO_AUTO_TEMPLATES
#include <casacore/scimath/Mathematics/ChauvenetCriterionStatistics.tcc>
#endif //# CASACORE_NO_AUTO_TEMPLATES

#endif
