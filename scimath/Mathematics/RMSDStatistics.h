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

#ifndef SCIMATH_RMSDSTATISTICS_H
#define SCIMATH_RMSDSTATISTICS_H

#include <casacore/casa/aips.h>

#include <casacore/scimath/Mathematics/ConstrainedRangeStatistics.h>

#include <set>
#include <vector>
#include <utility>

namespace casacore {

// Class to calculate statistics using the "RMSD" algorithm. This method
// iteratively calculates statistics by discarding outliers that have
// values outside the range (mean +/- f*rms). f is a positive constant
// specified by the user.

template <class AccumType, class DataIterator, class MaskIterator=const Bool*, class WeightsIterator=DataIterator>
class RMSDStatistics
    : public ConstrainedRangeStatistics<CASA_STATP> {

public:

    // f is multiplied by the rms value of the remaining data set. Points with
    // values outside the range (mean +/- f*rms) will be discarded in the next iteration.
    // f must be positive or an exception is thrown. maxIterations is the maximum number
    // of iterations to carry out before exiting. A non-positive value means to perform
    // a maximum of some very large, unspecified number of iterations (see code for current
    // setting).
    RMSDStatistics(Double f=3, Int maxIterations=0);

    virtual ~RMSDStatistics();

    // copy semantics
    RMSDStatistics<CASA_STATP>& operator=(
        const RMSDStatistics<CASA_STATP>& other
    );

    // get the algorithm that this object uses for computing stats
    virtual StatisticsData::ALGORITHM algorithm() const {
        return StatisticsData::RMSD;
    };

    // reset object to initial state. Clears all private fields including data,
    // accumulators, global range. It does not affect the fence factor (_f), which was
    // set at object construction.
    virtual void reset();

    // This class does not allow statistics to be calculated as datasets are added, so
    // an exception will be thrown if <src>c</src> is True.
    void setCalculateAsAdded(Bool c);

    // get the number of iterations
    uInt getNiter() const { return _niter; }

private:

    Double _f;
    Int _maxIterations;
    Bool _rangeIsSet;
    uInt _niter;

    void _setRange();
};

}

#ifndef CASACORE_NO_AUTO_TEMPLATES
#include <casacore/scimath/Mathematics/RMSDStatistics.tcc>
#endif //# CASACORE_NO_AUTO_TEMPLATES

#endif
