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

#ifndef SCIMATH_ITERATIVERANGESTATISTICS_H
#define SCIMATH_ITERATIVERANGESTATISTICS_H

#include <casacore/casa/aips.h>

#include <casacore/scimath/Mathematics/ConstrainedRangeStatistics.h>

namespace casacore {

// Abstract class for allowing easy addition of derived classes to compute
// statistics based on algorithms that operate by iteratively changing
// the range of acceptable points in order to iteratively exclude outliers.

template <class AccumType, class DataIterator, class MaskIterator=const Bool*, class WeightsIterator=DataIterator>
class IterativeRangeStatistics
    : public ConstrainedRangeStatistics<CASA_STATP> {
public:

    virtual ~IterativeRangeStatistics();

    // copy semantics
    IterativeRangeStatistics<CASA_STATP>& operator=(
        const IterativeRangeStatistics<CASA_STATP>& other
    );

    // reset object to initial state. Clears all private fields including data,
    // accumulators, global range. It does not affect the fence factor (_f), which was
    // set at object construction.
    virtual void reset();

    // This class does not allow statistics to be calculated as datasets are added, so
    // an exception will be thrown if <src>c</src> is True.
    void setCalculateAsAdded(Bool c);

    // get the number of iterations
    uInt getNiter() const { return _niter; }

protected:

    // if maxIterations is positive, only iterate for a maximum of that
    // many iterations. Else, iterate up to a max of a very large number
    // of iterations (see code for current setting).
    IterativeRangeStatistics(Int maxIterations);

    // Derived classes must implement. Describes how to compute new good range
    // after each iteration. sd is the set of stats computed in the current iteration.
    virtual CountedPtr<std::pair<AccumType, AccumType> > _setNewRange(const StatsData<AccumType>& sd) = 0;

    // does the iterative range setting. Calls _setNewRange().
    void _setRange();

private:

    Int _maxIterations;
    Bool _rangeIsSet;
    uInt _niter;

};

}

#ifndef CASACORE_NO_AUTO_TEMPLATES
#include <casacore/scimath/Mathematics/IterativeRangeStatistics.tcc>
#endif //# CASACORE_NO_AUTO_TEMPLATES

#endif
