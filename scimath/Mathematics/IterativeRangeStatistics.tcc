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

#ifndef SCIMATH_ITERATIVERANGESTATISTICS_TCC
#define SCIMATH_ITERATIVERANGESTATISTICS_TCC

#include <casacore/scimath/Mathematics/IterativeRangeStatistics.h>

namespace casacore {

CASA_STATD
IterativeRangeStatistics<CASA_STATP>::IterativeRangeStatistics(
    Int maxIterations
)
  : ConstrainedRangeStatistics<CASA_STATP>(),
    _maxIterations(maxIterations), _rangeIsSet(False), _niter(0) {}

CASA_STATD
IterativeRangeStatistics<CASA_STATP>::~IterativeRangeStatistics() {}

CASA_STATD
IterativeRangeStatistics<CASA_STATP>&
IterativeRangeStatistics<CASA_STATP>::operator=(
    const IterativeRangeStatistics<CASA_STATP>& other
) {
    if (this == &other) {
        return *this;
    }
    ClassicalStatistics<CASA_STATP>::operator=(other);
    _maxIterations = other._maxIterations;
    _niter = other._niter;
    _rangeIsSet = other._rangeIsSet;
    return *this;
}

CASA_STATD
void IterativeRangeStatistics<CASA_STATP>::reset() {
    ConstrainedRangeStatistics<CASA_STATP>::reset();
    _rangeIsSet = False;
}

CASA_STATD
void IterativeRangeStatistics<CASA_STATP>::setCalculateAsAdded(
    Bool c
) {
    ThrowIf(
        c, "This algorithm does not support calculating statistics "
        "incrementally as data sets are added"
    );
}

CASA_STATD
void IterativeRangeStatistics<CASA_STATP>::_setRange() {
    if (_rangeIsSet) {
        return;
    }
    uInt maxI = _maxIterations >= 0 ? _maxIterations : 1000;
    uInt prevNpts = 0;
    StatsData<AccumType> sd;
    while (_niter <= maxI) {
        if (_niter == 0) {
            ClassicalStatistics<CASA_STATP> cs(*this);
            sd = cs.getStatistics();
        }
        else {
            sd = this->getStatistics();
            if ((uInt64)sd.npts == prevNpts) {
                break;
            }
        }
        ConstrainedRangeStatistics<CASA_STATP>::_setRange(_setNewRange(sd));
        // _rangeIsSet is set here to prevent infinite recursion on next loop iteration
        _rangeIsSet = True;
        prevNpts = (uInt64)sd.npts;
        ++_niter;
    }
    --_niter;
}

}

#endif
