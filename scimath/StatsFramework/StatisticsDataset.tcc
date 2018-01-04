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

#ifndef SCIMATH_STATISTICSDATASET_TCC
#define SCIMATH_STATISTICSDATASET_TCC

#include <casacore/scimath/StatsFramework/StatisticsDataset.h>

#include <casacore/casa/Utilities/PtrHolder.h>
#include <casacore/scimath/StatsFramework/ClassicalStatisticsData.h>

namespace casacore {

CASA_STATD StatisticsDataset<CASA_STATP>::StatisticsDataset()
    : _data(), _weights(), _masks(), _counts(), _dataStrides(), _maskStrides(),
      _isIncludeRanges(), _dataRanges(), _dataProvider(NULL), _idataset(0)
{}

CASA_STATD
StatisticsDataset<CASA_STATP>::StatisticsDataset(const StatisticsDataset& other)
    : _data(other._data), _weights(other._weights), _masks(other._masks),
      _counts(other._counts), _dataStrides(other._dataStrides),
      _maskStrides(other._maskStrides),
      _isIncludeRanges(other._isIncludeRanges), _dataRanges(other._dataRanges),
      // WARN reference semantics
      _dataProvider(other._dataProvider), _idataset(0) {}

CASA_STATD StatisticsDataset<CASA_STATP>::~StatisticsDataset() {}

CASA_STATD StatisticsDataset<CASA_STATP>&
StatisticsDataset<CASA_STATP>::operator=(
    const StatisticsDataset<CASA_STATP>& other
) {
     if (this == &other) {
         return *this;
     }
     _data = other._data;
     _weights = other._weights;
     _masks = other._masks;
     _counts = other._counts;
     _dataStrides = other._dataStrides;
     _maskStrides = other._maskStrides;
     _isIncludeRanges = other._isIncludeRanges;
     _dataRanges = other._dataRanges;
     // WARN reference semantics
     _dataProvider = other._dataProvider;
     _idataset = _idataset;
     return *this;
}

CASA_STATD void StatisticsDataset<CASA_STATP>::addData(
    const DataIterator& first, uInt nr, uInt dataStride, Bool nrAccountsForStride
) {
    _throwIfDataProviderDefined();
    _data.push_back(first);
    // internally we store the number of strided points
    _counts.push_back(
        nrAccountsForStride ? nr
            : nr % dataStride == 0
              ? nr/dataStride
                : nr/dataStride + 1
    );
    _dataStrides.push_back(dataStride);
}

CASA_STATD void StatisticsDataset<CASA_STATP>::addData(
    const DataIterator& first, uInt nr,
    const DataRanges& dataRanges, Bool isInclude, uInt dataStride,
    Bool nrAccountsForStride
) {
    _throwIfDataProviderDefined();
    typename DataRanges::const_iterator riter = dataRanges.begin();
    typename DataRanges::const_iterator rend = dataRanges.end();
    while (riter != rend) {
        ThrowIf(
            (*riter).first > (*riter).second,
            "The first value in a range pair cannot be greater than the second"
        );
        ++riter;
    }
    uInt n = _data.size();
    _isIncludeRanges[n] = isInclude;
    _dataRanges[n] = dataRanges;
    addData(first, nr, dataStride, nrAccountsForStride);
}

CASA_STATD void StatisticsDataset<CASA_STATP>::addData(
    const DataIterator& first, const MaskIterator& maskFirst,
    uInt nr, uInt dataStride, Bool nrAccountsForStride, uInt maskStride
) {
    _throwIfDataProviderDefined();
    uInt key = _data.size();
    _maskStrides[key] = maskStride;
    _masks[key] = maskFirst;
    addData(first, nr, dataStride, nrAccountsForStride);
}

CASA_STATD void StatisticsDataset<CASA_STATP>::addData(
    const DataIterator& first, const MaskIterator& maskFirst,
    uInt nr, const DataRanges& dataRanges,
    Bool isInclude, uInt dataStride, Bool nrAccountsForStride,
    uInt maskStride
) {
    _throwIfDataProviderDefined();
    uInt key = _data.size();
    _maskStrides[key] = maskStride;
    _masks[key] = maskFirst;
    addData(
        first, nr, dataRanges, isInclude,
        dataStride, nrAccountsForStride
    );
}

CASA_STATD void StatisticsDataset<CASA_STATP>::addData(
    const DataIterator& first, const WeightsIterator& weightFirst,
    uInt nr, uInt dataStride, Bool nrAccountsForStride
) {
    _throwIfDataProviderDefined();
    _weights[_data.size()] = weightFirst;
    addData(first, nr, dataStride, nrAccountsForStride);
}

CASA_STATD void StatisticsDataset<CASA_STATP>::addData(
    const DataIterator& first, const WeightsIterator& weightFirst,
    uInt nr, const DataRanges& dataRanges,
    Bool isInclude, uInt dataStride, Bool nrAccountsForStride
) {
    _throwIfDataProviderDefined();
    _weights[_data.size()] = weightFirst;
    addData(
        first, nr, dataRanges, isInclude, dataStride, nrAccountsForStride
    );
}

CASA_STATD void StatisticsDataset<CASA_STATP>::addData(
    const DataIterator& first, const WeightsIterator& weightFirst,
    const MaskIterator& maskFirst, uInt nr, uInt dataStride,
    Bool nrAccountsForStride, uInt maskStride
) {
    _throwIfDataProviderDefined();
    _weights[_data.size()] = weightFirst;
    addData(
        first, maskFirst, nr, dataStride, nrAccountsForStride, maskStride
    );
}

CASA_STATD void StatisticsDataset<CASA_STATP>::addData(
    const DataIterator& first, const WeightsIterator& weightFirst,
    const MaskIterator& maskFirst, uInt nr, const DataRanges& dataRanges,
    Bool isInclude, uInt dataStride, Bool nrAccountsForStride,
    uInt maskStride
) {
    _throwIfDataProviderDefined();
    _weights[_data.size()] = weightFirst;
    addData(
        first, maskFirst, nr, dataRanges, isInclude, dataStride,
        nrAccountsForStride, maskStride
    );
}

CASA_STATD Bool StatisticsDataset<CASA_STATP>::empty() const {
    return ! _dataProvider && _data.empty();
}

CASA_STATD
Bool StatisticsDataset<CASA_STATP>::increment(Bool includeIDataset) {
    if (includeIDataset) {
        ++_idataset;
    }
    if (_dataProvider) {
        ++(*_dataProvider);
        if (_dataProvider->atEnd()) {
            _dataProvider->finalize();
            return True;
        }
    }
    else {
        ++_diter;
        if (_diter == _dend) {
            return True;
        }
        ++_citer;
        ++_dsiter;
        ++_dataCount;
    }
    return False;
}

CASA_STATD
void StatisticsDataset<CASA_STATP>::incrementThreadIters(
    DataIterator& dataIter, MaskIterator& maskIter,
    WeightsIterator& weightsIter, uInt64& offset, uInt nthreads
) const {
    uInt increment = nthreads*ClassicalStatisticsData::BLOCK_SIZE*_chunkStride;
    if (offset+increment >= _chunkCount*_chunkStride) {
        // necessary because in some cases std::advance will segfault
        // if advanced past the end of the data structure
        return;
    }
    std::advance(dataIter, increment);
    if (_chunkHasWeights) {
        std::advance(weightsIter, increment);
    }
    if (_chunkHasMask) {
        std::advance(maskIter, nthreads*ClassicalStatisticsData::BLOCK_SIZE*_chunkMaskStride);
    }
    offset += increment;
}

CASA_STATD
void StatisticsDataset<CASA_STATP>::initIterators() {
    ThrowIf(empty(), "No data sets have been added");
    if (_dataProvider) {
        _dataProvider->reset();
    }
    else {
        _dataCount = 0;
        // const std::vector<DataIterator>& data = this->_getDataset().getData();
        _diter = _data.begin();
        _dend = _data.end();
        // const std::vector<uInt>& dataStrides = this->_getDataset().getDataStrides();
        _dsiter = _dataStrides.begin();
        // const std::vector<Int64>& counts = this->_getDataset().getCounts();
        _citer = _counts.begin();
        //_masks = this->_getDataset().getMasks();
        //_weights = this->_getDataset().getWeights();
        //_ranges = this->_getDataset().getRanges();
        //_isIncludeRanges = this->_getDataset().getIsIncludeRanges();
    }
    _chunkHasRanges = False;
    _chunkRanges.clear();
    _chunkIsIncludeRanges = False;
    _chunkHasMask = False;
    _chunkHasWeights = False;
}

CASA_STATD
void StatisticsDataset<CASA_STATP>::initLoopVars(uInt64& chunkCount, Bool& chunkHasWeights) {
    if (_dataProvider) {
        _chunkData = _dataProvider->getData();
        _chunkCount = _dataProvider->getCount();
        _chunkStride = _dataProvider->getStride();
        _chunkHasRanges = _dataProvider->hasRanges();
        if (_chunkHasRanges) {
            _chunkRanges = _dataProvider->getRanges();
            _chunkIsIncludeRanges = _dataProvider->isInclude();
        }
        _chunkHasMask = _dataProvider->hasMask();
        if (_chunkHasMask) {
            _chunkMask = _dataProvider->getMask();
            _chunkMaskStride = _dataProvider->getMaskStride();
        }
        _chunkHasWeights = _dataProvider->hasWeights();
        if (_chunkHasWeights) {
            _chunkWeights = _dataProvider->getWeights();
        }
    }
    else {
        _chunkData = *_diter;
        _chunkCount = *_citer;
        _chunkStride = *_dsiter;
        typename std::map<uInt, DataRanges>::const_iterator rangeI = _chunkRanges.find(_dataCount);
        _chunkHasRanges = rangeI != _chunkRanges.end();
        if (_chunkHasRanges) {
            _chunkRanges = rangeI->second;
            _chunkIsIncludeRanges = _isIncludeRanges.find(_dataCount)->second;
        }
        typename std::map<uInt, MaskIterator>::const_iterator maskI = _masks.find(_dataCount);
        _chunkHasMask = maskI != _masks.end();
        if (_chunkHasMask) {
            _chunkMask = maskI->second;
            _chunkMaskStride = _maskStrides.find(_dataCount)->second;
        }
        _chunkHasWeights = _weights.find(_dataCount) != _weights.end();
        if (_chunkHasWeights) {
            _chunkWeights = _weights.find(_dataCount)->second;
        }
    }
    chunkCount = _chunkCount;
    chunkHasWeights = _chunkHasWeights;
}

CASA_STATD
void StatisticsDataset<CASA_STATP>::initThreadVars(
    uInt& nBlocks, uInt64& extra, uInt& nthreads, PtrHolder<DataIterator>& dataIter,
    PtrHolder<MaskIterator>& maskIter, PtrHolder<WeightsIterator>& weightsIter,
    PtrHolder<uInt64>& offset, uInt nThreadsMax
) const {
    uInt n = ClassicalStatisticsData::CACHE_PADDING*nThreadsMax;
    dataIter.set(new DataIterator[n], True);
    maskIter.set(new MaskIterator[n], True);
    weightsIter.set(new WeightsIterator[n], True);
    offset.set(new uInt64[n], True);
    nBlocks = _chunkCount/ClassicalStatisticsData::BLOCK_SIZE;
    extra = _chunkCount % ClassicalStatisticsData::BLOCK_SIZE;
    if (extra > 0) {
        ++nBlocks;
    }
    nthreads = min(nThreadsMax, nBlocks);
    for (uInt tid=0; tid<nthreads; ++tid) {
        // advance the per-thread iterators to their correct starting
        // locations
        uInt idx8 = ClassicalStatisticsData::CACHE_PADDING*tid;
        dataIter[idx8] = _chunkData;
        offset[idx8] = tid*ClassicalStatisticsData::BLOCK_SIZE*_chunkStride;
        std::advance(dataIter[idx8], offset[idx8]);
        if (_chunkHasWeights) {
            weightsIter[idx8] = _chunkWeights;
            std::advance(weightsIter[idx8], offset[idx8]);
        }
        if (_chunkHasMask) {
            maskIter[idx8] = _chunkMask;
            std::advance(maskIter[idx8], tid*ClassicalStatisticsData::BLOCK_SIZE*_chunkMaskStride);
        }
    }
}

CASA_STATD void StatisticsDataset<CASA_STATP>::reset() {
    _data.clear();
    _counts.clear();
    _masks.clear();
    _weights.clear();
    _dataRanges.clear();
    _dataStrides.clear();
    _maskStrides.clear();
    _dataProvider = NULL;
}

CASA_STATD void StatisticsDataset<CASA_STATP>::setData(
    const DataIterator& first, uInt nr, uInt dataStride, Bool nrAccountsForStride
) {
    reset();
    addData(first, nr, dataStride, nrAccountsForStride);
}

CASA_STATD void StatisticsDataset<CASA_STATP>::setData(
    const DataIterator& first, uInt nr,
    const DataRanges& dataRanges, Bool isInclude, uInt dataStride,
    Bool nrAccountsForStride
) {
    reset();
    addData(
        first, nr, dataRanges, isInclude, dataStride, nrAccountsForStride
    );
}

CASA_STATD void StatisticsDataset<CASA_STATP>::setData(
    const DataIterator& first, const MaskIterator& maskFirst,
    uInt nr, uInt dataStride, Bool nrAccountsForStride, uInt maskStride
) {
    reset();
    addData(
        first, maskFirst, nr, dataStride, nrAccountsForStride, maskStride
    );
}

CASA_STATD void StatisticsDataset<CASA_STATP>::setData(
    const DataIterator& first, const MaskIterator& maskFirst,
    uInt nr, const DataRanges& dataRanges,
    Bool isInclude, uInt dataStride, Bool nrAccountsForStride,
    uInt maskStride
) {
    reset();
    addData(
        first, maskFirst, nr, dataRanges, isInclude, dataStride,
        nrAccountsForStride, maskStride
    );
}

CASA_STATD void StatisticsDataset<CASA_STATP>::setData(
    const DataIterator& first, const WeightsIterator& weightFirst,
    uInt nr, uInt dataStride, Bool nrAccountsForStride
) {
    reset();
    addData(first, weightFirst, nr, dataStride, nrAccountsForStride);
}

CASA_STATD void StatisticsDataset<CASA_STATP>::setData(
    const DataIterator& first, const WeightsIterator& weightFirst,
    uInt nr, const DataRanges& dataRanges,
    Bool isInclude, uInt dataStride, Bool nrAccountsForStride
) {
    reset();
    addData(
        first, weightFirst, nr, dataRanges, isInclude,
        dataStride, nrAccountsForStride
    );
}

CASA_STATD void StatisticsDataset<CASA_STATP>::setData(
    const DataIterator& first, const WeightsIterator& weightFirst,
    const MaskIterator& maskFirst, uInt nr, uInt dataStride,
    Bool nrAccountsForStride, uInt maskStride
) {
    reset();
    addData(
        first, weightFirst, maskFirst, nr, dataStride,
        nrAccountsForStride, maskStride
    );
}

CASA_STATD void StatisticsDataset<CASA_STATP>::setData(
    const DataIterator& first, const WeightsIterator& weightFirst,
    const MaskIterator& maskFirst, uInt nr, const DataRanges& dataRanges,
    Bool isInclude, uInt dataStride, Bool nrAccountsForStride, uInt maskStride
) {
    reset();
    addData(
        first, weightFirst, maskFirst, nr, dataRanges, isInclude,
        dataStride, nrAccountsForStride, maskStride
    );
}

CASA_STATD void StatisticsDataset<CASA_STATP>::setDataProvider(
    StatsDataProvider<CASA_STATP> *dataProvider
) {
    ThrowIf(! dataProvider, "Logic Error: data provider cannot be NULL");
    reset();
    _dataProvider = dataProvider;
}

CASA_STATD void StatisticsDataset<CASA_STATP>::_throwIfDataProviderDefined() const {
    ThrowIf(
        _dataProvider,
        "Logic Error: Cannot add data after a data provider has been set. Call setData() to clear "
        "the existing data provider and to add this new data set"
    );
}

}

#endif
