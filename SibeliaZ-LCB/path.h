#ifndef _PATH_H_
#define _PATH_H_

#include <set>
#include <cassert>
#include <algorithm>
#include "distancekeeper.h"


#ifdef _MSC_VER
	#include <intrin.h>
#else

#endif

namespace Sibelia
{
	struct Instance
	{
		bool hasNext;
		bool parallelEnd;
		int32_t idx;
		int32_t chrId;
		int32_t endIdx[2];
		int32_t startIdx[2];

		Instance() : hasNext(false), parallelEnd(false)//, score(1)
		{

		}

		Instance(int32_t chrId, int32_t idx) : chrId(chrId), idx(idx)
		{

		}

		bool Valid(int64_t minBlockSize) const
		{
			for (size_t l = 0; l < 2; l++)
			{
				if (abs(endIdx[l] - startIdx[l]) < minBlockSize)
				{
					return false;
				}
			}

			return true;
		}

		Instance(const Instance & inst, JunctionStorage::Iterator & it, JunctionStorage::Iterator & jt) : hasNext(false)
		{
			startIdx[0] = inst.startIdx[0];
			startIdx[1] = inst.startIdx[1];
			endIdx[0] = it.GetPosition();
			endIdx[1] = jt.GetPosition();
			parallelEnd = it.GetChar() == jt.GetChar();
			
			idx = jt.GetIndex();
			chrId = jt.IsPositiveStrand() ? (jt.GetChrId() + 1) : -(jt.GetChrId() + 1);
		}

		Instance(JunctionStorage::Iterator & it, JunctionStorage::Iterator & jt) : hasNext(false)
		{
			startIdx[0] = endIdx[0] = it.GetPosition();
			startIdx[1] = endIdx[1] = jt.GetPosition();
			parallelEnd = it.GetChar() == jt.GetChar();

			idx = jt.GetIndex();
			chrId = jt.IsPositiveStrand() ? (jt.GetChrId() + 1) : -(jt.GetChrId() + 1);
		}

		bool IsPositiveStrand() const
		{
			return endIdx[1] > 0;
		}

		bool operator < (const Instance & cmp) const
		{
			auto id1 = abs(chrId);
			auto id2 = abs(cmp.chrId);
			if (id1 != id2)
			{
				return id1 < id2;
			}

			return idx < cmp.idx;
		}

		bool operator == (const Instance & cmp) const
		{
			for (size_t i = 1; i < 2; i++)
			{
				if (startIdx[i] != cmp.startIdx[i] || endIdx[i] != cmp.endIdx[i])
				{
					return false;
				}
			}

			return true;
		}

		bool operator != (const Instance & cmp) const
		{
			return !(*this == cmp);
		}
	};

	class InstanceSet
	{
	public:
		InstanceSet()
		{

		}

		void Init(size_t chr1, bool isPositiveStrand, size_t chrSize)
		{
			chr1_ = chr1;
			isPositiveStrand_ = isPositiveStrand;
		//	instance_.resize(chrSize, 0);
			isActive_.resize((chrSize >> 6) + 1, false);
		}

		void Add(Instance * inst, size_t idx)
		{
		//	instance_[idx] = inst;
			uint64_t bit;
			uint64_t element;
			GetCoord(idx, element, bit);
			isActive_[element] |= uint64_t(1) << uint64_t(bit);
		}

		Instance* Retreive(const JunctionStorage & storage, std::vector<std::vector<Instance>* >& lastPosEntry, std::vector<std::vector<Instance>* > &lastNegEntry, int32_t maxBranchSize, size_t chr0, int32_t chr1idx)
		{
			uint64_t bit;
			uint64_t element;
			GetCoord(chr1idx, element, bit);
			int32_t maxBranchSizeElement = (maxBranchSize >> 6) + 1;
			if (isPositiveStrand_)
			{
				int32_t elementLimit = max(0, int32_t(element) - maxBranchSizeElement);
				for (int32_t e = element; e >= elementLimit; e--)
				{
					auto mask = isActive_[e];
					if (e == element && bit < 63)
					{
						auto application = (uint64_t(1) << (bit + uint64_t(1)));
						mask = mask & (application - uint64_t(1));
					}

					if (mask != 0)
					{
						return GetInstanceBefore(storage, lastPosEntry, lastNegEntry, chr0, e, mask);
					}
				}
			}
			else
			{
				int32_t elementLimit = min(int32_t(isActive_.size()), int32_t(element) + maxBranchSizeElement);
				for (int32_t e = element; e < elementLimit; e++)
				{
					auto mask = isActive_[e];
					if (e == element)
					{
						auto application = (uint64_t(1) << bit) - uint64_t(1);
						mask = mask & (~application);
					}

					if (mask != 0)
					{
						return GetInstanceAfter(storage, lastPosEntry, lastNegEntry, chr0, e, mask);
					}
				}
			}

			return 0;
		}

		void Erase(Instance * inst, const JunctionStorage & storage, std::vector<std::vector<Instance>* >& lastPosEntry, std::vector<std::vector<Instance>* > &lastNegEntry, int32_t maxBranchSize, size_t chr0, int32_t chr1idx)
		{
			auto * currentInst = Retreive(storage, lastPosEntry, lastNegEntry, maxBranchSize, chr0, chr1idx);
			if (currentInst == inst)
			{
				//instance_[idx] = 0;
				uint64_t bit;
				uint64_t element;
				GetCoord(chr1idx, element, bit);
				isActive_[element] &= ~(uint64_t(1) << bit);
			}
		}

		

	private:
		size_t chr1_;
		bool isPositiveStrand_;
		std::vector<uint64_t> isActive_;
		//std::vector<Instance*> instance_;
		

		Instance* GetMagicIndex(const JunctionStorage & storage, std::vector<std::vector<Instance>* > & lastPosEntry, std::vector<std::vector<Instance>* > & lastNegEntry, size_t chr0, size_t chr1Idx) const
		{
			int64_t vid = storage.GetVertexId(chr1_, chr1Idx);
			if (!isPositiveStrand_)
			{
				vid = -vid;
			}

			std::vector<Instance>* instance = vid > 0 ? lastPosEntry[vid] : lastNegEntry[-vid];
			if (instance == 0)
			{
				return 0;
			}

			auto & val = *std::lower_bound(instance->begin(), instance->end(), Instance(chr1_ + 1, chr1Idx));
			return &val;
		}

		Instance* GetInstanceBefore(const JunctionStorage & storage, std::vector<std::vector<Instance>* >& lastPosEntry, std::vector<std::vector<Instance>* > &lastNegEntry, size_t chr0, uint64_t element, uint64_t mask)
		{
#ifdef _MSC_VER
			uint64_t bit = 64 - __lzcnt64(mask);
#else
			uint64_t bit = 64 - __builtin_clzll(mask);
#endif
			auto idx = (element << 6) | (bit - 1);
			return GetMagicIndex(storage, lastPosEntry, lastNegEntry, chr0, idx);
		}

		Instance* GetInstanceAfter(const JunctionStorage & storage, std::vector<std::vector<Instance>* > &lastPosEntry, std::vector<std::vector<Instance>* >& lastNegEntry, size_t chr0, uint64_t element, uint64_t mask)
		{
#ifdef _MSC_VER
			uint64_t bit = _tzcnt_u64(mask);
#else
			uint64_t bit = __builtin_ctzll(mask);
#endif
			auto idx = (element << 6) | bit;
			return GetMagicIndex(storage, lastPosEntry, lastNegEntry, chr0, idx);
		}

		void GetCoord(uint64_t idx, uint64_t & element, uint64_t & bit) const
		{
			bit = idx & ((uint64_t(1) << uint64_t(6)) - 1);
			element = idx >> 6;
		}


	};
}

#endif