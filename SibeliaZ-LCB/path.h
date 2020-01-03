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
			return endIdx[1] < cmp.endIdx[1];
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
			instance_.resize(chrSize, 0);
			isActive_.resize((chrSize >> 6) + 1, false);
		}

		void Add(Instance * inst, size_t idx)
		{
		instance_[idx] = inst;
		uint64_t bit;
		uint64_t element;
		GetCoord(idx, element, bit);
		isActive_[element] |= uint64_t(1) << uint64_t(bit);
		}

		Instance* Retreive(const JunctionStorage & storage, std::deque<std::pair<size_t, std::vector<Instance>* > > & q, int32_t maxBranchSize, size_t chr0, size_t chr0Bound, int32_t idx, bool isPositive)
		{
			uint64_t bit;
			uint64_t element;
			GetCoord(idx, element, bit);
			int32_t maxBranchSizeElement = (maxBranchSize >> 6) + 1;
			if (isPositive)
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
						return GetInstanceBefore(storage, q, chr0, chr0Bound, e, mask);
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
						return GetInstanceAfter(storage, q, chr0, chr0Bound, e, mask);
					}
				}

			}

			return 0;
		}

		void Erase(Instance * inst, size_t idx)
		{
			if (instance_[idx] == inst)
			{
				instance_[idx] = 0;
				uint64_t bit;
				uint64_t element;
				GetCoord(idx, element, bit);
				isActive_[element] &= ~(uint64_t(1) << bit);
			}
		}

	private:
		size_t chr1_;
		bool isPositiveStrand_;
		std::vector<uint64_t> isActive_;
		std::vector<Instance*> instance_;

		Instance* GetMagicIndex(const JunctionStorage & storage, std::deque<std::pair<size_t, std::vector<Instance>* > > & q, size_t chr0, size_t chr0Bound, size_t chr1Idx) const
		{
			int64_t vid = storage.GetVertexId(chr1_, chr1Idx);
			if (isPositiveStrand_)
			{
				vid = -vid;
			}

			size_t pointerChr0Idx;
			size_t pointerChr1Idx;
			auto & pointer = storage.Pointers(vid);
			for (size_t i = 0; i < pointer.size(); i++)
			{
				if (abs(pointer[i].chrId) == chr0 && pointer[i].idx <= chr0Bound)
				{
					pointerChr0Idx = i;
				}

				if (abs(pointer[i].chrId) == chr1_ && pointer[i].idx == chr1Idx)
				{
					pointerChr1Idx = i;
				}
			}

			auto qIndex = pointer[pointerChr0Idx].idx - q.front().first;
			auto iIndex = pointerChr1Idx - pointerChr0Idx - 1;
			if (iIndex >= q[qIndex].second->size())
			{
				return 0;
			}

			auto & val = (*q[qIndex].second)[iIndex];
			if (val.IsPositiveStrand() != isPositiveStrand_ || abs(val.chrId) - 1 != chr1_)
			{
				return 0;
			}

			return &val;
		}

		Instance* GetInstanceBefore(const JunctionStorage & storage, std::deque<std::pair<size_t, std::vector<Instance>* > > & q, size_t chr0, size_t chr0Bound, uint64_t element, uint64_t mask)
		{
#ifdef _MSC_VER
			uint64_t bit = 64 - __lzcnt64(mask);
#else
			uint64_t bit = 64 - __builtin_clzll(mask);
#endif
			auto idx = (element << 6) | (bit - 1);
			/*
			if (GetMagicIndex(storage, q, chr0, chr0Bound, idx) != instance_[idx])
			{
				GetMagicIndex(storage, q, chr0, chr0Bound, idx);
//				abort();
			}

			//assert(GetMagicIndex(storage, q, chr0, chr0Bound, idx) == instance_[idx]);

			//return instance_[idx];
			*/
			return GetMagicIndex(storage, q, chr0, chr0Bound, idx);
		}

		Instance* GetInstanceAfter(const JunctionStorage & storage, std::deque<std::pair<size_t, std::vector<Instance>* > > & q, size_t chr0, size_t chr0Bound, uint64_t element, uint64_t mask)
		{
#ifdef _MSC_VER
			uint64_t bit = _tzcnt_u64(mask);
#else
			uint64_t bit = __builtin_ctzll(mask);
#endif
			auto idx = (element << 6) | bit;
			return GetMagicIndex(storage, q, chr0, chr0Bound, idx);
		}

		void GetCoord(uint64_t idx, uint64_t & element, uint64_t & bit) const
		{
			bit = idx & ((uint64_t(1) << uint64_t(6)) - 1);
			element = idx >> 6;
		}


	};
}

#endif