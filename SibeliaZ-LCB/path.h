#ifndef _PATH_H_
#define _PATH_H_

#include <set>
#include <cassert>
#include <algorithm>
#include "distancekeeper.h"

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

		void Init(size_t chrSize)
		{
			instance_.resize(chrSize, 0);
			isActive_.resize(chrSize, false);
		}

		void Add(Instance * inst, size_t idx)
		{
			instance_[idx] = inst;
			isActive_[idx] = true;
		}

		Instance* Retreive(const JunctionStorage & storage, int32_t maxBranchSize, int32_t idx, bool isPositive)
		{
			if (isPositive)
			{
				int32_t limit = max(0, idx - maxBranchSize);
				for (; idx >= limit; idx--)
				{
					if (isActive_[idx])
					{
						return instance_[idx];
					}
				}
			}
			else
			{
				int32_t limit = min(isActive_.size(), idx + maxBranchSize);
				for (; idx < isActive_.size(); idx++)
				{
					if (isActive_[idx])
					{
						return instance_[idx];
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
				isActive_[idx] = false;
			}
		}

	private:
		std::vector<bool> isActive_;
		std::vector<Instance*> instance_;
	};
}

#endif