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
		int32_t idx;
		int32_t chrId;
		uint32_t score;
		int32_t endPosition[2];
		int32_t startPosition[2];

		Instance() : hasNext(false),  score(1)
		{

		}

		Instance(int32_t chrId, int32_t idx) : chrId(chrId), idx(idx)
		{

		}

		bool Valid(int64_t minBlockSize) const
		{
			for (size_t l = 0; l < 2; l++)
			{
				if (abs(endPosition[l] - startPosition[l]) < minBlockSize)
				{
					return false;
				}
			}

			return true;
		}

		Instance(const Instance & inst, JunctionStorage::Iterator & it, JunctionStorage::Iterator & jt) : hasNext(false), score(inst.score + 1)
		{
			startPosition[0] = inst.startPosition[0];
			startPosition[1] = inst.startPosition[1];
			endPosition[0] = it.GetPosition();
			endPosition[1] = jt.GetPosition();
			
			idx = static_cast<int32_t>(jt.GetIndex());
			chrId = jt.IsPositiveStrand() ? (jt.GetChrId() + 1) : -(jt.GetChrId() + 1);
		}

		Instance(JunctionStorage::Iterator & it, JunctionStorage::Iterator & jt) : hasNext(false)
		{
			startPosition[0] = endPosition[0] = it.GetPosition();
			startPosition[1] = endPosition[1] = jt.GetPosition();

			idx = static_cast<int32_t>(jt.GetIndex());
			chrId = jt.IsPositiveStrand() ? (jt.GetChrId() + 1) : -(jt.GetChrId() + 1);
		}

		bool IsPositiveStrand() const
		{
			return endPosition[1] > 0;
		}

		bool operator == (const Instance & cmp) const
		{
			for (size_t i = 1; i < 2; i++)
			{
				if (startPosition[i] != cmp.startPosition[i] || endPosition[i] != cmp.endPosition[i])
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

	struct VertexEntry
	{
		int64_t vertexId;
		uint16_t pointerIdx;
		std::vector<Instance>* instance;

		VertexEntry() {}
		VertexEntry(int64_t vertexId, uint16_t pointerIdx, std::vector<Instance>* instance) : vertexId(vertexId), pointerIdx(pointerIdx), instance(instance)
		{

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
			isActive_.resize((chrSize >> 6) + 1, false);
		}

		void Add(Instance * inst, size_t idx)
		{
			uint64_t bit;
			uint64_t element;
			GetCoord(idx, element, bit);
			isActive_[element] |= uint64_t(1) << uint64_t(bit);
		}

		Instance* Retreive(const JunctionStorage & storage, std::vector<VertexEntry* >& lastPosEntry, std::vector<VertexEntry* > &lastNegEntry, int32_t maxBranchSize, size_t chr0, int32_t chr1idx)
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

		uint32_t Compatible(const Instance & inst, const JunctionStorage::Iterator succ[2], int64_t maxBranchSize) const
		{
			bool withinBubble = true;
			for (size_t i = 0; i < 2; i++)
			{
				if (abs(inst.endPosition[i] - succ[i].GetPosition()) >= maxBranchSize)
				{
					withinBubble = false;
				}
			}

			if (withinBubble)
			{
				if (succ[0].GetChrId() == succ[1].GetChrId())
				{
					size_t startPosition2 = min(abs(inst.startPosition[1]), abs(inst.endPosition[1]));
					size_t endPosition2 = max(abs(inst.startPosition[1]), abs(inst.endPosition[1]));
					if ((startPosition2 >= inst.startPosition[0] && startPosition2 <= inst.endPosition[0]) || (inst.startPosition[0] >= startPosition2 && inst.startPosition[0] <= endPosition2))
					{
						return false;
					}
				}

				if (withinBubble)
				{
					return 1;
				}
			}

			return 0;
		}


		uint32_t CompatibleExact(const Instance & inst, const JunctionStorage::Iterator succ[2], int64_t maxBranchSize) const
		{
			if (succ[0].GetChrId() == succ[1].GetChrId())
			{
				size_t startPosition2 = min(abs(inst.startPosition[1]), abs(inst.endPosition[1]));
				size_t endPosition2 = max(abs(inst.startPosition[1]), abs(inst.endPosition[1]));
				if ((startPosition2 >= inst.startPosition[0] && startPosition2 <= inst.endPosition[0]) || (inst.startPosition[0] >= startPosition2 && inst.startPosition[0] <= endPosition2))
				{
					return 0;
				}
			}

			return abs(inst.endPosition[0] - succ[0].GetPosition());
		}

		std::pair<Instance*, uint32_t> TryRetreiveExact(const JunctionStorage & storage,
			std::vector<VertexEntry* >& lastPosEntry,
			std::vector<VertexEntry* >&lastNegEntry,
			const JunctionStorage::Iterator succ[2],
			JunctionStorage::Iterator chr0Prev)
		{

			Instance* ret = 0;
			int64_t bestScore = 0;
			if (chr0Prev.Valid())
			{
				auto chr1Prev = succ[1];
				chr1Prev.DecInSequence();
				if (chr1Prev.Valid() && chr0Prev.GetChar() == chr1Prev.GetChar())
				{
					auto inst = GetMagicIndex(storage, lastPosEntry, lastNegEntry, chr1Prev.GetIndex());
					if (inst != 0)
					{
						auto gapScore = CompatibleExact(*inst, succ, 0);
						if (gapScore > 0)
						{
							return std::make_pair(inst, inst->score + gapScore);
						}
					}
				}
			}

			return std::pair<Instance*, uint32_t>(0, 0);
		}

		std::pair<Instance*, uint32_t> RetreiveBest(const JunctionStorage & storage,
			std::vector<VertexEntry* >& lastPosEntry,
			std::vector<VertexEntry* >&lastNegEntry,
			int32_t maxBranchSize,
			const JunctionStorage::Iterator succ[2])
		{
			
			uint64_t bit;
			uint64_t element;
			Instance* ret = 0;
			uint32_t bestScore = 0;
			GetCoord(succ[1].GetIndex(), element, bit);
			int32_t maxBranchSizeElement = (maxBranchSize >> 6) + 1;
			if (isPositiveStrand_)
			{
				bool go = true;
				int32_t elementLimit = max(0, int32_t(element) - maxBranchSizeElement);
				for (int32_t e = element; e >= elementLimit && go; e--)
				{
					auto mask = isActive_[e];
					if (e == element && bit < 63)
					{
						auto application = (uint64_t(1) << (bit + uint64_t(1)));
						mask = mask & (application - uint64_t(1));
					}

					while (mask > 0)
					{
#ifdef _MSC_VER
						uint64_t bit = 64 - __lzcnt64(mask);
#else
						uint64_t bit = 64 - __builtin_clzll(mask);
#endif
						auto idx = (e << 6) | (bit - 1);
						mask &= ~(uint64_t(1) << (bit - 1));

						auto inst = GetMagicIndex(storage, lastPosEntry, lastNegEntry, idx);
						if (inst != 0)
						{
							if (succ[1].GetPosition() - inst->endPosition[1] >= maxBranchSize)
							{
								go = false;
								break;
							}

							auto gapScore = Compatible(*inst, succ, maxBranchSize);
							if (gapScore > 0 && (ret == 0 || (inst->score + gapScore > bestScore)))
							{
								ret = inst;
								bestScore = inst->score + gapScore;
							}
						}
					}
				}
			}
			else
			{
				bool go = true;
				int32_t elementLimit = min(int32_t(isActive_.size()), int32_t(element) + maxBranchSizeElement);
				for (int32_t e = element; e < elementLimit && go; e++)
				{
					auto mask = isActive_[e];
					if (e == element)
					{
						auto application = (uint64_t(1) << bit) - uint64_t(1);
						mask = mask & (~application);
					}

					while (mask != 0)
					{
#ifdef _MSC_VER
						uint64_t bit = _tzcnt_u64(mask);
#else
						uint64_t bit = __builtin_ctzll(mask);
#endif
						auto idx = (e << 6) | bit;

						mask &= ~(uint64_t(1) << bit);

						auto inst = GetMagicIndex(storage, lastPosEntry, lastNegEntry, idx);
						if (inst != 0)
						{
							if (inst->endPosition[1] - succ[1].GetPosition() >= maxBranchSize)
							{
								go = false;
								break;
							}

							auto gapScore = Compatible(*inst, succ, maxBranchSize);
							if (gapScore > 0 && (ret == 0 || (inst->score + gapScore > bestScore)))
							{
								ret = inst;
								bestScore = inst->score + gapScore;
							}							
						}
					}
				}
			}

			return std::make_pair(ret, bestScore);
		}

		void Erase(Instance * inst, const JunctionStorage & storage, std::vector<VertexEntry* >& lastPosEntry, std::vector<VertexEntry* > &lastNegEntry, int32_t maxBranchSize, size_t chr0, int32_t chr1idx)
		{
			auto * currentInst = Retreive(storage, lastPosEntry, lastNegEntry, maxBranchSize, chr0, chr1idx);
			if (currentInst == inst)
			{
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
		

		Instance* GetMagicIndex(const JunctionStorage & storage, std::vector<VertexEntry* > & lastPosEntry, std::vector<VertexEntry* > & lastNegEntry, size_t chr1Idx) const
		{
			int64_t vid = storage.GetVertexId(chr1_, chr1Idx);
			if (!isPositiveStrand_)
			{
				vid = -vid;
			}

			VertexEntry* e = vid > 0 ? lastPosEntry[vid] : lastNegEntry[-vid];
			if (e == 0)
			{
				return 0;
			}

			auto chr1PointerIdx = storage.GetPointerIndex(chr1_, chr1Idx);
			auto magicIdx = chr1PointerIdx - e->pointerIdx - 1;
			if (magicIdx >= 0 && magicIdx < (*e->instance).size())
			{
				return &(*e->instance)[magicIdx];
			}

			return 0;
		}

		Instance* GetInstanceBefore(const JunctionStorage & storage, std::vector<VertexEntry* >& lastPosEntry, std::vector<VertexEntry* > &lastNegEntry, size_t chr0, uint64_t element, uint64_t mask)
		{
#ifdef _MSC_VER
			uint64_t bit = 64 - __lzcnt64(mask);
#else
			uint64_t bit = 64 - __builtin_clzll(mask);
#endif
			auto idx = (element << 6) | (bit - 1);
			return GetMagicIndex(storage, lastPosEntry, lastNegEntry, idx);
		}

		Instance* GetInstanceAfter(const JunctionStorage & storage, std::vector<VertexEntry* > &lastPosEntry, std::vector<VertexEntry* >& lastNegEntry, size_t chr0, uint64_t element, uint64_t mask)
		{
#ifdef _MSC_VER
			uint64_t bit = _tzcnt_u64(mask);
#else
			uint64_t bit = __builtin_ctzll(mask);
#endif
			auto idx = (element << 6) | bit;
			return GetMagicIndex(storage, lastPosEntry, lastNegEntry, idx);
		}

		void GetCoord(uint64_t idx, uint64_t & element, uint64_t & bit) const
		{
			bit = idx & ((uint64_t(1) << uint64_t(6)) - 1);
			element = idx >> 6;
		}
	};
}

#endif