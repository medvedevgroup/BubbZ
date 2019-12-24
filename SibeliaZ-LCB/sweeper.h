#ifndef _SWEEPER_H_
#define _SWEEPER_H_

#include <set>
#include <deque>
#include <queue>
#include <cassert>
#include <algorithm>

#include <omp.h>

#include "path.h"

namespace Sibelia
{
	class BlockInstance
	{
	public:
		BlockInstance() {}
		BlockInstance(int id, const size_t chr, size_t start, size_t end) : id_(id), chr_(chr), start_(start), end_(end) {}
		void Reverse();
		int GetSignedBlockId() const;
		bool GetDirection() const;
		int GetBlockId() const;
		int GetSign() const;
		size_t GetChrId() const;
		size_t GetStart() const;
		size_t GetEnd() const;
		size_t GetLength() const;
		size_t GetConventionalStart() const;
		size_t GetConventionalEnd() const;
		std::pair<size_t, size_t> CalculateOverlap(const BlockInstance & instance) const;
		bool operator < (const BlockInstance & toCompare) const;
		bool operator == (const BlockInstance & toCompare) const;
		bool operator != (const BlockInstance & toCompare) const;
	private:
		int id_;
		size_t start_;
		size_t end_;
		size_t chr_;
	};

	class Sweeper
	{
	public:

		struct Instance
		{
			bool hasNext;
			bool parallelEnd;
			int32_t endIdx[2];
			int32_t startIdx[2];

			Instance(): hasNext(false), parallelEnd(false)//, score(1)
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
			
			Instance(const Instance & inst, JunctionStorage::Iterator & it, JunctionStorage::Iterator & jt): hasNext(false)
			{
				startIdx[0] = inst.startIdx[0];
				startIdx[1] = inst.startIdx[1];
				endIdx[0] = it.GetPosition();
				endIdx[1] = jt.GetPosition();
				parallelEnd = it.GetChar() == jt.GetChar();
			}

			Instance(JunctionStorage::Iterator & it, JunctionStorage::Iterator & jt): hasNext(false)
			{
				startIdx[0] = endIdx[0] = it.GetPosition();
				startIdx[1] = endIdx[1] = jt.GetPosition();
				parallelEnd = it.GetChar() == jt.GetChar();
			}

			Instance(JunctionStorage::Iterator & dummy)
			{
				endIdx[1] = dummy.GetPosition();
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

		typedef std::multiset<Instance>::iterator InstanceIt;

		Sweeper(JunctionStorage::Iterator start) : start_(start)
		{

		}

		void Purge(int32_t lastPos, int64_t k, std::atomic<int64_t> & blocksFound, std::vector<BlockInstance> & blocksInstance, int32_t minBlockSize, int32_t maxBranchSize, std::vector<std::vector<std::multiset<Sweeper::Instance> > > & instance)
		{
			while (purge_.size() > 0)
			{
				int64_t chrId = abs(purge_.front().first) - 1;
				size_t strand = purge_.front().first > 0 ? 0 : 1;
				int32_t diff = lastPos - (purge_.front().second->endIdx[0]);
				if (diff >= maxBranchSize)
				{
					auto it = purge_.front().second;
					bool hasNext = it->hasNext;
					if (it->Valid(minBlockSize) && !hasNext)
					{
						ReportBlock(blocksInstance, chrId, k, blocksFound, *it);
					}
					
					purge_.pop_front();
					instance[strand][chrId].erase(it);
				}
				else
				{
					break;
				}
			}
		}

		void Sweep(JunctionStorage & storage,
			int64_t minBlockSize,
			int64_t maxBranchSize,
			int64_t k,
			std::atomic<int64_t> & blocksFound,
			std::vector<BlockInstance> & blocksInstance,
			std::vector<std::vector<std::multiset<Sweeper::Instance> > > & instance)
		{
			JunctionStorage::Iterator successor[2];
			for (auto it = start_; it.Valid(); it.Inc())
			{
				auto jt = it;
				for (jt.Next(); jt.Valid(); jt.Next())
				{
					Instance bestPrev;
					bool found = false;
					int64_t chrId = jt.GetChrId();
					size_t strand = jt.IsPositiveStrand() ? 0 : 1;
					auto kt = instance[strand][chrId].upper_bound(Instance(jt));
					if (kt != instance[strand][chrId].begin())
					{
						--kt;
						successor[0] = it;							
						successor[1] = jt;
						if (Compatible(*kt, chrId, successor, maxBranchSize))
						{
							const_cast<Instance&>(*kt).hasNext = true;
							bestPrev = *kt;
							found = true;
						}
					}

					if (!found)
					{
						auto lt = instance[strand][chrId].insert(Instance(it, jt));
						purge_.push_back(std::make_pair(strand == 0 ? (chrId + 1) : -(chrId + 1), lt));
					}
					else
					{
						Instance newUpdate(bestPrev, it, jt);
						auto lt = instance[strand][chrId].insert(newUpdate);
						purge_.push_back(std::make_pair(strand == 0 ? (chrId + 1) : -(chrId + 1), lt));
					}
				}

				Purge(it.GetPosition(), k, blocksFound, blocksInstance, minBlockSize, maxBranchSize, instance);
			}

			Purge(INT32_MAX, k, blocksFound, blocksInstance, minBlockSize, maxBranchSize, instance);
		}

	private:
		
		std::deque<std::pair<int64_t, InstanceIt> > purge_;
		JunctionStorage::Iterator start_;

		bool Compatible(const Instance & inst, int64_t chrId1, const JunctionStorage::Iterator succ[2], int64_t maxBranchSize) const
		{
			bool withinBubble = true;
			int64_t chrId0 = start_.GetChrId();
			bool validSuccessor = inst.parallelEnd;
			for (size_t i = 0; i < 2; i++)
			{
				if (inst.endIdx[i] != succ[i].PreviousPosition())
				{					
					validSuccessor = false;
				}

				if (abs(inst.endIdx[i] - succ[i].GetPosition()) >= maxBranchSize)
				{					
					withinBubble = false;
				}
			}
			
			if (withinBubble || validSuccessor)
			{
				if(chrId0 == chrId1)
				{
					size_t startIdx2 = min(abs(inst.startIdx[1]), abs(inst.endIdx[1]));
					size_t endIdx2 = max(abs(inst.startIdx[1]), abs(inst.endIdx[1]));
					if ((startIdx2 >= inst.startIdx[0] && startIdx2 <= inst.endIdx[0]) || (inst.startIdx[0] >= startIdx2 && inst.startIdx[0] <= endIdx2))
					{					
						return false;
					}
				}

				return true;
			}
			
			return false;
		}

		void ReportBlock(std::vector<BlockInstance> & blocksInstance, int64_t chrId1, int64_t k, std::atomic<int64_t> & blocksFound, const Instance & inst)
		{
			int64_t chrId[] = { start_.GetChrId(), chrId1 };
			int64_t currentBlock = ++blocksFound;
			for (size_t l = 0; l < 2; l++)
			{
				if (inst.endIdx[l] >= 0)
				{
					blocksInstance.push_back(BlockInstance(+currentBlock, chrId[l], inst.startIdx[l], inst.endIdx[l] + k));
				}
				else
				{
					blocksInstance.push_back(BlockInstance(-currentBlock, chrId[l], -(inst.endIdx[l] - k), -inst.startIdx[l]));
				}
			}
		}
		
	};
}

#endif

