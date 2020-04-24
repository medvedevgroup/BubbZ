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

		typedef std::multiset<Instance>::iterator InstanceIt;

		Sweeper(JunctionStorage::Iterator start, std::vector<VertexEntry* > & lastPosEntry_, std::vector<VertexEntry* > & lastNegEntry_) :
			start_(start), lastPosEntry_(lastPosEntry_), lastNegEntry_(lastNegEntry_)
		{

		}

		void Purge(JunctionStorage & storage,
			int32_t lastPos,
			int32_t k,
			std::atomic<int64_t> & blocksFound,
			std::vector<BlockInstance> & blocksInstance,
			int32_t minBlockSize,
			int32_t maxBranchSize,
			std::vector<std::vector<InstanceSet> > & instance,
			int64_t currentVid)
		{
			while (purge_.size() > 0)
			{
				if (purge_.front().instance->size() > 0)
				{
					int32_t diff = lastPos - (purge_.front().instance->front().endIdx[0]);
					if (diff >= maxBranchSize)
					{
						for (auto & it : *purge_.front().instance)
						{
							int64_t chrId = abs(it.chrId) - 1;
							size_t strand = it.chrId > 0 ? 0 : 1;
							bool hasNext = it.hasNext;
							if (it.Valid(minBlockSize) && !hasNext)
							{
								ReportBlock(blocksInstance, chrId, k, blocksFound, it);
							}

							instance[strand][chrId].Erase(&it, storage, lastPosEntry_, lastNegEntry_, maxBranchSize, start_.GetChrId(), it.idx);
						}
						
						NotifyPop(purge_.front());
						purge_.front().instance->clear();
						pool_.push_back(purge_.front().instance);
						purge_.pop_front();
					}
					else
					{
						break;
					}
				}
				else 
				{
					if (purge_.front().vertexId != currentVid)
					{
						NotifyPop(purge_.front());
						pool_.push_back(purge_.front().instance);
						purge_.pop_front();
					}
					else
					{
						break;
					}
				}								
			}
		}

		void Sweep(JunctionStorage & storage,
			int64_t minBlockSize,
			int64_t maxBranchSize,
			int64_t k,
			std::atomic<int64_t> & blocksFound,
			std::vector<BlockInstance> & blocksInstance,
			std::vector<std::vector<InstanceSet> > & instance)
		{
			
			for (size_t i = 0; i < maxBranchSize + 1; i++)
			{
				pool_.push_back(new std::vector<Instance>());
				pool_.back()->reserve(storage.GetAbundance());
			}

			JunctionStorage::Iterator itPrev;
			JunctionStorage::Iterator successor[2];
			for (auto it = start_; it.Valid(); it.Inc())
			{
				auto jt = it;
				purge_.push_back(VertexEntry(it.GetVertexId(), it.GetPointerIndex(), pool_.back()));
				pool_.pop_back();
				for (jt.Next(); jt.Valid(); jt.Next())
				{
					size_t idx = jt.GetIndex();
					int32_t chrId = jt.GetChrId();
					size_t strand = jt.IsPositiveStrand() ? 0 : 1;
					successor[0] = it;
					successor[1] = jt;

					auto kt = instance[strand][chrId].TryRetreiveExact(storage, lastPosEntry_, lastNegEntry_, successor, itPrev);
					if (kt.first == 0)
					{
						kt = instance[strand][chrId].RetreiveBest(storage, lastPosEntry_, lastNegEntry_, maxBranchSize, successor);
					}
				
					if (kt.first != 0)
					{
						const_cast<Instance&>(*kt.first).hasNext = true;
						Instance newUpdate(*kt.first, it, jt);
						newUpdate.score = kt.second;
						purge_.back().instance->push_back(newUpdate);
						instance[strand][chrId].Add(&purge_.back().instance->back(), idx);
					}
					else
					{
						purge_.back().instance->push_back(Instance(it, jt));
						instance[strand][chrId].Add(&purge_.back().instance->back(), idx);
					}
				
				}

				NotifyPush(purge_.back());
				Purge(storage, it.GetPosition(), k, blocksFound, blocksInstance, minBlockSize, maxBranchSize, instance, 0);
				itPrev = it;
			}

			Purge(storage, INT32_MAX, k, blocksFound, blocksInstance, minBlockSize, maxBranchSize, instance, 0);
			for (auto pt : pool_)
			{
				delete pt;
			}
		}


	private:
		JunctionStorage::Iterator start_;
		std::deque<VertexEntry> purge_;
		std::vector<std::vector<Instance>* > pool_;
		std::vector<VertexEntry* > & lastPosEntry_;
		std::vector<VertexEntry* > & lastNegEntry_;

		void NotifyPush(VertexEntry & e)
		{
			if (e.vertexId > 0)
			{
				lastPosEntry_[e.vertexId] = &e;
			}
			else
			{
				lastNegEntry_[-e.vertexId] = &e;
			}
		}

		void NotifyPop(VertexEntry & e)
		{
			if (e.vertexId > 0)
			{
				if (lastPosEntry_[e.vertexId] == &e)
				{
					lastPosEntry_[e.vertexId] = 0;
				}
			}
			else
			{
				if (lastNegEntry_[-e.vertexId] == &e)
				{
					lastNegEntry_[-e.vertexId] = 0;
				}
			}
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
					blocksInstance.push_back(BlockInstance(-currentBlock, chrId[l], -(inst.endIdx[l]) - k, -inst.startIdx[l]));
				}
			}
		}
		
	};
}

#endif

