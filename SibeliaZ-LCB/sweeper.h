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
	typedef std::pair<JunctionStorage::JunctionSequentialIterator, JunctionStorage::JunctionSequentialIterator> Template;

	bool operator < (const Template & a, const Template & b);
	

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
			//int score;
			bool hasNext;
			int32_t startIdx[2];
			int32_t endIdx[2];

			Instance(): hasNext(false)//, score(1)
			{

			}
			
			bool Valid(const JunctionStorage & storage, int64_t minBlockSize, size_t chrId0, size_t chrId1) const
			{
				for (size_t l = 0; l < 2; l++)
				{
					auto it = Start(storage, l == 0 ? chrId0 : chrId1, l);
					auto jt = End(storage, l == 0 ? chrId0 : chrId1, l);
					if (abs(jt.GetPosition() - it.GetPosition()) < minBlockSize)
					{
						return false;
					}
				}

				return true;
			}

			int32_t MainPosition(size_t chrId, const JunctionStorage & storage) const
			{
				return storage.GetIterator(chrId, endIdx[0] - 1).GetPosition();
			}

			static int32_t GetIdx(JunctionStorage::JunctionSequentialIterator & dummy)
			{
				if (dummy.IsPositiveStrand())
				{
					return dummy.GetIndex() + 1;
				}

				return -(int32_t(dummy.GetIndex()) + 1);
			}

			Instance(JunctionStorage::JunctionSequentialIterator & main, JunctionStorage::JunctionSequentialIterator & secondary)
			{
				startIdx[0] = endIdx[0] = GetIdx(main);
				startIdx[1] = endIdx[1] = GetIdx(secondary);
			}

			Instance(JunctionStorage::JunctionSequentialIterator & dummy)
			{
				endIdx[1] = GetIdx(dummy);
			}

			bool IsPositiveStrand() const
			{
				return endIdx[1] > 0;
			}

			bool operator < (const Instance & cmp) const
			{
				if (IsPositiveStrand() != cmp.IsPositiveStrand())
				{
					return IsPositiveStrand() < cmp.IsPositiveStrand();
				}

				if (endIdx[1] > 0)
				{
					return endIdx[1] < cmp.endIdx[1];
				}

				return endIdx[1] > cmp.endIdx[1];
			}

			JunctionStorage::JunctionSequentialIterator End(const JunctionStorage & storage, size_t chrId, size_t idx) const
			{
				if (idx == 0)
				{
					return storage.GetIterator(chrId, endIdx[0] - 1);
				}

				if (IsPositiveStrand())
				{
					return storage.GetIterator(chrId, endIdx[1] - 1);
				}

				return storage.GetIterator(chrId, -endIdx[1] - 1, false);
			}

			JunctionStorage::JunctionSequentialIterator Start(const JunctionStorage & storage, size_t chrId, size_t idx) const
			{
				if (idx == 0)
				{
					return storage.GetIterator(chrId, startIdx[0] - 1);
				}

				if (IsPositiveStrand())
				{
					return storage.GetIterator(chrId, startIdx[1] - 1);
				}

				return storage.GetIterator(chrId, -startIdx[1] - 1, false);
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

		Sweeper(JunctionStorage::JunctionSequentialIterator start) : start_(start)
		{

		}

		void Purge(int32_t lastPos, const JunctionStorage & storage, int64_t k, std::atomic<int64_t> & blocksFound, std::vector<BlockInstance> & blocksInstance, int32_t minBlockSize, int32_t maxBranchSize, std::vector<std::multiset<Sweeper::Instance> > & instance)
		{
			while (purge_.size() > 0)
			{
				int64_t chrId = purge_.front().first;
				int32_t diff = lastPos - (purge_.front().second->MainPosition(start_.GetChrId(), storage));
				if (diff >= maxBranchSize)
				{
					auto it = purge_.front().second;
					if (it->Valid(storage, minBlockSize, start_.GetChrId(), chrId) && !it->hasNext)
					{
						ReportBlock(blocksInstance, storage, chrId, k, blocksFound, *it);
					}
					
					purge_.pop_front();
					instance[chrId].erase(it);
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
			std::vector<std::multiset<Sweeper::Instance> > & instance)
		{
			size_t totalErases = 0;
			JunctionStorage::JunctionSequentialIterator successor[2];
			for (auto it = start_; it.Valid(); ++it)
			{
				for (JunctionStorage::JunctionIterator vit(it.GetVertexId()); vit.Valid(); ++vit)
				{
					Instance bestPrev;
					bool found = false;
					auto jt = vit.SequentialIterator();
					if (jt != it && jt.GetChrId() >= it.GetChrId())
					{
						int64_t chrId = jt.GetChrId();
						auto kt = instance[chrId].upper_bound(Instance(jt));
						if (kt != instance[chrId].begin())
						{
							if ((--kt)->IsPositiveStrand() == jt.IsPositiveStrand())
							{
								successor[0] = it;							
								successor[1] = jt;
								if (Compatible(storage, *kt, chrId, successor, maxBranchSize))
								{
									const_cast<Instance&>(*kt).hasNext = true;
									{
										bestPrev = *kt;
										found = true;
									}
								}
							}
						}

						if (!found)
						{
							purge_.push_back(std::make_pair(chrId, instance[chrId].insert(Instance(it, jt))));
						}
						else
						{
							Instance newUpdate(bestPrev);
							newUpdate.hasNext = false;
							newUpdate.endIdx[0] = newUpdate.GetIdx(it);
							newUpdate.endIdx[1] = newUpdate.GetIdx(jt);
							purge_.push_back(std::make_pair(chrId, instance[chrId].insert(newUpdate)));
						}
					}
				}

				Purge(it.GetPosition(), storage, k, blocksFound, blocksInstance, minBlockSize, maxBranchSize, instance);
			}

			Purge(INT32_MAX, storage, k, blocksFound, blocksInstance, minBlockSize, maxBranchSize, instance);
		}

	private:
		typedef std::multiset<Instance>::iterator InstanceIt;
		
		std::deque<std::pair<int64_t, InstanceIt> > purge_;
		JunctionStorage::JunctionSequentialIterator start_;

		bool Compatible(const JunctionStorage & storage, const Instance & inst, int64_t chrId1, const JunctionStorage::JunctionSequentialIterator succ[2], int64_t maxBranchSize) const
		{
			bool withinBubble = true;
			int64_t chrId0 = start_.GetChrId();
			JunctionStorage::JunctionSequentialIterator end[2] = { inst.End(storage, chrId0, 0), inst.End(storage, chrId1, 1) };
			bool validSuccessor = end[0].GetChar() == end[1].GetChar();
			for (size_t i = 0; i < 2; i++)
			{
				if (end[i] + 1 != succ[i])
				{
					validSuccessor = false;
				}

				if (abs(end[i].GetPosition() - succ[i].GetPosition()) >= maxBranchSize)
				{
					withinBubble = false;
				}
			}

			if (withinBubble || validSuccessor)
			{
				JunctionStorage::JunctionSequentialIterator start[2] = { inst.Start(storage, chrId0, 0), inst.Start(storage, chrId1, 1) };
				if(start[0].GetChrId() == start[1].GetChrId())
				{
					size_t startIdx1 = start[0].GetIndex();
					size_t endIdx1 = succ[0].GetIndex();
					size_t startIdx2 = min(start[1].GetIndex(), succ[1].GetIndex());
					size_t endIdx2 = max(start[1].GetIndex(), succ[1].GetIndex());
					if ((startIdx2 >= startIdx1 && startIdx2 <= endIdx1) || (startIdx1 >= startIdx2 && startIdx1 <= endIdx2))
					{
						return false;
					}
				}

				return true;
			}

			return false;
		}

		void ReportBlock(std::vector<BlockInstance> & blocksInstance, const JunctionStorage & storage, int64_t chrId1, int64_t k, std::atomic<int64_t> & blocksFound, const Instance & inst)
		{
			int64_t chrId[] = { start_.GetChrId(), chrId1 };
			int64_t currentBlock = ++blocksFound;
			for (size_t l = 0; l < 2; l++)
			{
				auto it = inst.Start(storage, chrId[l], l);
				auto jt = inst.End(storage, chrId[l], l);
				{
					if (jt.IsPositiveStrand())
					{
						blocksInstance.push_back(BlockInstance(+currentBlock, jt.GetChrId(), it.GetPosition(), jt.GetPosition() + k));
					}
					else
					{
						blocksInstance.push_back(BlockInstance(-currentBlock, jt.GetChrId(), jt.GetPosition() - k, it.GetPosition()));
					}
				}
			}
		}
		
	};
}

#endif

