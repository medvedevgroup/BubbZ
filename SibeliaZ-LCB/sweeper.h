#ifndef _SWEEPER_H_
#define _SWEEPER_H_

#include <set>
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
			int score;
			bool hasNext;
			JunctionStorage::JunctionSequentialIterator start[2];
			JunctionStorage::JunctionSequentialIterator end[2];

			Instance(): hasNext(false), score(1)
			{

			}

			int64_t Score(JunctionStorage::JunctionSequentialIterator it, JunctionStorage::JunctionSequentialIterator jt) const
			{
				return abs(start[0].GetPosition() - it.GetPosition()) + abs(start[1].GetPosition() - jt.GetPosition());
			}

			bool Valid(int64_t minBlockSize) const
			{
				for (size_t i = 0; i < 2; i++)
				{
					if (abs(end[i].GetPosition() - start[i].GetPosition()) < minBlockSize)
					{
						return false;
					}
				}

				return true;
			}

			Instance(JunctionStorage::JunctionSequentialIterator dummy)
			{
				end[1] = dummy;
			}

			bool operator < (const Instance & cmp) const
			{
				if (end[1].IsPositiveStrand() != cmp.end[1].IsPositiveStrand())
				{
					return end[1].IsPositiveStrand() < cmp.end[1].IsPositiveStrand();
				}

				int64_t idx1 = end[1].GetIndex();
				int64_t idx2 = cmp.end[1].GetIndex();

				if (end[1].IsPositiveStrand())
				{
					return idx1 < idx2;
				}

				return idx1 > idx2;
			}

			bool operator == (const Instance & cmp) const
			{
				for (size_t i = 0; i < 2; i++)
				{
					if (start[i] != cmp.start[i] || end[i] != cmp.end[i])
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
					bool found = false;
					int64_t bestScore = -1;
					std::multiset<Sweeper::Instance>::iterator bestPrev;
					auto jt = vit.SequentialIterator();
					if (jt != it && jt.GetChrId() >= it.GetChrId())
					{
						int64_t chrId = jt.GetChrId();
						auto kt = instance[chrId].upper_bound(Instance(jt));
						if (kt != instance[chrId].begin())
						{
							do
							{
								--kt;
								if (kt->end[1].GetChrId() == jt.GetChrId() && kt->end[1].IsPositiveStrand() == jt.IsPositiveStrand())
								{
									successor[0] = it;
									successor[1] = jt;
									if (Compatible(*kt, successor, maxBranchSize))
									{
										const_cast<Instance&>(*kt).hasNext = true;
										if (kt->score > bestScore)
										{
											found = true;
											bestPrev = kt;
											bestScore = kt->score;
										}
									}

									if (jt.GetPosition() - kt->end[1].GetPosition() > maxBranchSize)
									{
										break;
									}
								}
								else
								{
									break;
								}

							} while (kt != instance[chrId].begin());
						}

						if (!found)
						{
							Instance newInstance;
							newInstance.start[0] = newInstance.end[0] = it;
							newInstance.start[1] = newInstance.end[1] = jt;
							purge_.push(instance[chrId].insert(newInstance));
						}
						else
						{
							Instance newUpdate(*bestPrev);
							newUpdate.hasNext = false;
							newUpdate.end[0] = it;
							newUpdate.end[1] = jt;
							newUpdate.score += 1;
							purge_.push(instance[chrId].insert(newUpdate));
						}
					}
				}

				while (purge_.size() > 0)
				{
					int64_t diff = it.GetPosition() - (purge_.top())->end[0].GetPosition();
					if (diff >= maxBranchSize)
					{
						auto it = purge_.top();
						if (it->Valid(minBlockSize) && !it->hasNext)
						{
							ReportBlock(blocksInstance, k, blocksFound, *it);
						}

						int64_t chrId = it->start[1].GetChrId();
						purge_.pop();
						instance[chrId].erase(it);
					}
					else
					{
						break;
					}
				}
			}

			for(;purge_.size() > 0; purge_.pop())
			{
				if (purge_.top()->Valid(minBlockSize) && !purge_.top()->hasNext)
				{
					ReportBlock(blocksInstance, k, blocksFound, *purge_.top());
				}

				int64_t chrId = purge_.top()->start[1].GetChrId();
				instance[chrId].erase(purge_.top());
			}
		}

	private:
		typedef std::multiset<Instance>::iterator InstanceIt;
		struct IteratorCmp
		{
			bool operator() (const InstanceIt & a, const InstanceIt & b) const
			{
				return a->end[0].GetPosition() > b->end[0].GetPosition();
			}
		};

		std::priority_queue<InstanceIt, std::vector<InstanceIt>, IteratorCmp> purge_;
		JunctionStorage::JunctionSequentialIterator start_;

		bool Compatible(const Instance & inst, const JunctionStorage::JunctionSequentialIterator succ[2], int64_t maxBranchSize) const
		{
			bool withinBubble = true;
			bool validSuccessor = inst.end[0].GetChar() == inst.end[1].GetChar();
			for (size_t i = 0; i < 2; i++)
			{
				if (inst.end[i] + 1 != succ[i])
				{
					validSuccessor = false;
				}

				if (abs(inst.end[i].GetPosition() - succ[i].GetPosition()) >= maxBranchSize)
				{
					withinBubble = false;
				}
			}

			if (withinBubble || validSuccessor)
			{
				if(inst.start[0].GetChrId() == inst.start[1].GetChrId())
				{
					size_t startIdx1 = inst.start[0].GetIndex();
					size_t endIdx1 = succ[0].GetIndex();
					size_t startIdx2 = min(inst.start[1].GetIndex(), succ[1].GetIndex());
					size_t endIdx2 = max(inst.start[1].GetIndex(), succ[1].GetIndex());
					if ((startIdx2 >= startIdx1 && startIdx2 <= endIdx1) || (startIdx1 >= startIdx2 && startIdx1 <= endIdx2))
					{
						return false;
					}
				}

				return true;
			}

			return false;
		}

		void ReportBlock(std::vector<BlockInstance> & blocksInstance, int64_t k, std::atomic<int64_t> & blocksFound, const Instance & inst)
		{
			int64_t currentBlock = ++blocksFound;
			for (size_t l = 0; l < 2; l++)
			{
				auto it = inst.start[l];
				auto jt = inst.end[l];
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

