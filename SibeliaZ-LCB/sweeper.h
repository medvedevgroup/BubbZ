#ifndef _SWEEPER_H_
#define _SWEEPER_H_

#include <set>
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
			bool hasNext;
			JunctionStorage::JunctionSequentialIterator start[2];
			JunctionStorage::JunctionSequentialIterator end[2];

			Instance(): hasNext(false)
			{

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
				if (end[1].GetChrId() != cmp.end[1].GetChrId())
				{
					return end[1].GetChrId() < cmp.end[1].GetChrId();
				}

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

		bool InstanceFound(const Instance & newInstance)
		{
			auto er = instance_.equal_range(newInstance);
			for (auto it = er.first; it != er.second; ++it)
			{
				if (*it == newInstance)
				{
					return true;
				}
			}

			return false;
		}

		void Sweep(int64_t minBlockSize, int64_t maxBranchSize, int64_t k, std::atomic<int64_t> & blocksFound, std::vector<BlockInstance> & blocksInstance, omp_lock_t & outMutex)
		{
			JunctionStorage::JunctionSequentialIterator successor[2];
			for (auto it = start_; it.Valid(); ++it)
			{
				{
					//std::cout << instance_.size() << ' ';
				}

				for (JunctionStorage::JunctionIterator vit(it.GetVertexId()); vit.Valid(); ++vit)
				{
					bool found = false;
					auto jt = vit.SequentialIterator();
					std::vector<Instance> newInstance;
					if (jt != it && jt.GetChrId() >= it.GetChrId())
					{
						auto kt = instance_.upper_bound(Instance(jt));
						if (kt != instance_.begin())
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
										found = true;
										Instance newUpdate(*kt);
										newUpdate.end[0] = it;
										newUpdate.end[1] = jt;
										if (!InstanceFound(newUpdate) && std::find(newInstance.begin(), newInstance.end(), newUpdate) == newInstance.end())
										{
											newInstance.push_back(newUpdate);
											(const_cast<Instance&>(*kt)).hasNext = true;
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

							} while (kt != instance_.begin());
						}

						if (!found)
						{
							Instance novel;
							novel.start[0] = novel.end[0] = it;
							novel.start[1] = novel.end[1] = jt;
							purge_.insert(instance_.insert(novel));
						}
						
						while (purge_.size() > 0)
						{
							if (it.GetPosition() - (*purge_.begin())->end[0].GetPosition() >= maxBranchSize)
							{
								auto it = *purge_.begin();
								if (it->Valid(minBlockSize) && !it->hasNext)
								{
									ReportBlock(blocksInstance, outMutex, k, blocksFound, *it);
								}

								purge_.erase(purge_.begin());
								instance_.erase(*it);
							}
							else
							{
								break;
							}
						}

						for (auto & inst : newInstance)
						{
							purge_.insert(instance_.insert(inst));
						}

						newInstance.clear();
					}
				}
			}

			for (auto inst = instance_.begin(); inst != instance_.end(); ++inst)
			{
				if (inst->Valid(minBlockSize) && !inst->hasNext)
				{
					ReportBlock(blocksInstance, outMutex, k, blocksFound, *inst);
				}
			}
		}

	private:
		typedef std::multiset<Instance>::iterator InstanceIt;
		struct IteratorCmp
		{
			bool operator() (const InstanceIt & a, const InstanceIt & b) const
			{
				return a->end[0].GetPosition() < b->end[0].GetPosition();
			}
		};

		std::multiset<Instance> instance_;
		std::multiset<InstanceIt, IteratorCmp> purge_;
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

		void ReportBlock(std::vector<BlockInstance> & blocksInstance, omp_lock_t & outMutex, int64_t k, std::atomic<int64_t> & blocksFound, const Instance & inst)
		{
			int64_t currentBlock = ++blocksFound;
			omp_set_lock(&outMutex);
			for (size_t l = 0; l < 2; l++)
			{
				auto it = inst.start[l];
				auto jt = inst.end[l];
				if (jt.IsPositiveStrand())
				{
					blocksInstance.push_back(BlockInstance(+currentBlock, jt.GetChrId(), it.GetPosition(), jt.GetPosition() + k));
				}
				else
				{
					blocksInstance.push_back(BlockInstance(-currentBlock, jt.GetChrId(), jt.GetPosition() - k, it.GetPosition()));
				}
			}
			omp_unset_lock(&outMutex);
		}
		
	};
}

#endif

