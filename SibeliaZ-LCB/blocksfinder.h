#ifndef _TRASERVAL_H_
#define _TRAVERSAL_H_

//#define _DEBUG_OUT_

#include <set>
#include <map>
#include <list>
#include <ctime>
#include <queue>
#include <iterator>
#include <cassert>
#include <numeric>
#include <sstream>
#include <iostream>
#include <functional>
#include <unordered_map>

#include <tbb/parallel_for.h>

#include "sweeper.h"

namespace Sibelia
{
	extern const std::string DELIMITER;
	extern const std::string VERSION;

	namespace
	{
		const bool COVERED = true;
		typedef std::vector<BlockInstance> BlockList;
		typedef std::pair<size_t, std::vector<BlockInstance> > GroupedBlock;
		typedef std::vector<GroupedBlock> GroupedBlockList;
		bool ByFirstElement(const GroupedBlock & a, const GroupedBlock & b)
		{
			return a.first < b.first;
		}

		std::string IntToStr(size_t x)
		{
			std::stringstream ss;
			ss << x;
			return ss.str();
		}

		template<class Iterator1, class Iterator2>
		void CopyN(Iterator1 it, size_t count, Iterator2 out)
		{
			for (size_t i = 0; i < count; i++)
			{
				*out++ = *it++;
			}
		}

		template<class Iterator>
		Iterator AdvanceForward(Iterator it, size_t step)
		{
			std::advance(it, step);
			return it;
		}

		template<class Iterator>
		Iterator AdvanceBackward(Iterator it, size_t step)
		{
			for (size_t i = 0; i < step; i++)
			{
				--it;
			}

			return it;
		}


		typedef std::pair<size_t, size_t> IndexPair;
		template<class T, class F, class It>
		void GroupBy(std::vector<T> & store, F pred, It out)
		{
			sort(store.begin(), store.end(), pred);
			for (size_t now = 0; now < store.size();)
			{
				size_t prev = now;
				for (; now < store.size() && !pred(store[prev], store[now]); now++);
				*out++ = std::make_pair(prev, now);
			}
		}

		template<class F>
		bool CompareBlocks(const BlockInstance & a, const BlockInstance & b, F f)
		{
			return (a.*f)() < (b.*f)();
		}

		template<class F>
		bool EqualBlocks(const BlockInstance & a, const BlockInstance & b, F f)
		{
			return f(a) == f(b);
		}

		template<class Iterator, class F, class ReturnType>
		struct FancyIterator : public std::iterator<std::forward_iterator_tag, ReturnType>
		{
		public:
			FancyIterator& operator++()
			{
				++it;
				return *this;
			}

			FancyIterator operator++(int)
			{
				FancyIterator ret(*this);
				++(*this);
				return ret;
			}

			bool operator == (FancyIterator toCompare) const
			{
				return it == toCompare.it;
			}

			bool operator != (FancyIterator toCompare) const
			{
				return !(*this == toCompare);
			}

			ReturnType operator * ()
			{
				return f(*it);
			}

			FancyIterator() {}
			FancyIterator(Iterator it, F f) : it(it), f(f) {}

		private:
			F f;
			Iterator it;
		};

		template<class Iterator, class F, class ReturnType>
		FancyIterator<Iterator, F, ReturnType> CFancyIterator(Iterator it, F f, ReturnType)
		{
			return FancyIterator<Iterator, F, ReturnType>(it, f);
		}

	}

	bool compareById(const BlockInstance & a, const BlockInstance & b);
	bool compareByChrId(const BlockInstance & a, const BlockInstance & b);
	bool compareByStart(const BlockInstance & a, const BlockInstance & b);

	void CreateOutDirectory(const std::string & path);

	class BlocksFinder
	{
	public:

		BlocksFinder(JunctionStorage & storage, size_t k) : storage_(storage), k_(k)
		{
			progressCount_ = 50;
			scoreFullChains_ = true;
		}

		static bool DegreeCompare(const JunctionStorage & storage, int64_t v1, int64_t v2)
		{
			return storage.GetInstancesCount(v1) > storage.GetInstancesCount(v2);
		}

		void Split(std::string & source, std::vector<std::string> & result)
		{
			std::stringstream ss;
			ss << source;
			result.clear();
			while (ss >> source)
			{
				result.push_back(source);
			}
		}

		void FindBlocks(int64_t minBlockSize, int64_t maxBranchSize, int64_t maxFlankingSize, int64_t lookingDepth, int64_t sampleSize, int64_t threads, const std::string & debugOut)
		{
			blocksFound_ = 0;
			minBlockSize_ = minBlockSize;
			maxBranchSize_ = maxBranchSize;
			maxFlankingSize_ = maxFlankingSize;

			using namespace std::placeholders;
			tbb::task_scheduler_init init(static_cast<int>(threads));

			std::vector<int64_t> shuffle;
			for (int64_t v = -storage_.GetVerticesNumber() + 1; v < storage_.GetVerticesNumber(); v++)
			{
				for (JunctionStorage::JunctionIterator it(v); it.Valid(); ++it)
				{
					if (it.IsPositiveStrand())
					{
						shuffle.push_back(v);
						break;
					}
				}
			}

			time_t mark = time(0);
			count_ = 0;
			
			starter_ = 0;
			progressCount_ = 0;

			event_.resize(storage_.GetChrNumber());
			std::cout << storage_.GetChrVerticesCount(0) << std::endl;
			time_t start = clock();
			starter_ = 0;
			tbb::parallel_for(tbb::blocked_range<size_t>(0, shuffle.size()), CheckIfSource(*this, shuffle));
			std::cout << "Time: " << time(0) - mark << std::endl;
			std::cout << source_.size() << std::endl;
			CheckSmart();
		}

		void CheckSmart()
		{
			for (size_t chr = 0; chr < event_.size(); chr++)
			{
				std::sort(event_[chr].begin(), event_[chr].end());
				for (int64_t i = 1; i < event_[chr].size(); i++)
				{
					if (!event_[chr][i].fork.branch[0].IsPositiveStrand())
					{
						std::cerr << i;
					}

					if (!event_[chr][i].isSource)
					{
						int64_t bestJ = i;
						int64_t minDiff = INT_MAX;
						int64_t minLength = INT_MAX;
						auto & it = event_[chr][i].fork.branch[1];
						for (int64_t j = i - 1; j >= 0; j--)
						{
							auto & jt = event_[chr][j].fork.branch[1];
							if (it.GetChrId() == jt.GetChrId() && it.IsPositiveStrand() == jt.IsPositiveStrand() &&
								((it.IsPositiveStrand() && jt.GetPosition() < it.GetPosition()) || (!it.IsPositiveStrand() && jt.GetPosition() > it.GetPosition())))
							{
								if (!event_[chr][j].isSource)
								{
									continue;
								}

								int64_t length[2];
								for (size_t l = 0; l < 2; l++)
								{
									length[l] = abs(event_[chr][i].fork.branch[l].GetPosition() - event_[chr][j].fork.branch[l].GetPosition());
								}

								auto diff = abs(length[0] - length[1]);
								if (diff < minDiff && diff < 1000)
								{
									minLength = Min(length[0], length[1]);
									minDiff = diff;
									bestJ = j;
								}
							}
						}

						if (bestJ < i && minLength > minBlockSize_)
						{
							std::vector<int64_t> length;
							int64_t currentBlock = ++blocksFound_;
							for (size_t l = 0; l < 2; l++)
							{
								auto it = event_[chr][bestJ].fork.branch[l];
								auto jt = event_[chr][i].fork.branch[l];
								if (jt.IsPositiveStrand())
								{
									auto start = it.GetPosition();
									auto end = jt.GetPosition() + k_;
									length.push_back(end - start);									
									blocksInstance_.push_back(BlockInstance(+currentBlock, jt.GetChrId(), it.GetPosition(), jt.GetPosition() + k_));
								}
								else
								{
									auto start = jt.GetPosition();
									auto end = it.GetPosition() - k_;
								//	if (start > end)
								//		std::cerr << l << ' ' << start << ' ' << end << std::endl;
									blocksInstance_.push_back(BlockInstance(-currentBlock, jt.GetChrId(), jt.GetPosition() - k_, it.GetPosition()));
								}
							}

							//std::cerr << length[0] - length[1] << std::endl;
						}
					}
				}
			}
		}

		struct ChrSweep
		{
		public:
			BlocksFinder & finder;

			ChrSweep(BlocksFinder & finder) : finder(finder)
			{

			}

			void operator()(tbb::blocked_range<size_t> & range) const
			{
				for (size_t r = range.begin(); r != range.end(); r++)
				{
					Sweeper sweeper(finder.storage_.Begin(r));
					sweeper.Sweep(finder.minBlockSize_, finder.maxBranchSize_, finder.k_, finder.blocksFound_, finder.template_, finder.globalMutex_);
					{
						tbb::mutex::scoped_lock lock(finder.globalMutex_);
						std::cout << ++finder.progressCount_ << ' ' << finder.storage_.GetChrNumber() << std::endl;
					}
				}
			}
		};

		struct RunTemplate
		{
		public:
			BlocksFinder & finder;

			RunTemplate(BlocksFinder & finder) : finder(finder)
			{

			}

			void operator()(tbb::blocked_range<size_t> & range) const
			{
				std::vector<size_t> data;
				std::vector<uint32_t> count(finder.storage_.GetVerticesNumber() * 2 + 1, 0);
				Path currentPath(finder.storage_, finder.maxBranchSize_, finder.minBlockSize_, finder.minBlockSize_, finder.maxFlankingSize_, true);
				for (size_t i = range.begin(); i < range.end(); i++)
				{
					//std::cout << i << ' ' << range.end() << std::endl;
					auto run = finder.template_[i];
					while (run.first < run.second)
					{
						for (; run.first.IsUsed() && run.first < run.second; ++run.first);
						auto end = run.first;
						for (; !end.IsUsed() && end < run.second; ++end);
						if (end.GetPosition() - run.first.GetPosition() >= finder.minBlockSize_)
						{
							auto init = run.first;
							auto length = end.GetPosition() - run.first.GetPosition();
							currentPath.Init(run.first.GetVertexId(), 'N');
							for (; run.first != end; ++run.first)
							{
								currentPath.PointPushBack(run.first.OutgoingEdge(), (run.first.GetPosition() - init.GetPosition()) <= finder.maxBranchSize_);
							}
							
							//currentPath.DumpInstances(std::cerr);
							if (currentPath.GoodInstancesList().size() > 1)
							{
								int64_t currentBlock = ++finder.blocksFound_;
							//	std::cerr << "g " << length << std::endl;
								for (auto it : currentPath.GoodInstancesList())
								{
									auto & jt = *it;
									if (jt.Front().IsPositiveStrand())
									{
										finder.blocksInstance_.push_back(BlockInstance(+currentBlock, jt.Front().GetChrId(), jt.Front().GetPosition(), jt.Back().GetPosition() + finder.k_));
									}
									else
									{
										finder.blocksInstance_.push_back(BlockInstance(-currentBlock, jt.Front().GetChrId(), jt.Back().GetPosition() - finder.k_, jt.Front().GetPosition()));
									}

							//		std::cerr << jt.RealLength() << std::endl;
									for (auto it = jt.Front(); it != jt.Back(); ++it)
									{
										it.MarkUsed();
									}
								}
							}

							currentPath.Clear();
						}

						run.first = end;
					}
				}
			}
		};
		
		

		void ListBlocksSequences(const BlockList & block, const std::string & directory) const
		{
			std::vector<IndexPair> group;
			BlockList blockList = block;
			GroupBy(blockList, compareById, std::back_inserter(group));
			for (std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
			{
				std::ofstream out;
				std::stringstream ss;
				ss << directory << "/" << blockList[it->first].GetBlockId() << ".fa";
				TryOpenFile(ss.str(), out);
				for (size_t block = it->first; block < it->second; block++)
				{
					int64_t length = blockList[block].GetLength();
					size_t chr = blockList[block].GetChrId();
					int64_t chrSize = storage_.GetChrSequence(chr).size();
					out << ">" << blockList[block].GetBlockId() << "_" << block - it->first << " ";
					out << storage_.GetChrDescription(chr) << ";";
					if (blockList[block].GetSignedBlockId() > 0)
					{
						out << blockList[block].GetStart() << ";" << length << ";" << "+;" << chrSize << std::endl;
						OutputLines(storage_.GetChrSequence(chr).begin() + blockList[block].GetStart(), length, out);
					}
					else
					{
						int64_t start = chrSize - blockList[block].GetEnd();
						out << start << ";" << length << ";" << "-;" << chrSize << std::endl;
						std::string::const_reverse_iterator it(storage_.GetChrSequence(chr).begin() + blockList[block].GetEnd());
						OutputLines(CFancyIterator(it, TwoPaCo::DnaChar::ReverseChar, ' '), length, out);
					}

					out << std::endl;
				}

				//std::cout << std::endl;
			}
		}

		struct SortByMultiplicity
		{
			SortByMultiplicity(const std::vector<int> & multiplicityOrigin) : multiplicity(multiplicityOrigin)
			{
			}

			bool operator () (const BlockInstance & a, const BlockInstance & b) const
			{
				auto mlp1 = multiplicity[a.GetBlockId()];
				auto mlp2 = multiplicity[b.GetBlockId()];
				if (mlp1 != mlp2)
				{
					return mlp1 > mlp2;
				}

				return a.GetBlockId() < b.GetBlockId();
			}

			const std::vector<int> & multiplicity;
		};


		void GenerateOutput(const std::string & outDir, bool genSeq)
		{
			const auto & trimmedBlocks = blocksInstance_;
			std::vector<std::vector<bool> > covered(storage_.GetChrNumber());
			for (size_t i = 0; i < covered.size(); i++)
			{
				covered[i].assign(storage_.GetChrSequence(i).size() + 1, false);
			}

			for (auto & b : blocksInstance_)
			{
				for (size_t i = b.GetStart(); i < b.GetEnd(); i++)
				{
					covered[b.GetChrId()][i] = true;
				}
			}

			size_t total = 0;
			size_t totalCovered = 0;
			for (auto & chr : covered)
			{
				total += chr.size();
				totalCovered += std::count(chr.begin(), chr.end(), true);
			}

			std::cout.setf(std::cout.fixed);
			std::cout.precision(2);
			std::cout << "Blocks found: " << blocksFound_ << std::endl;
			std::cout << "Coverage: " << double(totalCovered) / total << std::endl;

			CreateOutDirectory(outDir);
			std::string blocksDir = outDir + "/blocks";
			ListBlocksIndicesGFF(trimmedBlocks, outDir + "/" + "blocks_coords.gff");
			if (genSeq)
			{
				CreateOutDirectory(blocksDir);
				ListBlocksSequences(trimmedBlocks, blocksDir);
			}

			ListBlocksIndices(trimmedBlocks, outDir + "/" + "blocks_coords.txt");

		}

		struct BranchData
		{
			std::vector<size_t> branchId;
		};

		typedef std::vector<std::vector<size_t> > BubbledBranches;

		struct Fork
		{
			Fork(JunctionStorage::JunctionSequentialIterator it, JunctionStorage::JunctionSequentialIterator jt)
			{
				if (it < jt)
				{
					branch[0] = it;
					branch[1] = jt;
				}
				else
				{
					branch[0] = jt;
					branch[1] = it;
				}
			}

			bool operator == (const Fork & f) const
			{
				return branch[0] == f.branch[0] && branch[1] == f.branch[1];
			}

			bool operator != (const Fork & f) const
			{
				return !(*this == f);
			}

			std::string ToString() const
			{
				std::stringstream ss;
				for (size_t l = 0; l < 2; l++)
				{
					ss << (branch[l].IsPositiveStrand() ? '+' : '-') << ' ' << branch[l].GetChrId() << ' ' << branch[l].GetPosition() << ' ' << branch[l].GetChar() << ' ' << branch[l].GetVertexId() << "; ";
				}

				return ss.str();
			}

			bool operator < (const Fork & f) const
			{
				//return branch[0] < f.branch[0];
				return std::make_pair(branch[0], branch[1]) < std::make_pair(f.branch[0], f.branch[1]);
			}

			JunctionStorage::JunctionSequentialIterator branch[2];
		};

		struct Event
		{
			bool isSource;
			Fork fork;

			Event(bool isSource, Fork fork) : isSource(isSource), fork(fork)
			{

			}

			bool operator < (const Event & e) const
			{
				return fork.branch[0] < e.fork.branch[0];
			}
		};

		int64_t ChainLength(const Fork & now, const Fork & next) const
		{
			return min(abs(now.branch[0].GetPosition() - next.branch[0].GetPosition()), abs(now.branch[1].GetPosition() - next.branch[1].GetPosition()));
		}

		void BubbledBranchesForward(int64_t vertexId, const std::vector<JunctionStorage::JunctionSequentialIterator> & instance, BubbledBranches & bulges) const
		{
			std::vector<size_t> parallelEdge[5];
			std::map<int64_t, BranchData> visit;
			bulges.assign(instance.size(), std::vector<size_t>());
			for (size_t i = 0; i < instance.size(); i++)
			{
				auto vertex = instance[i];
				if ((vertex + 1).Valid())
				{
					parallelEdge[TwoPaCo::DnaChar::MakeUpChar(vertex.GetChar())].push_back(i);
				}

				for (int64_t startPosition = vertex++.GetPosition(); vertex.Valid() && abs(startPosition - vertex.GetPosition()) <= maxBranchSize_; ++vertex)
				{
					int64_t nowVertexId = vertex.GetVertexId();
					auto point = visit.find(nowVertexId);
					if (point == visit.end())
					{
						BranchData bData;
						bData.branchId.push_back(i);
						visit[nowVertexId] = bData;
					}
					else
					{
						point->second.branchId.push_back(i);
					}
				}
			}

			for (size_t i = 0; i < 5; i++)
			{
				for (size_t j = 0; j < parallelEdge[i].size(); j++)
				{
					for (size_t k = j + 1; k < parallelEdge[i].size(); k++)
					{
						size_t smallBranch = parallelEdge[i][j];
						size_t largeBranch = parallelEdge[i][k];
						bulges[smallBranch].push_back(largeBranch);
					}
				}
			}

			for (auto point = visit.begin(); point != visit.end(); ++point)
			{
				std::sort(point->second.branchId.begin(), point->second.branchId.end());
				for (size_t j = 0; j < point->second.branchId.size(); j++)
				{
					for (size_t k = j + 1; k < point->second.branchId.size(); k++)
					{
						size_t smallBranch = point->second.branchId[j];
						size_t largeBranch = point->second.branchId[k];
						if (smallBranch != largeBranch && std::find(bulges[smallBranch].begin(), bulges[smallBranch].end(), largeBranch) == bulges[smallBranch].end())
						{
							bulges[smallBranch].push_back(largeBranch);
						}
					}
				}
			}
		}

		void BubbledBranchesBackward(int64_t vertexId, const std::vector<JunctionStorage::JunctionSequentialIterator> & instance, BubbledBranches & bulges) const
		{
			std::vector<size_t> parallelEdge[5];
			std::map<int64_t, BranchData> visit;
			bulges.assign(instance.size(), std::vector<size_t>());
			for (size_t i = 0; i < instance.size(); i++)
			{
				auto iPrev = instance[i] - 1;
				if (iPrev.Valid())
				{
					for (size_t j = i + 1; j < instance.size(); j++)
					{
						auto jPrev = instance[j] - 1;
						if (jPrev.Valid() && iPrev.GetVertexId() == jPrev.GetVertexId() && iPrev.GetChar() == jPrev.GetChar())
						{
							bulges[i].push_back(j);
						}
					}
				}
			}


			for (size_t i = 0; i < instance.size(); i++)
			{
				auto vertex = instance[i];
				auto prev = vertex - 1;

				for (int64_t startPosition = vertex--.GetPosition(); vertex.Valid() && abs(startPosition - vertex.GetPosition()) <= maxBranchSize_; --vertex)
				{
					int64_t nowVertexId = vertex.GetVertexId();
					auto point = visit.find(nowVertexId);
					if (point == visit.end())
					{
						BranchData bData;
						bData.branchId.push_back(i);
						visit[nowVertexId] = bData;
					}
					else
					{
						point->second.branchId.push_back(i);
					}
				}
			}

			for (auto point = visit.begin(); point != visit.end(); ++point)
			{
				std::sort(point->second.branchId.begin(), point->second.branchId.end());
				for (size_t j = 0; j < point->second.branchId.size(); j++)
				{
					for (size_t k = j + 1; k < point->second.branchId.size(); k++)
					{
						size_t smallBranch = point->second.branchId[j];
						size_t largeBranch = point->second.branchId[k];
						if (smallBranch != largeBranch && std::find(bulges[smallBranch].begin(), bulges[smallBranch].end(), largeBranch) == bulges[smallBranch].end())
						{
							bulges[smallBranch].push_back(largeBranch);
						}
					}
				}
			}
		}

		template<class T> T Min(const T & a, const T & b)
		{
			if (a < b)
			{
				return a;
			}

			return b;
		}

		struct CheckIfSource
		{
		public:
			BlocksFinder & finder;
			std::vector<int64_t> & shuffle;

			CheckIfSource(BlocksFinder & finder, std::vector<int64_t> & shuffle) : finder(finder), shuffle(shuffle)
			{
			}

			void operator()(tbb::blocked_range<size_t> & range) const
			{
				BubbledBranches forwardBubble;
				BubbledBranches backwardBubble;
				std::vector<JunctionStorage::JunctionSequentialIterator> instance;
				for (size_t r = range.begin(); r != range.end(); r++)
				{
					if (finder.count_++ % 10000 == 0)
					{
						tbb::mutex::scoped_lock lock(finder.globalMutex_);
						std::cout << finder.count_ << std::endl;
					}

					instance.clear();
					int64_t vertex = shuffle[r];
					for (auto it = JunctionStorage::JunctionIterator(vertex); it.Valid(); ++it)
					{
						instance.push_back(it.SequentialIterator());
					}

					finder.BubbledBranchesForward(vertex, instance, forwardBubble);
					finder.BubbledBranchesBackward(vertex, instance, backwardBubble);
					for (size_t i = 0; i < forwardBubble.size(); i++)
					{
						for (size_t j = 0; j < forwardBubble[i].size(); j++)
						{
							size_t k = forwardBubble[i][j];
							if (!instance[i].IsPositiveStrand() && !instance[k].IsPositiveStrand())
							{
								continue;
							}

							auto it = std::find(backwardBubble[i].begin(), backwardBubble[i].end(), k);
							if (it == backwardBubble[i].end() && (instance[i].IsPositiveStrand()))
							{
								tbb::mutex::scoped_lock lock(finder.globalMutex_);
								int64_t minChrId = finder.Min(instance[i].GetChrId(), instance[k].GetChrId());
								//finder.source_[minChrId].push_back(Fork(instance[i], instance[k]));
								finder.event_[minChrId].push_back(Event(true, Fork(instance[i], instance[k])));
							}
						}
					}
				
					for (size_t i = 0; i < backwardBubble.size(); i++)
					{
						for (size_t j = 0; j < backwardBubble[i].size(); j++)
						{
							size_t k = backwardBubble[i][j];
							if (std::find(forwardBubble[i].begin(), forwardBubble[i].end(), k) == forwardBubble[i].end() && (instance[i].IsPositiveStrand()))
							{
								tbb::mutex::scoped_lock lock(finder.globalMutex_);
								int64_t minChrId = finder.Min(instance[i].GetChrId(), instance[k].GetChrId());
								//finder.sink_[minChrId].push_back(Fork(instance[i], instance[k]));
								finder.event_[minChrId].push_back(Event(false, Fork(instance[i], instance[k])));
							}
						}
					}
				}
			}
		};


		template<class T>
		static void AddIfNotExists(std::vector<T> & adj, T value)
		{
			if (std::find(adj.begin(), adj.end(), value) == adj.end())
			{
				adj.push_back(value);
			}
		}

		struct NextVertex
		{
			int64_t diff;
			int64_t count;
			JunctionStorage::JunctionSequentialIterator origin;
			NextVertex() : count(0)
			{

			}

			NextVertex(int64_t diff, JunctionStorage::JunctionSequentialIterator origin) : origin(origin), diff(diff), count(1)
			{

			}
		};


	private:

		template<class Iterator>
		void OutputLines(Iterator start, size_t length, std::ostream & out) const
		{
			for (size_t i = 1; i <= length; i++, ++start)
			{
				out << *start;
				if (i % 80 == 0 && i != length)
				{
					out << std::endl;
				}
			}
		}

		void ListChrs(std::ostream & out) const;
		std::string OutputIndex(const BlockInstance & block) const;
		void OutputBlocks(const std::vector<BlockInstance>& block, std::ofstream& out) const;
		void ListBlocksIndices(const BlockList & block, const std::string & fileName) const;
		void ListBlocksIndicesGFF(const BlockList & blockList, const std::string & fileName) const;
		void TryOpenFile(const std::string & fileName, std::ofstream & stream) const;


		std::vector<std::vector<Fork> > sink_;
		std::vector<std::vector<Fork> > source_;
		std::vector<std::vector<Event> > event_;
		std::vector<size_t> origin_;
		std::vector<std::pair<Fork, Fork> > matchSource_;


		int64_t k_;
		size_t progressCount_;
		size_t progressPortion_;
		std::atomic<int64_t> count_;
		std::atomic<int64_t> starter_;
		std::atomic<int64_t> blocksFound_;
		std::vector<size_t> pointComponent_;

		bool scoreFullChains_;
		int64_t scalingFactor_;
		int64_t minBlockSize_;
		int64_t maxBranchSize_;
		int64_t maxFlankingSize_;
		JunctionStorage & storage_;
		tbb::mutex progressMutex_;
		tbb::mutex globalMutex_;
		std::ofstream debugOut_;
		std::vector<Template> template_;
		std::vector<BlockInstance> blocksInstance_;

		//std::ofstream forkLog;


#ifdef _DEBUG_OUT_
		bool debug_;
		std::set<int64_t> missingVertex_;
#endif
	};
}

#endif
