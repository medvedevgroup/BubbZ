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

		void FindBlocks(int32_t minBlockSize, int32_t maxBranchSize, int32_t threads, const std::string & debugOut)
		{
			blocksFound_ = 0;
			minBlockSize_ = minBlockSize;
			maxBranchSize_ = maxBranchSize;

			using namespace std::placeholders;

			count_ = 0;
			time_t start = clock();
			currentIndex_ = 0;
			workInstance_.resize(threads);

			std::cout << '[' << std::flush;
			progressPortion_ = storage_.GetChrNumber() / progressCount_;
			if (progressPortion_ == 0)
			{
				progressPortion_ = 1;
			}

			#pragma omp parallel num_threads(threads)
			{
				ChrSweep process(*this);
				process();
			}

			for (auto & outVector : workInstance_)
			{
				std::copy(outVector.begin(), outVector.end(), std::back_inserter(blocksInstance_));
				outVector.clear();
			}

			std::cout << ']' << std::endl;

			//std::cout << double(clock() - start) / CLOCKS_PER_SEC << std::endl;
		}

		struct ChrSweep
		{
		public:
			BlocksFinder & finder;

			ChrSweep(BlocksFinder & finder) : finder(finder)
			{
			}

			void operator()() const
			{
				std::vector<std::vector<InstanceSet > > instance(2, std::vector<InstanceSet>(finder.storage_.GetChrNumber()));
				for (size_t i = 0; i < 2; i++)
				{
					for (size_t j = 0; j < finder.storage_.GetChrNumber(); j++)
					{
						instance[i][j].Init(j, i == 0, finder.storage_.GeChrSize(j));
					}
				}

				std::vector<VertexEntry* > lastPosEntry_(finder.storage_.GetMaxVertexId() + 1, 0);
				std::vector<VertexEntry* > lastNegEntry_(finder.storage_.GetMaxVertexId() + 1, 0);

				size_t endIndex = finder.storage_.GetChrNumber();
				for(bool go = true; go;)
				{
					size_t nowChr;
					#pragma omp critical
					{
						if (finder.currentIndex_ < endIndex)
						{
							nowChr = finder.currentIndex_++;
						}
						else
						{
							go = false;
						}
					}

					if (go)
					{
						auto it = JunctionStorage::Iterator(nowChr);
						Sweeper sweeper(it, lastPosEntry_, lastNegEntry_);
						sweeper.Sweep(finder.storage_, finder.minBlockSize_, finder.maxBranchSize_, finder.k_, finder.blocksFound_, finder.workInstance_[omp_get_thread_num()], instance);
						{
							if (finder.count_++ % finder.progressPortion_ == 0)
							{
								std::cout << '.' << std::flush;
							}
						}
					}
					
				}
			}
		};
		

		void GenerateOutput(const std::string & outDir, bool genSeq, bool legacyOut)
		{
			const auto & trimmedBlocks = blocksInstance_;

			std::cout.setf(std::cout.fixed);
			std::cout.precision(2);
			std::cout << "Blocks found: " << blocksFound_ << std::endl;

			CreateOutDirectory(outDir);
			std::string blocksDir = outDir + "/blocks";
			ListBlocksIndicesGFF(trimmedBlocks, outDir + "/" + "blocks_coords.gff");

			if (legacyOut)
			{
				ListBlocksIndices(trimmedBlocks, outDir + "/" + "blocks_coords.txt");
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


		int64_t k_;
		size_t progressCount_;
		size_t progressPortion_;
		std::atomic<int64_t> count_;
		std::atomic<size_t> currentIndex_;
		std::atomic<int64_t> blocksFound_;

		int32_t minBlockSize_;
		int32_t maxBranchSize_;
		JunctionStorage & storage_;
		std::ofstream debugOut_;
		std::vector<BlockInstance> blocksInstance_;
		std::vector<std::vector<BlockInstance> > workInstance_;


		//std::ofstream forkLog;


#ifdef _DEBUG_OUT_
		bool debug_;
		std::set<int64_t> missingVertex_;
#endif
	};
}

#endif
