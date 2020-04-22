#ifndef _JUNCTION_STORAGE_H_
#define _JUNCTION_STORAGE_H_

#include <set>
#include <atomic>
#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include <stdexcept>
#include <algorithm>

#include <streamfastaparser.h>
#include <junctionapi.h>

namespace Sibelia
{
	using std::min;
	using std::max;

	
	class JunctionStorage
	{
	private:

		struct Position
		{
			uint32_t pos;
			uint32_t nextIdx;
			int32_t nextChr;
			int64_t vertexId;
			uint16_t pointerIdx;
			char ch;
			char revCh;
			bool invert;

			Position(const TwoPaCo::JunctionPosition & junction) : nextIdx(UINT_MAX), vertexId(junction.GetId())
			{
				pos = junction.GetPos();
			}
		};

		struct PrevPosition
		{
			PrevPosition() : prevId(0)
			{

			}

			int64_t prevId;
			uint32_t prevChr;
			uint32_t prevIdx;
		};

		typedef std::vector<Position> PositionVector;

	public:

		class Iterator
		{
		public:
			Iterator(): chrId_(SIZE_MAX)
			{

			}

			Iterator(size_t chrId) : chrId_(chrId), idx_(0), isPositive_(true)
			{

			}

			Iterator(size_t chrId, int32_t idx, bool isPositive = true) : chrId_(chrId), idx_(idx), isPositive_(isPositive)
			{

			}

			void Inc()
			{
				idx_ += 1;
			}

			void DecInSequence()
			{
				if (IsPositiveStrand())
				{
					idx_ -= 1;
				}
				else
				{
					idx_ += 1;
				}
			}

			void Next()
			{
				auto & pos = JunctionStorage::this_->position_[chrId_][idx_];
				if (pos.nextIdx != UINT_MAX)
				{
					if (pos.invert)
					{
						isPositive_ = !isPositive_;
					}

					idx_ = pos.nextIdx;
					chrId_ = pos.nextChr;
				}
				else
				{
					chrId_ = SIZE_MAX;
				}

			}

			int32_t GetPointerIndex() const
			{
				return JunctionStorage::this_->position_[chrId_][idx_].pointerIdx;
			}

			int32_t GetChrId() const
			{
				return chrId_;
			}

			int32_t PreviousPosition() const
			{
				if (IsPositiveStrand())
				{
					return JunctionStorage::this_->position_[GetChrId()][idx_ - 1].pos;
				}

				return -(JunctionStorage::this_->position_[GetChrId()][idx_ + 1].pos + JunctionStorage::this_->k_);
			}

			int32_t GetVertexId() const
			{
				if (IsPositiveStrand())
				{
					return JunctionStorage::this_->position_[GetChrId()][idx_].vertexId;
				}

				return -JunctionStorage::this_->position_[GetChrId()][idx_].vertexId;
			}

			int32_t GetPosition() const
			{
				if (IsPositiveStrand())
				{
					return JunctionStorage::this_->position_[GetChrId()][idx_].pos;
				}

				return -(JunctionStorage::this_->position_[GetChrId()][idx_].pos + JunctionStorage::this_->k_);
			}

			char GetChar() const
			{
				if (IsPositiveStrand())
				{
					return JunctionStorage::this_->position_[chrId_][idx_].ch;
				}

				return JunctionStorage::this_->position_[chrId_][idx_].revCh;
			}

			bool IsPositiveStrand() const
			{
				return isPositive_;
			}

			size_t GetIndex() const
			{
				return idx_;
			}

			bool Valid() const
			{
				return chrId_ < JunctionStorage::this_->position_.size() && idx_ < JunctionStorage::this_->position_[chrId_].size();
			}

			bool operator == (const Iterator & arg) const
			{
				return this->chrId_ == arg.chrId_ && this->idx_ == arg.idx_ && isPositive_ == arg.isPositive_;
			}

			bool operator != (const Iterator & arg) const
			{
				return !(*this == arg);
			}


		private:
			size_t chrId_;
			size_t idx_;
			bool isPositive_;
		};

		int64_t GetVertexId(size_t chr, size_t idx) const
		{
			return position_[chr][idx].vertexId;
		}

		int64_t GetMaxVertexId() const
		{
			return maxId_;
		}

		int64_t GetPosition(size_t chr, size_t idx) const
		{
			return position_[chr][idx].pos;
		}


		int32_t GetPointerIndex(size_t chr, size_t idx) const
		{
			return position_[chr][idx].pointerIdx;
		}

		int64_t GetChrNumber() const
		{
			return position_.size();
		}

		const std::string& GetChrDescription(uint64_t idx) const
		{
			return sequenceDescription_[idx];
		}

		size_t GeChrSequenceSize(size_t chr) const
		{
			return chrSeqSize_[chr];
		}

		size_t GeChrSize(size_t chr) const
		{
			return position_[chr].size();
		}

		void Init(const std::string & inFileName, const std::vector<std::string> & genomesFileName, int64_t threads, int64_t abundanceThreshold, int64_t loopThreshold)
		{
			this_ = this;
			size_t record = 0;
			maxId_ = 0;

			std::vector<std::string> sequence_(position_.size());
			for (const auto & fastaFileName : genomesFileName)
			{
				for (TwoPaCo::StreamFastaParser parser(fastaFileName); parser.ReadRecord(); record++)
				{
					sequenceDescription_.push_back(parser.GetCurrentHeader());
					sequenceId_[parser.GetCurrentHeader()] = sequenceDescription_.size() - 1;
					sequence_.push_back(std::string());
					for (char ch; parser.GetChar(ch); )
					{
						sequence_[record].push_back(ch);
					}

					chrSeqSize_.push_back(sequence_[record].size());
				}
			}

			std::vector<PrevPosition> prevPos;
			position_.resize(sequence_.size());
			TwoPaCo::JunctionPositionReader reader(inFileName);
			for (TwoPaCo::JunctionPosition junction; reader.NextJunctionPosition(junction);)
			{
				size_t absId = abs(junction.GetId());
				maxId_ = max(absId, maxId_);
				if (absId >= prevPos.size())
				{
					prevPos.resize(absId + 1);
				}
			
				{
					auto chr = junction.GetChr();
					auto pos = junction.GetPos();
					position_[chr].push_back(Position(junction));
					position_[chr].back().ch = sequence_[chr][pos + JunctionStorage::this_->k_];
					position_[chr].back().revCh = pos > 0 ? TwoPaCo::DnaChar::ReverseChar(sequence_[chr][pos - 1]) : 'N';
					if (prevPos[absId].prevId != 0)
					{
						auto & prev = prevPos[absId];
						position_[prev.prevChr][prev.prevIdx].nextChr = junction.GetChr();
						position_[prev.prevChr][prev.prevIdx].invert = prev.prevId != junction.GetId();
						position_[prev.prevChr][prev.prevIdx].nextIdx = position_[junction.GetChr()].size() - 1;
						position_[chr].back().pointerIdx = position_[prev.prevChr][prev.prevIdx].pointerIdx + 1;
					}
					else
					{
						position_[chr].back().pointerIdx = 0;
					}

					prevPos[absId].prevId = junction.GetId();
					prevPos[absId].prevChr = junction.GetChr();
					prevPos[absId].prevIdx = position_[junction.GetChr()].size() - 1;

		
				}
			}
		}

		struct Pointer
		{
			int32_t chrId;
			int32_t idx;

			Pointer() {}
			Pointer(int32_t chrId, int32_t idx) : chrId(chrId), idx(idx)
			{

			}

			bool operator < (const Pointer & p) const
			{
				if (abs(chrId) != abs(p.chrId))
				{
					return abs(chrId) < abs(p.chrId);
				}

				return idx < p.idx;
			}
		};

		
		JunctionStorage() {}
		JunctionStorage(const std::string & fileName, const std::vector<std::string> & genomesFileName, uint64_t k, int64_t threads, int64_t abundanceThreshold, int64_t loopThreshold) : k_(k), abundance_(abundanceThreshold)
		{
			Init(fileName, genomesFileName, threads, abundanceThreshold, loopThreshold);
		}

		size_t GetAbundance() const
		{
			return abundance_;
		}

	private:


		int64_t k_;
		size_t maxId_;
		size_t abundance_;
		std::map<std::string, size_t> sequenceId_;
		std::vector<size_t> chrSeqSize_;
		std::vector<std::string> sequenceDescription_;
		std::vector<std::vector<Position> > position_;
		static JunctionStorage * this_;
		friend class Iterator;
	};
}

#endif
