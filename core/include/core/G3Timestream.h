#ifndef _G3_TIMESTREAM_H
#define _G3_TIMESTREAM_H

#include <G3Frame.h>
#include <G3TimeStamp.h>
#include <G3Vector.h>
#include <vector>
#include <map>

class G3Timestream : public G3FrameObject {
public:
	enum TimestreamUnits {
		None = 0,
		Counts = 1,
		Current = 2,
		Power = 3,
		Resistance = 5,
		Tcmb = 4,
		Angle = 6,
		Distance = 7,
		Voltage = 8,
		Pressure = 9,
		FluxDensity = 10,
		Trj = 11,
		Frequency = 12,
	};

	G3Timestream(const G3Timestream &r);
	G3Timestream(std::vector<double>::size_type s = 0, double val = 0) :
	    units(None), use_flac_(0),
	    buffer_((s == 0) ? NULL : new std::vector<double>(s, val)),
	    data_((s == 0) ? NULL : &(*buffer_)[0]), len_(s),
	    data_type_(TS_DOUBLE) {}
	template <typename Iterator> G3Timestream(Iterator l, Iterator r) :
	    units(None), use_flac_(0), buffer_(new std::vector<double>(l, r)),
	    data_(&(*buffer_)[0]), len_(buffer_->size()), data_type_(TS_DOUBLE) {}
	virtual ~G3Timestream() {
		if (buffer_) delete buffer_;
	}

	// FLAC compression levels range from 0-9. 0 means do not use FLAC.
	void SetFLACCompression(int compression_level);

	TimestreamUnits units;
	G3Time start, stop;

	// Accessors
	typedef double value_type;
	double &operator[](size_t i) {
		if (data_type_ != TS_DOUBLE)
			throw std::runtime_error("Cannot access non-double "
			    "timestream read/write");
		return ((double *)data_)[i];
	}
	double operator[](size_t i) const noexcept {
		switch (data_type_) {
		case TS_DOUBLE:
			return ((double *)data_)[i];
		case TS_FLOAT:
			return ((float *)data_)[i];
		case TS_INT32:
			return ((int32_t *)data_)[i];
		case TS_INT64:
			return ((int64_t *)data_)[i];
		}
		__builtin_unreachable();
	}
	double &at(size_t i) {
		if (i >= len_)
			throw std::out_of_range("Out of range");
		return (*this)[i];
	}
	double at(size_t i) const {
		if (i >= len_)
			throw std::out_of_range("Out of range");
		return (*this)[i];
	}
	size_t size() const noexcept {
		return len_;
	}
	double *begin() {
		if (data_type_ != TS_DOUBLE)
			throw std::runtime_error("Cannot iterate non-double "
			    "timestream");
		return (double *)data_;
	}
	double *end() {
		if (data_type_ != TS_DOUBLE)
			throw std::runtime_error("Cannot iterate non-double "
			    "timestream");
		return (double *)data_ + len_;
	}
	const double *begin() const {
		if (data_type_ != TS_DOUBLE)
			throw std::runtime_error("Cannot iterate non-double "
			    "timestream");
		return (const double *)data_;
	}
	const double *end() const {
		if (data_type_ != TS_DOUBLE)
			throw std::runtime_error("Cannot iterate non-double "
			    "timestream");
		return (const double *)data_ + len_;
	}

	double GetSampleRate() const;
	uint8_t GetCompressionLevel() const{ return use_flac_; }
	
	template <class A> void load(A &ar, unsigned v);
	template <class A> void save(A &ar, unsigned v) const;

	std::string Description() const;
	std::string Summary() const { return Description(); };

	// Arithmetic operations
	G3Timestream operator +(const G3Timestream &r) const;
	G3Timestream &operator +=(const G3Timestream &r);
	G3Timestream &operator -=(const G3Timestream &r);
	G3Timestream operator -(const G3Timestream &r) const;
	G3Timestream operator *(const G3Timestream &r) const;
	G3Timestream operator /(const G3Timestream &r) const;
	G3Timestream operator +(double r) const;
	G3Timestream operator -(double r) const;
	G3Timestream operator *(double r) const;
	G3Timestream &operator *=(double r);
	G3Timestream operator /(double r) const;
	G3Timestream &operator /=(double r);

	class G3TimestreamPythonHelpers;

private:
	SET_LOGGER("G3Timestream");

	friend class G3TimestreamMap;
	friend class G3TimestreamPythonHelpers;

	uint8_t use_flac_;

	std::vector<double> *buffer_;
	std::shared_ptr<void> root_data_ref_;
	void *data_;
	size_t len_;
	enum {
		TS_DOUBLE,
		TS_FLOAT,
		TS_INT32,
		TS_INT64
	} data_type_;
	
	template<typename T>
	struct TimeStreamTypeResolver{
		static_assert(sizeof(T)!=sizeof(T), "Unsupported datatype for G3Timestream");
	};
};

G3Timestream operator +(double l, const G3Timestream &r);
G3Timestream operator -(double l, const G3Timestream &r);
G3Timestream operator *(double l, const G3Timestream &r);
G3Timestream operator /(double l, const G3Timestream &r);

G3_POINTERS(G3Timestream);

template<>
struct G3Timestream::TimeStreamTypeResolver<double>{
	constexpr static auto type_tag=G3Timestream::TS_DOUBLE;
};

template<>
struct G3Timestream::TimeStreamTypeResolver<float>{
	constexpr static auto type_tag=G3Timestream::TS_FLOAT;
};

template<>
struct G3Timestream::TimeStreamTypeResolver<int32_t>{
	constexpr static auto type_tag=G3Timestream::TS_INT32;
};

template<>
struct G3Timestream::TimeStreamTypeResolver<int64_t>{
	constexpr static auto type_tag=G3Timestream::TS_INT64;
};

class G3TimestreamMap : public G3FrameObject,
    public std::map<std::string, G3TimestreamPtr> {
public:
	G3TimestreamMap(){}
	G3TimestreamMap(const G3TimestreamMap& other):std::map<std::string, G3TimestreamPtr>(other){}
	G3TimestreamMap(G3TimestreamMap&& other):std::map<std::string, G3TimestreamPtr>(std::move(other)){}

	G3TimestreamMap& operator=(const G3TimestreamMap&)=default;
	G3TimestreamMap& operator=(G3TimestreamMap&&)=default;

	// Return true if all timestreams start and end at the same time and
	// contain the same number of samples.
	bool CheckAlignment() const;

	// Return statistics for contained timestreams. Results undefined if
	// CheckAlignment is false.
	G3Time GetStartTime() const;
	void SetStartTime(G3Time start);
	G3Time GetStopTime() const;
	void SetStopTime(G3Time stop);
	double GetSampleRate() const;
	size_t NSamples() const;
	G3Timestream::TimestreamUnits GetUnits() const;
	void SetUnits(G3Timestream::TimestreamUnits units);
	/// FLAC compression levels range from 0-9. 0 means do not use FLAC.
	uint8_t GetCompressionLevel() const;
	void SetFLACCompression(int compression_level);

	// Compact underlying data storage into a contiguous 2D block.
	// This invalidates any references to data inside any member
	// timestreams (though not the timestream objects themselves)
	// and requires timestream alignment (throws exception if
	// CheckAlignment is false).
	void Compactify();
	
	/// Construct a map whose underlying data storage is a contiguous 2D block.
	/// \param keys the timestream keys for which the map should be constructed
	/// \param start the start time which will be shared by all time streams
	/// \param stop the stop time which will be shared by all time streams
	/// \param n_samples the number of samples each timestream will have
	/// \pre keys must be in sorted order
	template<typename SampleType>
	static G3TimestreamMap
	MakeCompact(const std::vector<std::string>& keys, std::size_t n_samples,
	    G3Time start, G3Time stop, G3Timestream::TimestreamUnits units=G3Timestream::None,
	    int compression_level=0) {
		std::shared_ptr<SampleType[]> data(new SampleType[n_samples*keys.size()]);
		return MakeCompact(keys, n_samples, data, start, stop, units, compression_level);
	}

	/// Construct a map using an existing contiguous 2D block of data as the underlying storage.
	/// Within the data block, the samples for each timestream must be laid out in the same order
	/// as the (sorted) keys, without gaps.
	/// \param keys the timestream keys for which the map should be constructed
	/// \param start the start time which will be shared by all time streams
	/// \param stop the stop time which will be shared by all time streams
	/// \param n_samples the number of samples each timestream will have
	/// \param data existing data into which the new timestreams should be views
	/// \pre keys must be in sorted order
	template<typename SampleType>
	static G3TimestreamMap
	MakeCompact(const std::vector<std::string>& keys, std::size_t n_samples,
	    std::shared_ptr<SampleType[]> data, G3Time start, G3Time stop,
	    G3Timestream::TimestreamUnits units=G3Timestream::None, int compression_level=0) {
		G3TimestreamMap map;
		map.FromBuffer(keys, n_samples, data, start, stop, units, compression_level);
		return map;
	}

	template<typename SampleType>
	void FromBuffer(const std::vector<std::string>& keys, std::size_t n_samples,
	    std::shared_ptr<SampleType[]> data, G3Time start, G3Time stop,
	    G3Timestream::TimestreamUnits units=G3Timestream::None, int compression_level=0) {
		const auto data_type=G3Timestream::TimeStreamTypeResolver<SampleType>::type_tag;
		std::size_t offset = 0;
		for (const auto& key : keys) {
			auto ts = std::make_shared<G3Timestream>(0);
			ts->start = start;
			ts->stop = stop;
			ts->units = units;
			ts->use_flac_ = compression_level;
			ts->root_data_ref_ = data;
			ts->data_ = data.get() + offset;
			ts->data_type_ = data_type;
			ts->len_ = n_samples;
			emplace(key, std::move(ts));
			offset += n_samples;
		}
	}

	template <class A> void serialize(A &ar, unsigned v);
	std::string Description() const;
};

G3_POINTERS(G3TimestreamMap);

namespace cereal {
	template <class A> struct specialize<A, G3Timestream, cereal::specialization::member_load_save> {};
	template <class A> struct specialize<A, G3TimestreamMap, cereal::specialization::member_serialize> {};
}

G3_SERIALIZABLE(G3Timestream, 3);
G3_SERIALIZABLE(G3TimestreamMap, 3);

#endif

