#include <vector>

#include <cereal/types/vector.hpp>

class SparseMapData;
class DenseMapData;

class SparseMapData {
public:
	SparseMapData(size_t xlen, size_t ylen) :
	    xlen_(xlen), ylen_(ylen), offset_(0) {}
	SparseMapData(const DenseMapData &dense_map);

	bool in_bounds(size_t x, size_t y) const {
		return !(x < 0 || x >= xlen_ || y < 0 || y >= ylen_);
	}

	double set(size_t x, size_t y, double val) {
		assert(x >= 0);
		assert(x < xlen_);
		assert(y >= 0);
		assert(y < ylen_);

		if (x < offset_) {
			data_.insert(data_.begin(), offset_-x, data_element());
			offset_ = x;
		} else if (x >= offset_ + data_.size()) {
			data_.resize(x - offset_ + 1);
		}
		data_element &column = data_[x-offset_];

		if (y < column.first) {
			column.second.insert(column.second.begin(),
			    y-column.first, double(0));
			column.first = y;
		} else if (y >= column.first + column.second.size()) {
			column.second.resize(y - column.first + 1, double(0));
		}
		column.second[y - column.first] = val;
		return val;
	}

	double get(size_t x, size_t y) const {
		if (x < offset_ || x >= offset_ + data_.size())
			return 0;
		const data_element &column = data_[x-offset_];
		if (y < column.first ||
		    y >= column.first + column.second.size())
			return 0;
		return column.second[y-column.first];
	}

	size_t xdim() const { return xlen_; }
	size_t ydim() const { return ylen_; }

	// dynamically deallocate the sparse map while filling?
	DenseMapData *to_dense() const;

	template <class A> void serialize(A &ar, unsigned v) {
		using namespace cereal;
		// XXX: size_t vs. uint64_t
		ar & make_nvp("xlen", xlen_);
		ar & make_nvp("xlen", ylen_);
		ar & make_nvp("offset", offset_);
		ar & make_nvp("data", data_);
	}

private:
	size_t xlen_, ylen_;
	typedef std::pair<int, std::vector<double> > data_element;
	std::vector<data_element> data_;
	size_t offset_;
};


class DenseMapData {
public:
	DenseMapData(size_t xlen, size_t ylen) :
	    xlen_(xlen), ylen_(ylen) {
		data_.resize(xlen*ylen);
	}

	bool in_bounds(size_t x, size_t y) const {
		return !(x < 0 || x >= xlen_ || y < 0 || y >= ylen_);
	}

	size_t xdim() const { return xlen_; }
	size_t ydim() const { return ylen_; }

	double get(size_t x, size_t y) const {
		if (!in_bounds(x, y))
			return 0;
		return data_[idxat(x, y)];
	}

	double set(size_t x, size_t y, double val) {
		assert(x >= 0);
		assert(x < xlen_);
		assert(y >= 0);
		assert(y < ylen_);

		data_[idxat(x, y)] = val;
		return val;
	}

	template <class A> void serialize(A &ar, unsigned v) {
		using namespace cereal;
		// XXX: size_t vs. uint64_t
		ar & make_nvp("xlen", xlen_);
		ar & make_nvp("xlen", ylen_);
		ar & make_nvp("data", data_);
	}

private:
	size_t xlen_, ylen_;
	std::vector<double> data_;

	inline size_t idxat(size_t x, size_t y) const {
		return x + y * xlen_;
	}
};

DenseMapData *
SparseMapData::to_dense() const
{
	DenseMapData *rv = new DenseMapData(xlen_, ylen_);
	for (size_t ix = 0; ix < data_.size(); ix++) {
		size_t x = offset_ + ix;
		const data_element &column = data_[ix];
		for (size_t iy = 0; iy < column.second.size(); iy++) {
			size_t y = column.first + iy;
			rv->set(x, y, column.second[iy]);
		}
	}
	return rv;
}

SparseMapData::SparseMapData(const DenseMapData &dense_map) :
    xlen_(dense_map.xdim()), ylen_(dense_map.ydim())
{
	double val;

	for (size_t ix = 0; ix < xlen_; ix++) {
		for (size_t iy = 0; iy < ylen_; iy++) {
			val = dense_map.get(ix, iy);
			if (val != 0)
				set(ix, iy, val);
		}
	}
}

