#include <vector>

#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp>

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

	size_t xdim() const { return xlen_; }
	size_t ydim() const { return ylen_; }
	size_t nonzero() const {
		size_t nz = 0;
		for (size_t i = 0; i < data_.size(); i++)
			nz += data_[i].second.size();
		return nz;
	}

	double at(size_t x, size_t y) const {
		if (x < offset_ || x >= offset_ + data_.size())
			return 0;
		const data_element &column = data_[x-offset_];
		if (y < column.first ||
		    y >= column.first + column.second.size())
			return 0;
		return column.second[y-column.first];
	}

	double operator()(size_t x, size_t y) const { return at(x, y); }

	double &operator()(size_t x, size_t y) {
		assert(x >= 0);
		assert(x < xlen_);
		assert(y >= 0);
		assert(y < ylen_);

		if (data_.size() == 0) {
			data_.resize(1);
			offset_ = x;
		} else if (x < offset_) {
			data_.insert(data_.begin(), offset_-x, data_element());
			offset_ = x;
		} else if (x >= offset_ + data_.size()) {
			data_.resize(x - offset_ + 1);
		}
		data_element &column = data_[x-offset_];

		if (column.second.size() == 0) {
			column.first = y;
			column.second.resize(1, double(0));
		} else if (y < column.first) {
			column.second.insert(column.second.begin(),
			    column.first-y, double(0));
			column.first = y;
		} else if (y >= column.first + column.second.size()) {
			column.second.resize(y - column.first + 1, double(0));
		}
		return column.second[y - column.first];
	}

	SparseMapData &operator+=(const SparseMapData &r);
	SparseMapData &operator-=(const SparseMapData &r);
	SparseMapData &operator*=(const SparseMapData &r);
	SparseMapData &operator/=(const SparseMapData &r);

	SparseMapData &operator+=(const DenseMapData &r);
	SparseMapData &operator-=(const DenseMapData &r);
	SparseMapData &operator*=(const DenseMapData &r);
	SparseMapData &operator/=(const DenseMapData &r);

	SparseMapData &operator*=(double r);
	SparseMapData &operator/=(double r);

	SparseMapData *clone(bool copy_data = true) {
		SparseMapData *m = new SparseMapData(xlen_, ylen_);
		if (copy_data) {
			m->data_ = data_;
			m->offset_ = offset_;
		}
		return m;
	}

	DenseMapData *to_dense() const;

	template <class A> void serialize(A &ar, unsigned v) {
		using namespace cereal;
		ar & make_nvp("xlen", xlen_);
		ar & make_nvp("ylen", ylen_);
		ar & make_nvp("offset", offset_);
		ar & make_nvp("data", data_);
	}

	class const_iterator {
	public:
		const_iterator(const SparseMapData &sparse, size_t x_, size_t y_) :
		    x(x_), y(y_), sparse_(sparse) {}

		size_t x, y;

		bool operator==(const const_iterator & other) const {
			return ((x == other.x) && (y == other.y));
		}
		bool operator!=(const const_iterator & other) const {
			return ((x != other.x) || (y != other.y));
		}

		double operator*() const { return sparse_.at(x, y); }

		const_iterator operator++();
		const_iterator operator++(int) { const_iterator i = *this; ++(*this); return i; }

	private:
		const SparseMapData & sparse_;
	};

	const_iterator begin() const {
		if (data_.size() == 0)
			return const_iterator(*this, 0, 0);
		return const_iterator(*this, offset_, data_[0].first);
	}

	const_iterator end() const {
		if (data_.size() == 0)
			return const_iterator(*this, 0, 0);
		size_t x = offset_ + data_.size() - 1;
		const data_element &column = data_[x - offset_];
		size_t y = column.first + column.second.size();
		return const_iterator(*this, x, y);
	}

private:
	uint64_t xlen_, ylen_;
	typedef std::pair<int, std::vector<double> > data_element;
	std::vector<data_element> data_;
	uint64_t offset_;
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

	double operator()(size_t x, size_t y) const {
		if (!in_bounds(x, y))
			return 0;
		return data_[idxat(x, y)];
	}

	double &operator()(size_t x, size_t y) {
		assert(x >= 0);
		assert(x < xlen_);
		assert(y >= 0);
		assert(y < ylen_);

		return data_[idxat(x, y)];
	}

	DenseMapData &operator=(const DenseMapData &r)
	{
		xlen_ = r.xlen_;
		ylen_ = r.ylen_;
		data_ = r.data_;
		return *this;
	}
	DenseMapData &operator=(const std::vector<double> &d)
	{
		assert(data_.size() == d.size());
		data_ = d;
		return *this;
	}

	DenseMapData &operator+=(const DenseMapData &r);
	DenseMapData &operator-=(const DenseMapData &r);
	DenseMapData &operator*=(const DenseMapData &r);
	DenseMapData &operator/=(const DenseMapData &r);

	DenseMapData &operator+=(const SparseMapData &r);
	DenseMapData &operator-=(const SparseMapData &r);
	DenseMapData &operator*=(const SparseMapData &r);
	DenseMapData &operator/=(const SparseMapData &r);

	DenseMapData &operator+=(double r);
	DenseMapData &operator-=(double r);
	DenseMapData &operator*=(double r);
	DenseMapData &operator/=(double r);

	DenseMapData *clone(bool copy_data = true) {
		DenseMapData *m = new DenseMapData(xlen_, ylen_);
		if (copy_data) {
			m->data_ = data_;
		}
		return m;
	}

	template <class A> void serialize(A &ar, unsigned v) {
		using namespace cereal;
		ar & make_nvp("xlen", xlen_);
		ar & make_nvp("ylen", ylen_);
		ar & make_nvp("data", data_);
	}

	class const_iterator {
	public:
		const_iterator(const DenseMapData &dense, size_t x_, size_t y_) :
		    x(x_), y(y_), dense_(dense) {}

		size_t x, y;

		bool operator==(const const_iterator & other) const {
			return ((x == other.x) && (y == other.y));
		}
		bool operator!=(const const_iterator & other) const {
			return ((x != other.x) || (y != other.y));
		}

		double operator*() const {
			return dense_.data_[dense_.idxat(x, y)];
		}

		const_iterator operator++() {
			size_t idx = dense_.idxat(x, y);
			if (idx < dense_.data_.size()) {
				++idx;
				x = idx % dense_.xlen_;
				y = idx / dense_.xlen_;
			}
			return *this;
		}

		const_iterator operator++(int) { const_iterator i = *this; ++(*this); return i; }

	private:
		const DenseMapData & dense_;
	};

	const_iterator begin() const {
		return const_iterator(*this, 0, 0);
	}

	const_iterator end() const {
		return const_iterator(*this, 0, ylen_);
	}

private:
	uint64_t xlen_, ylen_;
	std::vector<double> data_;

	inline size_t idxat(size_t x, size_t y) const {
		return x + y * xlen_;
	}
};

CEREAL_CLASS_VERSION(SparseMapData, 1);
CEREAL_CLASS_VERSION(DenseMapData, 1);
