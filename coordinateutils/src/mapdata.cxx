#include <boost/make_shared.hpp>
#include <serialization.h>

#include "mapdata.h"

DenseMapData *
SparseMapData::to_dense() const
{
	DenseMapData *rv = new DenseMapData(xlen_, ylen_);
	for (size_t ix = 0; ix < data_.size(); ix++) {
		size_t x = offset_ + ix;
		const data_element &column = data_[ix];
		for (size_t iy = 0; iy < column.second.size(); iy++) {
			size_t y = column.first + iy;
			(*rv)(x, y) = column.second[iy];
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
			val = dense_map(ix, iy);
			if (val != 0)
				(*this)(ix, iy) = val;
		}
	}
}


SparseMapData &
SparseMapData::operator+=(const SparseMapData &r)
{
	assert(xlen_ == r.xlen_);
	assert(ylen_ == r.ylen_);

	for (size_t ix = 0; ix < r.data_.size(); ix++) {
		size_t x = r.offset_ + ix;
		const data_element &column = r.data_[ix];
		for (size_t iy = 0; iy < column.second.size(); iy++) {
			size_t y = column.first + iy;
			double val = column.second[iy];
			if (val != 0)
				(*this)(x, y) += val;
		}
	}

	return *this;
}

SparseMapData &
SparseMapData::operator+=(const DenseMapData &r)
{
	assert(xlen_ == r.xdim());
	assert(ylen_ == r.ydim());

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			double val = r(x, y);
			if (val != 0)
				(*this)(x, y) += val;
		}
	}

	return *this;
}

SparseMapData &
SparseMapData::operator-=(const SparseMapData &r)
{
	assert(xlen_ == r.xlen_);
	assert(ylen_ == r.ylen_);

	for (size_t ix = 0; ix < r.data_.size(); ix++) {
		size_t x = r.offset_ + ix;
		const data_element &column = r.data_[ix];
		for (size_t iy = 0; iy < column.second.size(); iy++) {
			size_t y = column.first + iy;
			double val = column.second[iy];
			if (val != 0)
				(*this)(x, y) -= val;
		}
	}

	return *this;
}

SparseMapData &
SparseMapData::operator-=(const DenseMapData &r)
{
	assert(xlen_ == r.xdim());
	assert(ylen_ == r.ydim());

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			double val = r(x, y);
			if (val != 0)
				(*this)(x, y) -= val;
		}
	}

	return *this;
}

SparseMapData &
SparseMapData::operator*=(const SparseMapData &r)
{
	assert(xlen_ == r.xlen_);
	assert(ylen_ == r.ylen_);

	for (size_t ix = 0; ix < data_.size(); ix++) {
		size_t x = offset_ + ix;
		data_element &column = data_[ix];
		for (size_t iy = 0; iy < column.second.size(); iy++) {
			size_t y = column.first + iy;
			column.second[iy] *= r(x, y);
		}
	}

	return *this;
}

SparseMapData &
SparseMapData::operator*=(const DenseMapData &r)
{
	assert(xlen_ == r.xdim());
	assert(ylen_ == r.ydim());

	for (size_t ix = 0; ix < data_.size(); ix++) {
		size_t x = offset_ + ix;
		data_element &column = data_[ix];
		for (size_t iy = 0; iy < column.second.size(); iy++) {
			size_t y = column.first + iy;
			column.second[iy] *= r(x, y);
		}
	}

	return *this;
}

SparseMapData &
SparseMapData::operator*=(double r)
{
	for (size_t ix = 0; ix < data_.size(); ix++) {
		size_t x = offset_ + ix;
		data_element &column = data_[ix];
		for (size_t iy = 0; iy < column.second.size(); iy++) {
			size_t y = column.first + iy;
			column.second[iy] *= r;
		}
	}

	return *this;
}

SparseMapData &
SparseMapData::operator/=(const SparseMapData &r)
{
	assert(xlen_ == r.xlen_);
	assert(ylen_ == r.ylen_);

	for (size_t ix = 0; ix < data_.size(); ix++) {
		size_t x = offset_ + ix;
		data_element &column = data_[ix];
		for (size_t iy = 0; iy < column.second.size(); iy++) {
			size_t y = column.first + iy;
			column.second[iy] /= r(x, y);
		}
	}

	return *this;
}

SparseMapData &
SparseMapData::operator/=(const DenseMapData &r)
{
	assert(xlen_ == r.xdim());
	assert(ylen_ == r.ydim());

	for (size_t ix = 0; ix < data_.size(); ix++) {
		size_t x = offset_ + ix;
		data_element &column = data_[ix];
		for (size_t iy = 0; iy < column.second.size(); iy++) {
			size_t y = column.first + iy;
			column.second[iy] /= r(x, y);
		}
	}

	return *this;
}

SparseMapData &
SparseMapData::operator/=(double r)
{
	for (size_t ix = 0; ix < data_.size(); ix++) {
		size_t x = offset_ + ix;
		data_element &column = data_[ix];
		for (size_t iy = 0; iy < column.second.size(); iy++) {
			size_t y = column.first + iy;
			column.second[iy] /= r;
		}
	}

	return *this;
}

DenseMapData &
DenseMapData::operator+=(const DenseMapData &r)
{
	assert(xlen_ == r.xlen_);
	assert(ylen_ == r.ylen_);

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			(*this)(x, y) += r(x, y);
		}
	}

	return *this;
}

DenseMapData &
DenseMapData::operator+=(const SparseMapData &r)
{
	assert(xlen_ == r.xdim());
	assert(ylen_ == r.ydim());

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			(*this)(x, y) += r(x, y);
		}
	}

	return *this;
}

DenseMapData &
DenseMapData::operator+=(double r)
{
	if (r == 0)
		return *this;

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			(*this)(x, y) += r;
		}
	}

	return *this;
}

DenseMapData &
DenseMapData::operator-=(const DenseMapData &r)
{
	assert(xlen_ == r.xlen_);
	assert(ylen_ == r.ylen_);

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			(*this)(x, y) -= r(x, y);
		}
	}

	return *this;
}

DenseMapData &
DenseMapData::operator-=(const SparseMapData &r)
{
	assert(xlen_ == r.xdim());
	assert(ylen_ == r.ydim());

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			(*this)(x, y) -= r(x, y);
		}
	}

	return *this;
}

DenseMapData &
DenseMapData::operator-=(double r)
{
	if (r == 0)
		return *this;

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			(*this)(x, y) -= r;
		}
	}

	return *this;
}

DenseMapData &
DenseMapData::operator*=(const DenseMapData &r)
{
	assert(xlen_ == r.xlen_);
	assert(ylen_ == r.ylen_);

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			(*this)(x, y) *= r(x, y);
		}
	}

	return *this;
}

DenseMapData &
DenseMapData::operator*=(const SparseMapData &r)
{
	assert(xlen_ == r.xdim());
	assert(ylen_ == r.ydim());

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			(*this)(x, y) *= r(x, y);
		}
	}

	return *this;
}

DenseMapData &
DenseMapData::operator*=(double r)
{
	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			(*this)(x, y) *= r;
		}
	}

	return *this;
}

DenseMapData &
DenseMapData::operator/=(const DenseMapData &r)
{
	assert(xlen_ == r.xlen_);
	assert(ylen_ == r.ylen_);

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			(*this)(x, y) /= r(x, y);
		}
	}

	return *this;
}

DenseMapData &
DenseMapData::operator/=(const SparseMapData &r)
{
	assert(xlen_ == r.xdim());
	assert(ylen_ == r.ydim());

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			(*this)(x, y) /= r(x, y);
		}
	}

	return *this;
}

DenseMapData &
DenseMapData::operator/=(double r)
{
	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			(*this)(x, y) /= r;
		}
	}

	return *this;
}

