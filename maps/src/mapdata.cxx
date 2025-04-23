#include "mapdata.h"

template <typename T>
typename SparseMapData<T>::const_iterator
SparseMapData<T>::const_iterator::operator++() {
	const_iterator end(sparse_.end());

	if (x > end.x || sparse_.data_.size() == 0) {
		x = end.x;
		y = end.y;
		return *this;
	}

	if (x < sparse_.offset_) {
		x = sparse_.offset_;
		y = sparse_.data_[0].first;
		return *this;
	}

	const SparseMapData<T>::data_element &column = sparse_.data_[x - sparse_.offset_];
	if (column.second.size() > 0 && y < (size_t) column.first) {
		y = column.first;
		return *this;
	}
	if (column.second.size() > 0 && y < (size_t) column.first + column.second.size() - 1) {
		y++;
		return *this;
	}

	do {
		x++;
		if (x > end.x) {
			x = end.x;
			y = end.y;
			return *this;
		}
	} while (sparse_.data_[x - sparse_.offset_].second.size() == 0);

	y = sparse_.data_[x - sparse_.offset_].first;

	return *this;
}

template <typename T>
DenseMapData *
SparseMapData<T>::to_dense() const
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

template <typename T>
SparseMapData<T>::SparseMapData(const DenseMapData &dense_map) :
    xlen_(dense_map.xdim()), ylen_(dense_map.ydim()), offset_(0)
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

template <typename T>
void
SparseMapData<T>::compact()
{
	if (data_.size() == 0)
		return;

	for (size_t ix = 0; ix < data_.size(); ix++) {
		data_element &column = data_[ix];
		if (column.second.size() == 0)
			continue;
		while (column.second.size() > 0 && column.second[column.second.size() - 1] == 0)
			column.second.pop_back();
		while (column.second.size() > 0 && column.second[0] == 0) {
			column.second.erase(column.second.begin());
			column.first++;
		}
		if (column.second.size() == 0)
			column.first = 0;
	}

	while (data_.size() > 0 && data_[data_.size() - 1].second.size() == 0)
		data_.pop_back();
	while (data_.size() > 0 && data_[0].second.size() == 0) {
		data_.erase(data_.begin());
		offset_++;
	}

	if (data_.size() == 0)
		offset_ = 0;
}


template <>
SparseMapData<double> &
SparseMapData<double>::operator+=(const SparseMapData<double> &r)
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

template <>
SparseMapData<double> &
SparseMapData<double>::operator+=(const DenseMapData &r)
{
	assert(xlen_ == r.xdim());
	assert(ylen_ == r.ydim());

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			double val = r.at(x, y);
			if (val != 0)
				(*this)(x, y) += val;
		}
	}

	return *this;
}

template <>
SparseMapData<double> &
SparseMapData<double>::operator-=(const SparseMapData<double> &r)
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

template <>
SparseMapData<double> &
SparseMapData<double>::operator-=(const DenseMapData &r)
{
	assert(xlen_ == r.xdim());
	assert(ylen_ == r.ydim());

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			double val = r.at(x, y);
			if (val != 0)
				(*this)(x, y) -= val;
		}
	}

	return *this;
}

template <>
SparseMapData<double> &
SparseMapData<double>::operator*=(const SparseMapData<double> &r)
{
	assert(xlen_ == r.xlen_);
	assert(ylen_ == r.ylen_);

	for (size_t ix = 0; ix < data_.size(); ix++) {
		size_t x = offset_ + ix;
		data_element &column = data_[ix];
		for (size_t iy = 0; iy < column.second.size(); iy++) {
			size_t y = column.first + iy;
			column.second[iy] *= r.at(x, y);
		}
	}

	return *this;
}

template <>
SparseMapData<double> &
SparseMapData<double>::operator*=(const DenseMapData &r)
{
	assert(xlen_ == r.xdim());
	assert(ylen_ == r.ydim());

	for (size_t ix = 0; ix < data_.size(); ix++) {
		size_t x = offset_ + ix;
		data_element &column = data_[ix];
		for (size_t iy = 0; iy < column.second.size(); iy++) {
			size_t y = column.first + iy;
			column.second[iy] *= r.at(x, y);
		}
	}

	return *this;
}

template <>
SparseMapData<double> &
SparseMapData<double>::operator*=(double r)
{
	for (size_t ix = 0; ix < data_.size(); ix++) {
		data_element &column = data_[ix];
		for (size_t iy = 0; iy < column.second.size(); iy++) {
			column.second[iy] *= r;
		}
	}

	return *this;
}

template <>
SparseMapData<double> &
SparseMapData<double>::operator/=(const SparseMapData<double> &r)
{
	assert(xlen_ == r.xlen_);
	assert(ylen_ == r.ylen_);

	// Division by zero is not the same as doing nothing.
	// Have to do this the long and painful way
	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			double val = this->at(x, y);
			double valr = r.at(x, y);
			if (valr == 0 || valr != valr || val != 0)
				(*this)(x, y) /= valr;
		}
	}

	return *this;
}

template <>
SparseMapData<double> &
SparseMapData<double>::operator/=(const DenseMapData &r)
{
	assert(xlen_ == r.xdim());
	assert(ylen_ == r.ydim());

	// Division by zero is not the same as doing nothing.
	// Have to do this the long and painful way
	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			double val = this->at(x, y);
			double valr = r.at(x, y);
			if (valr == 0 || valr != valr || val != 0)
				(*this)(x, y) /= valr;
		}
	}

	return *this;
}

template <>
SparseMapData<double> &
SparseMapData<double>::operator/=(double r)
{
	assert(r != 0);

	for (size_t ix = 0; ix < data_.size(); ix++) {
		data_element &column = data_[ix];
		for (size_t iy = 0; iy < column.second.size(); iy++) {
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
			(*this)(x, y) += r.at(x, y);
		}
	}

	return *this;
}

DenseMapData &
DenseMapData::operator+=(const SparseMapData<double> &r)
{
	assert(xlen_ == r.xdim());
	assert(ylen_ == r.ydim());

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			(*this)(x, y) += r.at(x, y);
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
			(*this)(x, y) -= r.at(x, y);
		}
	}

	return *this;
}

DenseMapData &
DenseMapData::operator-=(const SparseMapData<double> &r)
{
	assert(xlen_ == r.xdim());
	assert(ylen_ == r.ydim());

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			(*this)(x, y) -= r.at(x, y);
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
			(*this)(x, y) *= r.at(x, y);
		}
	}

	return *this;
}

DenseMapData &
DenseMapData::operator*=(const SparseMapData<double> &r)
{
	assert(xlen_ == r.xdim());
	assert(ylen_ == r.ydim());

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			(*this)(x, y) *= r.at(x, y);
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
			(*this)(x, y) /= r.at(x, y);
		}
	}

	return *this;
}

DenseMapData &
DenseMapData::operator/=(const SparseMapData<double> &r)
{
	assert(xlen_ == r.xdim());
	assert(ylen_ == r.ydim());

	for (size_t x = 0; x < xlen_; x++) {
		for (size_t y = 0; y < ylen_; y++) {
			(*this)(x, y) /= r.at(x, y);
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

template class SparseMapData<double>;
template class SparseMapData<bool>;

