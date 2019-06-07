#include <vector>
#include <map>
#include <stdexcept>
#include <boost/python.hpp>


class MapData {
public:
	MapData(size_t xlen, size_t ylen) :
	    xlen_(xlen), ylen_(ylen) {}

	size_t xdim() const { return xlen_; }
	size_t ydim() const { return ylen_; }

	virtual double & operator() (size_t x, size_t y);
	virtual double operator() (size_t x, size_t y) const;

private:
	bool check_index(size_t x, size_t y, bool raise) {
		bool valid = !(x < 0 || x >= xlen_ || y < 0 || y >= ylen_);
		if (!valid && raise) {
			std::ostringstream os;
			os << "Map index (" << x << ", " << y << ") out of range";
			throw std::out_of_range(os.str());
		}
		return valid;
	}

	size_t xlen_, ylen_;
};

class FilledMapData;

class SparseMapData : public MapData {
public:
	SparseMapData(size_t xlen, size_t ylen) :
	    MapData(xlen, ylen), offset(0) {}

	double & operator() (size_t x, size_t y) {
		(void) check_index(x, y, true);

		if (x < offset) {
			data.insert(data.begin(), offset - x, data_element());
			offset = x;
		} else if (x >= offset + data.size()) {
			data.resize(x - offset + 1);
		}
		data_element &column = data[x - offset];

		if (y < column.first) {
			column.second.insert(column.second.begin(),
			    y - column.first, double(0));
			column.first = y;
		} else if (y >= column.first + column.second.size()) {
			column.second.resize(y - column.first + 1, double(0));
		}
		return column.second[y - column.first];
	}

	double operator() (size_t x, size_t y) const {
		if (x < offset || x >= offset + data.size())
			return 0;
		const data_element &column = data[x - offset];

		if (y < column.first ||
		    y >= column.first + column.second.size())
			return 0;
		return column.second[y - column.first];
	}

	// dynamically deallocate the sparse map while filling?
	boost::shared_ptr<FilledMapData> fill() const {
		FilledMapData arr(xlen_, ylen_);
		for (size_t ix = 0; ix < data.size(); ix++) {
			size_t x = offset + ix;
			data_element &column = data[ix];
			for (size_t iy = 0; iy < column.second.size(); iy++) {
				size_t y = column.first + iy;
				arr(x, y) = column.second[iy];
			}
		}
		return boost::make_shared<FilledMapData>(arr);
	}

private:
	typedef std::pair<int, std::vector<double> > data_element;
	std::vector<data_element> data;
	size_t offset;
};


class FilledMapData : public MapData {
public:
	FilledMapData(size_t xlen, size_t ylen) :
	    MapData(xlen, ylen) {}

	double & operator() (size_t x, size_t y) {
		(void) check_index(x, y, true);
		allocate();
		return data[idxat(x, y)];
	}

	double operator() (size_t x, size_t y) const {
		if (!check_index(x, y, false))
			return 0;
		if (!allocated())
			return 0;
		return data[idxat(x, y)];
	}

	// dynamically deallocate the filled map while compressing?
	boost::shared_ptr<SparseMapData> sparsen() const {
		SparseMapData arr(xlen_, ylen_);
		for (size_t x = 0; x < xlen_; x++) {
			for (size_t y = 0; y < ylen_; y++) {
				double val = data[idxat(x, y)];
				if (val != 0)
					arr(x, y) = val;
			}
		}
		return boost::make_shared<SparseMapData>(arr);
	}

private:
	std::vector<double> data;

	bool allocated() const {
		return (data.size() > 0);
	}

	size_t allocate() {
		if (!allocated())
			data.resize(xlen_ * ylen_);
		return data.size();
	}

	inline size_t idxat(size_t x, size_t y) const {
		return x + y * xlen_;
	}
};
