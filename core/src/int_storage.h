#ifndef _G3_INT_STORAGE_H
#define _G3_INT_STORAGE_H

#include <vector>
#include <string>
#include <map>

#include <cereal/cereal.hpp>

int bit_count(const std::vector<int64_t> &d);
int bit_count(const std::map<std::string, int64_t> &d);

template <class A, typename FROM_TYPE>
void load_as(A &ar, std::vector<int64_t> &dest) {
	std::vector<FROM_TYPE> temp;
	ar & cereal::make_nvp("vector", temp);
	dest.resize(temp.size());
	std::copy(temp.begin(), temp.end(), dest.begin());
}

template <class A, typename FROM_TYPE>
void load_as(A &ar, std::map<std::string, int64_t> &dest) {
	std::map<std::string, FROM_TYPE> temp;
	ar & cereal::make_nvp("map", temp);
	dest.insert(temp.begin(), temp.end());
}

template <class A, typename TO_TYPE>
void save_as(A &ar, const std::vector<int64_t> &src) {
	std::vector<TO_TYPE> temp(src.begin(), src.end());
	ar & cereal::make_nvp("vector", temp);
}

template <class A, typename TO_TYPE>
void save_as(A &ar, const std::map<std::string, int64_t> &src) {
	std::map<std::string, TO_TYPE> temp(src.begin(), src.end());
	ar & cereal::make_nvp("map", temp);
}

#endif
