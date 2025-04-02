#ifndef _G3_CONTAINERS_H
#define _G3_CONTAINERS_H

#include <vector>
#include <unordered_map>
#include <iterator>

template <typename K, typename T>
class OrderedMap {
public:
	using map_type = std::unordered_map<K, T>;
	using key_map_type = std::vector<K>;
	using key_type = K;
	using mapped_type = T;
	using value_type = typename map_type::value_type;
	using difference_type = typename map_type::difference_type;
	using size_type = typename map_type::size_type;

	OrderedMap() {};
	OrderedMap(const OrderedMap& other) :
	    map_(other.map_), keys_(other.keys_) {}
	OrderedMap(OrderedMap&& other) :
	    map_(std::move(other.map_)), keys_(std::move(other.keys_)) {}
	virtual ~OrderedMap() {};

	OrderedMap& operator=(const OrderedMap&) = default;
	OrderedMap& operator=(OrderedMap&&) = default;

	mapped_type& operator[](const key_type& key) {
		auto it = std::find(keys_.begin(), keys_.end(), key);
		if (it == keys_.end())
			keys_.push_back(key);
		return map_[key];
	}
	const mapped_type& operator[](const key_type& key) const { return map_.at(key); }
	mapped_type& at(const key_type& key) { return map_.at(key); }
	const mapped_type& at(const key_type& key) const { return map_.at(key); }
	bool contains(const key_type& key) const { return map_.count(key) > 0; }

	size_type size() const { return map_.size(); }
	size_type max_size() const { return map_.max_size(); }
	bool empty() const { return map_.empty(); }

	void clear() {
		map_.clear();
		keys_.clear();
	}

	size_type erase(const key_type& key) {
		auto it = std::find(keys_.begin(), keys_.end(), key);
		if (it == keys_.end())
			return 0;

		keys_.erase(it);
		return map_.erase(key);
	}

	class iterator {
	public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = typename map_type::value_type;
		using difference_type = typename map_type::difference_type;
		using pointer = value_type*;
		using reference = value_type&;

		typename key_map_type::iterator cur;
		map_type& map;

		iterator(typename key_map_type::iterator it, map_type& map) :
		    cur(it), map(map) {}

		iterator& operator++() {
			++cur;
			return *this;
		}

		iterator operator++(int) {
			iterator tmp = *this;
			++(*this);
			return tmp;
		}

		bool operator==(const iterator& other) const { return cur == other.cur; }
		bool operator!=(const iterator& other) const { return cur != other.cur; }

		pointer operator->() const { return &(*(map.find(*cur))); }
		reference operator*() const { return *(map.find(*cur)); }
	};

	class const_iterator {
	public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = typename map_type::value_type;
		using difference_type = typename map_type::difference_type;
		using pointer = const value_type*;
		using reference = const value_type&;

		typename key_map_type::const_iterator cur;
		const map_type& map;

		const_iterator(typename key_map_type::const_iterator it,
		    const map_type& map) : cur(it), map(map) {}

		const_iterator(const iterator& other) : cur(other.cur), map(other.map) {}

		const_iterator& operator++() {
			++cur;
			return *this;
		}

		const_iterator operator++(int) {
			const_iterator tmp = *this;
			++(*this);
			return tmp;
		}

		bool operator==(const const_iterator& other) const { return cur == other.cur; }
		bool operator!=(const const_iterator& other) const { return cur != other.cur; }

		pointer operator->() const { return &(*(map.find(*cur))); }
		reference operator*() const { return *(map.find(*cur)); }
	};

	iterator begin() { return iterator(keys_.begin(), map_); }
	iterator end() { return iterator(keys_.end(), map_); }
	const_iterator begin() const { return const_iterator(keys_.begin(), map_); }
	const_iterator end() const { return const_iterator(keys_.end(), map_); }
	const_iterator cbegin() const { return const_iterator(keys_.cbegin(), map_); }
	const_iterator cend() const { return const_iterator(keys_.cend(), map_); }

	iterator find(const key_type& key) {
		return iterator(std::find(keys_.begin(), keys_.end(), key), map_);
	}

	const_iterator find(const key_type& key) const {
		return const_iterator(std::find(keys_.cbegin(), keys_.cend(), key), map_);
	}

	std::pair<iterator, bool> insert(const value_type& pair) {
		const key_type& key = pair.first;

		bool check = false;
		auto it = std::find(keys_.begin(), keys_.end(), key);
		if (it == keys_.end()) {
			keys_.push_back(key);
			it = std::prev(keys_.end());
			check = true;
		}
		map_.insert(pair);
		return {iterator(it, map_), check};
	}

	std::pair<iterator, bool> insert(value_type&& pair) {
		const key_type& key = pair.first;

		bool check = false;
		auto it = std::find(keys_.begin(), keys_.end(), key);
		if (it == keys_.end()) {
			keys_.push_back(key);
			it = std::prev(keys_.end());
			check = true;
		}
		map_.insert(pair);
		return {iterator(it, map_), check};
	}

protected:
	map_type map_;
	key_map_type keys_;
};

#endif
