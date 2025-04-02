#ifndef G3_CONTAINERS_H
#define G3_CONTAINERS_H

#include <algorithm>
#include <list>
#include <unordered_map>

///This is an ordered associative container which supports unique keys.
///It supports bidirectional iterators.
///The ordering of entries in this container is unrelated to their keys;
///instead it is the order in which they were inserted.
template <typename K, typename T,
          typename Hash=std::hash<K>,
          typename Pred=std::equal_to<K>,
          typename DataAlloc = std::allocator<std::pair<const K, T>>,
          typename IndexAlloc = std::allocator<std::pair<const K,
                                    typename std::list<std::pair<const K, T>>::iterator>>>
class OrderedMap{
public:
	using key_type = K;
	using mapped_type = T;
	using value_type = std::pair<const key_type, mapped_type>;
	using hasher = Hash;
	using key_equal = Pred;
	using allocator_type = DataAlloc;
	using index_allocator_type = IndexAlloc;
	using reference = value_type&;
	using const_reference = const value_type&;
	using pointer = typename std::allocator_traits<allocator_type>::pointer;
	using const_pointer = typename std::allocator_traits<allocator_type>::const_pointer;
protected:
	using data_store_type = std::list<value_type, DataAlloc>;
public:
	using iterator = typename data_store_type::iterator;
	using const_iterator = typename data_store_type::const_iterator;
	using reverse_iterator = typename data_store_type::reverse_iterator;
	using const_reverse_iterator = typename data_store_type::const_reverse_iterator;
	//Not implemented: local iterators
	using size_type = typename data_store_type::size_type;
	using difference_type = typename data_store_type::difference_type;
protected:
	using index_type = std::unordered_map<K, iterator, Hash, Pred, IndexAlloc>;
public:
	
	OrderedMap(){};
	
	explicit OrderedMap(size_type n, const hasher& hf = hasher(),
	                    const key_equal& eql = key_equal(),
	                    const allocator_type& da = allocator_type(),
	                    const index_allocator_type& ia = index_allocator_type()):
	data_(n,hf,eql,da),index_(ia){}
	
	template<class InputIterator>
	OrderedMap(InputIterator f, InputIterator l, size_type n = 0,
	           const hasher& hf = hasher(),
	           const key_equal& eql = key_equal(),
	           const allocator_type& da = allocator_type(),
	           const index_allocator_type& ia = index_allocator_type()):
	data_(n,hf,eql,da),index_(ia){
		while(f!=l)
			insert(*f++);
	}
	
	template<class InputIterator>
	OrderedMap(std::initializer_list<value_type> il, size_type n = 0,
	           const hasher& hf = hasher(),
	           const key_equal& eql = key_equal(),
	           const allocator_type& da = allocator_type(),
	           const index_allocator_type& ia = index_allocator_type()):
	data_(n,hf,eql,da),index_(ia){
		for(value_type&& val : il)
			insert(std::move(val));
	}
	
	explicit OrderedMap(const allocator_type& da,
	                    const index_allocator_type& ia = index_allocator_type()):
	data_(da),index_(ia){}
	
	OrderedMap(const OrderedMap& other):
	data_(other.data_),index_(other.index_){}
	
	OrderedMap(OrderedMap&& other) :
	    data_(std::move(other.data_)),index_(std::move(other.index_)){}
	virtual ~OrderedMap(){};
	
	OrderedMap(const OrderedMap& other, const allocator_type& da,
	           const index_allocator_type& ia = index_allocator_type()):
	data_(other.data_, da),index_(other.index_, ia){}
	
	OrderedMap(size_type n, const allocator_type& da,
	           const index_allocator_type& ia = index_allocator_type()):
	OrderedMap(n, hasher(), key_equal(), da, ia){}
	
	OrderedMap(size_type n, const hasher& hf, const allocator_type& da,
	           const index_allocator_type& ia = index_allocator_type()):
	OrderedMap(n, hf, key_equal(), da, ia){}
	
	template<class InputIterator>
	OrderedMap(InputIterator f, InputIterator l, size_type n,
	           const allocator_type& da,
	           const index_allocator_type& ia = index_allocator_type()):
	OrderedMap(f, l, n, hasher(), key_equal(), da, ia){}
	
	template<class InputIterator>
	OrderedMap(InputIterator f, InputIterator l, size_type n,
	           const hasher& hf,
	           const allocator_type& da,
	           const index_allocator_type& ia = index_allocator_type()):
	OrderedMap(f, l, n, hf, key_equal(), da, ia){}
	
	OrderedMap(std::initializer_list<value_type> il, size_type n, const allocator_type& da,
	           const index_allocator_type& ia = index_allocator_type()):
	OrderedMap(il, n, hasher(), key_equal(), da, ia){}

	OrderedMap(std::initializer_list<value_type> il, size_type n, const hasher& hf,
	           const allocator_type& da,
	           const index_allocator_type& ia = index_allocator_type()):
	OrderedMap(il, n, hf, key_equal(), da, ia){}

	OrderedMap& operator=(const OrderedMap&) = default;
	OrderedMap& operator=(OrderedMap&&) = default;
	OrderedMap& operator=(std::initializer_list<value_type> il){
		clear();
		for(value_type&& val : il)
			insert(std::move(val));
		return *this;
	}
	
	allocator_type get_allocator() const noexcept{ return data_.get_allocator(); }
	index_allocator_type get_index_allocator() const noexcept{ return index_.get_allocator(); }
	
	iterator begin() noexcept{ return data_.begin(); }
	iterator end() noexcept{ return data_.end(); }
	const_iterator begin() const noexcept{ return data_.begin(); }
	const_iterator end() const noexcept{ return data_.end(); }
	const_iterator cbegin() const noexcept{ return data_.cbegin(); }
	const_iterator cend() const noexcept{ return data_.cend(); }
	reverse_iterator rbegin() noexcept{ return data_.rbegin(); }
	reverse_iterator rend() noexcept{ return data_.rend(); }
	const_reverse_iterator rbegin() const noexcept{ return data_.rbegin(); }
	const_reverse_iterator rend() const noexcept{ return data_.rend(); }
	const_reverse_iterator crbegin() const noexcept{ return data_.crbegin(); }
	const_reverse_iterator crend() const noexcept{ return data_.crend(); }
	
	bool empty() const noexcept{ return index_.empty(); }
	size_type size() const noexcept{ return index_.size(); }
	size_type max_size() const noexcept{ return data_.max_size(); }
	
	template <class... Args>
	std::pair<iterator,bool> emplace(Args&&... args){
		value_type val(std::forward<Args>(args)...);
		auto iit=index_.find(val.first);
		if(iit!=index_.end())
			return std::make_pair(iit->second, false);
		data_.emplace_back(std::move(val));
		auto it=data_.end();
		it--;
		index_.insert(make_pair(val.first, it));
		return std::make_pair(it, true);
	}
	
	template <class... Args>
	iterator emplace_hint(const_iterator, Args&&... args){
		return emplace(std::forward<Args>(args)...).first;
	}
	
	std::pair<iterator,bool> insert(const value_type& val){
		auto iit=index_.find(val.first);
		if(iit!=index_.end())
			return std::make_pair(iit->second, false);
		data_.push_back(val);
		auto it=data_.end();
		it--;
		index_.insert(make_pair(val.first, it));
		return std::make_pair(it, true);
	}
	
	std::pair<iterator,bool> insert(value_type&& val){
		auto iit=index_.find(val.first);
		if(iit!=index_.end())
			return std::make_pair(iit->second, false);
		data_.emplace_back(std::move(val));
		auto it=data_.end();
		it--;
		index_.insert(make_pair(val.first, it));
		return std::make_pair(it, true);
	}
	
	template<typename P>
	std::pair<iterator, bool> insert(P&& obj){
		return emplace(std::forward<P>(obj));
	}
	
	iterator insert(const_iterator hint, const value_type& obj){
		return emplace_hint(hint, obj);
	}
	
	iterator insert(const_iterator hint, value_type&& obj){
		return emplace_hint(hint, std::move(obj));
	}
	
	template<typename P>
	iterator insert(const_iterator hint, P&& obj){
		return emplace_hint(hint, std::forward<P>(obj));
	}
	
	template<class InputIterator>
	void insert(InputIterator first, InputIterator last){
		while(first!=last)
			insert(*first++);
	}
	
	void insert(std::initializer_list<value_type> il){
		for(auto&& v : il)
			insert(std::move(v));
	}
	
	//Not implemented: (c++17) extract
	//Not implemented: (c++17) insert(node_type)
	//Not implemented: (c++17) try_emplace
	//Not implemented: (c++17) insert_or_assign
	
	iterator erase(iterator position){
		const key_type& key=position->first;
		iterator result=position;
		result++;
		auto iit=index_.find(key);
		data_.erase(position);
		index_.erase(iit);
		return result;
	}
	
	iterator erase(const_iterator position){
		const key_type& key=position->first;
		iterator result=position;
		result++;
		auto iit=index_.find(key);
		data_.erase(position);
		index_.erase(iit);
		return result;
	}
	
	size_type erase(const key_type& key){
		auto iit=index_.find(key);
		if(iit==index_.end())
			return 0;
		data_.erase(iit->second);
		index_.erase(iit);
		return 1;
	}
	
	iterator erase(const_iterator first, const_iterator last){
		//figure out each affected index entry and remove them all
		for(auto it=first; it!=last; it++){
			const key_type& key=it->first;
			auto iit=index_.find(key);
			index_.erase(iit);
		}
		//erase main entries in bulk
		return data_.erase(first, last);
	}
	
	void swap(const OrderedMap& other){
		data_.swap(other.data_);
		index_.swap(other.index_);
	}
	
	//Not implemented: (c++17) merge
	
	hasher hash_function() const{ return index_.hash_function(); }
	key_equal key_eq() const{ return index_.key_eq(); }
	
	iterator find(const key_type& k){
		auto it=index_.find(k);
		if(it==index_.end())
			return end();
		return it->second;
	}
	
	const_iterator find(const key_type& k) const{
		auto it=index_.find(k);
		if(it==index_.end())
			return end();
		return it->second;
	}
	
	template<class K2>
	iterator find(const K2& k){
		auto it=index_.find(k);
		if(it==index_.end())
			return end();
		return it->second;
	}
	
	template<class K2>
	const_iterator find(const K2& k) const{
		auto it=index_.find(k);
		if(it==index_.end())
			return end();
		return it->second;
	}
	
	size_type count(const key_type& key) const{ return index_.count(key); }
	
	template<class K2>
	size_type count(const K2& key) const{ return index_.count(key); }
	
	bool contains(const key_type& key) const{ return index_.count(key); }
	
	template<class K2>
	bool contains(const K2& key) const{ return index_.count(key); }
	
	std::pair<iterator, iterator> equal_range(const key_type& key){
		iterator first=find(key);
		iterator last;
		if(first==end())
			last=end();
		else{
			last=first;
			last++;
		}
		return std::make_pair(first, last);
	}
	
	std::pair<const_iterator, const_iterator> equal_range(const key_type& key) const{
		const_iterator first=find(key);
		const_iterator last;
		if(first==end())
			last=end();
		else{
			last=first;
			last++;
		}
		return std::make_pair(first, last);
	}
	
	template<class K2>
	std::pair<iterator, iterator> equal_range(const K2& key){
		iterator first=find(key);
		iterator last;
		if(first==end())
			last=end();
		else{
			last=first;
			last++;
		}
		return std::make_pair(first, last);
	}
	
	template<class K2>
	std::pair<const_iterator, const_iterator> equal_range(const K2& key) const{
		const_iterator first=find(key);
		const_iterator last;
		if(first==end())
			last=end();
		else{
			last=first;
			last++;
		}
		return std::make_pair(first, last);
	}
	
	mapped_type& operator[](const key_type& key){
		auto it=index_.find(key);
		if(it==index_.end()){
			data_.push_back(std::make_pair(key,mapped_type{}));
			auto dit=data_.end();
			dit--;
			index_[key]=dit;
			return dit->second;
		}
		return it->second->second;
	}
	
	mapped_type& operator[](key_type&& key){
		auto it=index_.find(key);
		if(it==index_.end()){
			data_.push_back(std::make_pair(key,mapped_type{}));
			auto dit=data_.end();
			dit--;
			index_[key]=dit;
			return dit->second;
		}
		return it->second->second;
	}
	
	mapped_type& at(const key_type& k){
		auto it=index_.at(k);
		return it->second;
	}
	
	const mapped_type& at(const key_type& k) const{
		auto it=index_.at(k);
		return it->second;
	}
	
	size_type bucket_count() const noexcept{ return index_.bucket_count(); }
	size_type max_bucket_count() const noexcept{ return index_.max_bucket_count(); }
	size_type bucket_size(size_type n) const{ return index_.bucket_size(n); }
	size_type bucket(const key_type& k) const{ return index_.bucket(k); }
	
	float load_factor() const noexcept{ return index_.load_factor(); }
	float max_load_factor() const noexcept{ return index_.max_load_factor(); }
	void max_load_factor(float z){ index_.max_load_factor(z); }
	void rehash(size_type n){ index_.rehash(n); }
	void reserve(size_type n){ index_.reserve(n); /*no sensible way to reserve in data_*/}
	
	void clear(){
		data_.clear();
		index_.clear();
	}
	
protected:
	data_store_type data_;
	index_type index_;
};
	
//Not implemented: deduction guides

template<typename Key, typename T, typename Hash, typename Pred,
         typename DataAlloc, typename IndexAlloc>
bool operator==(const OrderedMap<Key, T, Hash, Pred, DataAlloc, IndexAlloc>& a,
                const OrderedMap<Key, T, Hash, Pred, DataAlloc, IndexAlloc>& b){
    if(a.size()!=b.size())
    	return false;
	return std::equal(a.cbegin(), a.cend(), b.cbegin());
}

template<typename Key, typename T, typename Hash, typename Pred,
         typename DataAlloc, typename IndexAlloc>
bool operator!=(const OrderedMap<Key, T, Hash, Pred, DataAlloc, IndexAlloc>& a,
                const OrderedMap<Key, T, Hash, Pred, DataAlloc, IndexAlloc>& b){
	if(a.size()==b.size())
    	return false;
	return !std::equal(a.cbegin(), a.cend(), b.cbegin());
}

template<typename Key, typename T, typename Hash, typename Pred,
         typename DataAlloc, typename IndexAlloc>
void swap(const OrderedMap<Key, T, Hash, Pred, DataAlloc, IndexAlloc>& x,
          const OrderedMap<Key, T, Hash, Pred, DataAlloc, IndexAlloc>& y){
	x.swap(y);
}

#endif //G3_CONTAINERS_H