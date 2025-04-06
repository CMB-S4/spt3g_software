#include <G3Test.h>

#include <memory>
#include <type_traits>

#include <core/containers.h>

TEST_GROUP(OrderedMap)

using TestMap=OrderedMap<std::string,std::shared_ptr<std::string>>;
using TestPairType=std::pair<std::string,std::shared_ptr<std::string>>;
using TestValueType=std::pair<const std::string,std::shared_ptr<std::string>>;

static_assert(std::is_default_constructible<TestMap>::value,
              "OrderedMap must be default constructible");

static_assert(std::is_constructible<TestMap, std::size_t>::value,
              "OrderedMap must be constructible from a number of hash buckets");

static_assert(std::is_constructible<TestMap, std::size_t, std::hash<std::string>>::value,
              "OrderedMap must be constructible from a number of hash buckets and a hasher");

static_assert(std::is_constructible<TestMap, 
                                    std::size_t, std::hash<std::string>, std::equal_to<std::string>
                                   >::value,
              "OrderedMap must be constructible from a number of hash buckets, a hasher, "
              "and a key comprator");

static_assert(std::is_constructible<TestMap, 
                                    std::size_t, std::hash<std::string>, std::equal_to<std::string>,
                                    std::allocator<std::pair<const std::string, 
                                                             std::shared_ptr<std::string>>>
                                   >::value,
              "OrderedMap must be constructible from a number of hash buckets, a hasher, "
              "a key comprator, and a data allocator");

static_assert(std::is_constructible<TestMap, 
                                    std::size_t, std::hash<std::string>, std::equal_to<std::string>,
                                    std::allocator<std::pair<const std::string, 
                                                             std::shared_ptr<std::string>>>,
                                    std::allocator<std::pair<const std::string, 
                                                             typename std::list<std::pair<const std::string, std::shared_ptr<std::string>>>::iterator>>
                                   >::value,
              "OrderedMap must be constructible from a number of hash buckets, a hasher, "
              "a key comprator, a data allocator, and an index allocator");

static_assert(std::is_constructible<TestMap, 
                                    std::list<TestPairType>::iterator, std::list<TestPairType>::iterator
                                   >::value,
              "OrderedMap must be constructible from a pair of data iterators");

static_assert(std::is_constructible<TestMap, 
                                    std::list<TestPairType>::iterator, std::list<TestPairType>::iterator,
                                    std::size_t
                                   >::value,
              "OrderedMap must be constructible from a pair of data iterators and a number of hash buckets");

static_assert(std::is_constructible<TestMap, 
                                    std::list<TestPairType>::iterator, std::list<TestPairType>::iterator,
                                    std::size_t, std::hash<std::string>
                                   >::value,
              "OrderedMap must be constructible from a pair of data iterators, a number of hash buckets, "
              "and a hasher");

static_assert(std::is_constructible<TestMap, 
                                    std::list<TestPairType>::iterator, std::list<TestPairType>::iterator,
                                    std::size_t, std::hash<std::string>, std::equal_to<std::string>
                                   >::value,
              "OrderedMap must be constructible from a pair of data iterators, a number of hash buckets, "
              "a hasher, and a comparator");

static_assert(std::is_constructible<TestMap, 
                                    std::list<TestPairType>::iterator, std::list<TestPairType>::iterator,
                                    std::size_t, std::hash<std::string>, std::equal_to<std::string>,
                                    std::allocator<std::pair<const std::string, 
                                                             std::shared_ptr<std::string>>>
                                   >::value,
              "OrderedMap must be constructible from a pair of data iterators, a number of hash buckets, "
              "a hasher, a comparator, and a data allocator");

static_assert(std::is_constructible<TestMap, 
                                    std::list<TestPairType>::iterator, std::list<TestPairType>::iterator,
                                    std::size_t, std::hash<std::string>, std::equal_to<std::string>,
                                    std::allocator<std::pair<const std::string, 
                                                             std::shared_ptr<std::string>>>,
                                    std::allocator<std::pair<const std::string, 
                                                             typename std::list<std::pair<const std::string, std::shared_ptr<std::string>>>::iterator>>
                                   >::value,
              "OrderedMap must be constructible from a pair of data iterators, a number of hash buckets, "
              "a hasher, a comparator, a data allocator, and an index allocator");

static_assert(std::is_constructible<TestMap, 
                                    std::initializer_list<TestValueType>
                                   >::value,
              "OrderedMap must be constructible from an initializer list");

static_assert(std::is_constructible<TestMap, 
                                    std::initializer_list<TestValueType>,
                                    std::size_t
                                   >::value,
              "OrderedMap must be constructible from an initializer list and a number of hash buckets");

static_assert(std::is_constructible<TestMap, 
                                    std::initializer_list<TestValueType>,
                                    std::size_t, std::hash<std::string>
                                   >::value,
              "OrderedMap must be constructible from an initializer list, a number of hash buckets, "
              "and a hasher");

static_assert(std::is_constructible<TestMap, 
                                    std::initializer_list<TestValueType>,
                                    std::size_t, std::hash<std::string>, std::equal_to<std::string>
                                   >::value,
              "OrderedMap must be constructible from an initializer list, a number of hash buckets, "
              "a hasher, and a comparator");

static_assert(std::is_constructible<TestMap, 
                                    std::initializer_list<TestValueType>,
                                    std::size_t, std::hash<std::string>, std::equal_to<std::string>,
                                    std::allocator<std::pair<const std::string, 
                                                             std::shared_ptr<std::string>>>
                                   >::value,
              "OrderedMap must be constructible from an initializer list, a number of hash buckets, "
              "a hasher, a comparator, and a data allocator");

static_assert(std::is_constructible<TestMap, 
                                    std::initializer_list<TestValueType>,
                                    std::size_t, std::hash<std::string>, std::equal_to<std::string>,
                                    std::allocator<std::pair<const std::string, 
                                                             std::shared_ptr<std::string>>>,
                                    std::allocator<std::pair<const std::string, 
                                                             typename std::list<std::pair<const std::string, std::shared_ptr<std::string>>>::iterator>>
                                   >::value,
              "OrderedMap must be constructible from an initializer list, a number of hash buckets, "
              "a hasher, a comparator, a data allocator, and an index allocator");

static_assert(std::is_constructible<TestMap, 
                                    std::allocator<std::pair<const std::string, 
                                                             std::shared_ptr<std::string>>>
                                   >::value,
              "OrderedMap must be constructible from a data allocator");

static_assert(std::is_constructible<TestMap, 
                                    std::allocator<std::pair<const std::string, 
                                                             std::shared_ptr<std::string>>>,
                                    std::allocator<std::pair<const std::string, 
                                                             typename std::list<std::pair<const std::string, std::shared_ptr<std::string>>>::iterator>>
                                   >::value,
              "OrderedMap must be constructible from a data allocator and an index allocator");

static_assert(std::is_constructible<TestMap,
                                    TestMap&,
                                    std::allocator<std::pair<const std::string,
                                                             std::shared_ptr<std::string>>>
                                   >::value,
              "OrderedMap must be constructible from a reference to another map and a data allocator");

static_assert(std::is_constructible<TestMap,
                                    TestMap&,
                                    std::allocator<std::pair<const std::string,
                                                             std::shared_ptr<std::string>>>,
                                    std::allocator<std::pair<const std::string, 
                                                             typename std::list<std::pair<const std::string, std::shared_ptr<std::string>>>::iterator>>
                                   >::value,
              "OrderedMap must be constructible from a reference to another map, a data allocator, "
              "and an index allocator");

static_assert(std::is_constructible<TestMap,
                                    std::size_t,
                                    std::allocator<std::pair<const std::string,
                                                             std::shared_ptr<std::string>>>
                                   >::value,
              "OrderedMap must be constructible from a number of hash buckets and a data allocator");

static_assert(std::is_constructible<TestMap,
                                    std::size_t,
                                    std::allocator<std::pair<const std::string,
                                                             std::shared_ptr<std::string>>>,
                                    std::allocator<std::pair<const std::string, 
                                                             typename std::list<std::pair<const std::string, std::shared_ptr<std::string>>>::iterator>>
                                   >::value,
              "OrderedMap must be constructible from a number of hash buckets, a data allocator, "
              "and an index allocator");

static_assert(std::is_constructible<TestMap,
                                    std::size_t,
                                    std::hash<std::string>
                                   >::value,
              "OrderedMap must be constructible from a number of hash buckets and a hasher");

static_assert(std::is_constructible<TestMap,
                                    std::size_t,
                                    std::hash<std::string>,
                                    std::allocator<std::pair<const std::string,
                                                             std::shared_ptr<std::string>>>
                                   >::value,
              "OrderedMap must be constructible from a number of hash buckets, a hasher, "
              "and a data allocator");

static_assert(std::is_constructible<TestMap,
                                    std::size_t,
                                    std::hash<std::string>,
                                    std::allocator<std::pair<const std::string,
                                                             std::shared_ptr<std::string>>>,
                                    std::allocator<std::pair<const std::string, 
                                                             typename std::list<std::pair<const std::string, std::shared_ptr<std::string>>>::iterator>>
                                   >::value,
              "OrderedMap must be constructible from a number of hash buckets, a hasher, "
              "a data allocator, and an index allocator");

static_assert(std::is_constructible<TestMap,
                                    std::list<TestPairType>::iterator, std::list<TestPairType>::iterator,
                                    std::size_t,
                                    std::allocator<std::pair<const std::string,
                                                             std::shared_ptr<std::string>>>
                                   >::value,
              "OrderedMap must be constructible from a pair of data iterators, "
              "a number of hash buckets, and a data allocator");

static_assert(std::is_constructible<TestMap,
                                    std::list<TestPairType>::iterator, std::list<TestPairType>::iterator,
                                    std::size_t,
                                    std::allocator<std::pair<const std::string,
                                                             std::shared_ptr<std::string>>>,
                                    std::allocator<std::pair<const std::string, 
                                                             typename std::list<std::pair<const std::string, std::shared_ptr<std::string>>>::iterator>>
                                   >::value,
              "OrderedMap must be constructible from a pair of data iterators, "
              "a number of hash buckets, a data allocator, and an index allocator");

static_assert(std::is_constructible<TestMap,
                                    std::list<TestPairType>::iterator, std::list<TestPairType>::iterator,
                                    std::size_t,
                                    std::hash<std::string>,
                                    std::allocator<std::pair<const std::string,
                                                             std::shared_ptr<std::string>>>
                                   >::value,
              "OrderedMap must be constructible from a pair of data iterators, "
              "a number of hash buckets, a hasher, and a data allocator");

static_assert(std::is_constructible<TestMap,
                                    std::list<TestPairType>::iterator, std::list<TestPairType>::iterator,
                                    std::size_t,
                                    std::hash<std::string>,
                                    std::allocator<std::pair<const std::string,
                                                             std::shared_ptr<std::string>>>,
                                    std::allocator<std::pair<const std::string, 
                                                             typename std::list<std::pair<const std::string, std::shared_ptr<std::string>>>::iterator>>
                                   >::value,
              "OrderedMap must be constructible from a pair of data iterators, "
              "a number of hash buckets, a hasher, a data allocator, and an index allocator");

static_assert(std::is_constructible<TestMap,
                                    std::initializer_list<TestValueType>,
                                    std::size_t,
                                    std::allocator<std::pair<const std::string,
                                                             std::shared_ptr<std::string>>>
                                   >::value,
              "OrderedMap must be constructible from an initializer list, "
              "a number of hash buckets, and a data allocator");

static_assert(std::is_constructible<TestMap,
                                    std::initializer_list<TestValueType>,
                                    std::size_t,
                                    std::allocator<std::pair<const std::string,
                                                             std::shared_ptr<std::string>>>,
                                    std::allocator<std::pair<const std::string, 
                                                             typename std::list<std::pair<const std::string, std::shared_ptr<std::string>>>::iterator>>
                                   >::value,
              "OrderedMap must be constructible from an initializer list, "
              "a number of hash buckets, a data allocator, and an index allocator");

static_assert(std::is_constructible<TestMap,
                                    std::initializer_list<TestValueType>,
                                    std::size_t,
                                    std::hash<std::string>,
                                    std::allocator<std::pair<const std::string,
                                                             std::shared_ptr<std::string>>>,
                                    std::allocator<std::pair<const std::string, 
                                                             typename std::list<std::pair<const std::string, std::shared_ptr<std::string>>>::iterator>>
                                   >::value,
              "OrderedMap must be constructible from an initializer list, "
              "a number of hash buckets, a hasher, a data allocator, and an index allocator");

static_assert(std::is_constructible<TestMap,
                                    std::initializer_list<TestValueType>,
                                    std::size_t,
                                    std::hash<std::string>,
                                    std::allocator<std::pair<const std::string,
                                                             std::shared_ptr<std::string>>>
                                   >::value,
              "OrderedMap must be constructible from an initializer list, "
              "a number of hash buckets, a hasher, and a data allocator");

static_assert(std::is_copy_constructible<TestMap>::value,
              "OrderedMap must be copy constructible");

static_assert(std::is_move_constructible<TestMap>::value,
              "OrderedMap must be move constructible");

static_assert(std::is_copy_assignable<TestMap>::value,
              "OrderedMap must be copy assignable");

static_assert(std::is_move_assignable<TestMap>::value,
              "OrderedMap must be move assignable");

static_assert(std::is_assignable<TestMap, std::initializer_list<TestValueType>>::value,
              "OrderedMap must be assignable from an initializer list");

void check_contents(const TestMap& map, const std::initializer_list<std::pair<std::string,std::string>>& exp){
	ENSURE_EQUAL(map.size(), exp.size());
	//check via iteration interface
	auto mit=map.begin();
	for(auto eit=exp.begin(), end=exp.end(); eit!=end; mit++,eit++){
		ENSURE_EQUAL(mit->first, eit->first, "Keys must match");
		ENSURE_EQUAL(*mit->second, eit->second, "Values must match");
	}
	//check via lookup interface
	for(auto eit=exp.begin(), end=exp.end(); eit!=end; mit++,eit++){
		auto mit=map.find(eit->first);
		ENSURE(mit!=map.end());
		ENSURE_EQUAL(mit->first, eit->first, "Keys must match");
		ENSURE_EQUAL(*mit->second, eit->second, "Values must match");
	}
}

TEST(IteratorRangeConstruction){
	std::list<TestPairType> raw={{"foo",std::make_shared<std::string>("bar")},
	                             {"baz",std::make_shared<std::string>("quux")},
	                             {"xen",std::make_shared<std::string>("hom")}};
	TestMap map(raw.begin(), raw.end());
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
}

TEST(InitializerListConstruction){
	std::list<TestPairType> raw={{"foo",std::make_shared<std::string>("bar")},
	                             {"baz",std::make_shared<std::string>("quux")},
	                             {"xen",std::make_shared<std::string>("hom")}};
	TestMap map({{"foo",std::make_shared<std::string>("bar")},
	             {"baz",std::make_shared<std::string>("quux")},
	             {"xen",std::make_shared<std::string>("hom")}});
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
}

TEST(InitializerListAssignement){
	std::initializer_list<TestValueType> il{{"foo",std::make_shared<std::string>("bar")},
	                                        {"baz",std::make_shared<std::string>("quux")},
	                                        {"xen",std::make_shared<std::string>("hom")}};
	TestMap map({{"drel",std::make_shared<std::string>("plugh")}});
	map = il;
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
}

TEST(Emplace){
	TestMap map;
	
	auto res=map.emplace("foo", std::make_shared<std::string>("bar"));
	ENSURE_EQUAL(res.first->first, "foo");
	ENSURE_EQUAL(*res.first->second, "bar");
	ENSURE_EQUAL(res.second, true);
	check_contents(map, {{"foo","bar"}});
	res=map.emplace(std::make_pair("baz", std::make_shared<std::string>("quux")));
	ENSURE_EQUAL(res.first->first, "baz");
	ENSURE_EQUAL(*res.first->second, "quux");
	ENSURE_EQUAL(res.second, true);
	check_contents(map, {{"foo","bar"},{"baz","quux"}});
	res=map.emplace(std::piecewise_construct, 
	                std::make_tuple("xen"), 
	                std::make_tuple(std::make_shared<std::string>("hom")));
	ENSURE_EQUAL(res.first->first, "xen");
	ENSURE_EQUAL(*res.first->second, "hom");
	ENSURE_EQUAL(res.second, true);
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
	
	//atempting to emplace with an existing key should do nothing
	res=map.emplace("foo", std::make_shared<std::string>("drel"));
	ENSURE_EQUAL(res.first->first, "foo");
	ENSURE_EQUAL(*res.first->second, "bar");
	ENSURE_EQUAL(res.second, false);
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
}

TEST(EmplaceHint){
	TestMap map;
	
	auto res=map.emplace_hint(map.begin(), "foo", std::make_shared<std::string>("bar"));
	ENSURE_EQUAL(res->first, "foo");
	ENSURE_EQUAL(*res->second, "bar");
	check_contents(map, {{"foo","bar"}});
	res=map.emplace_hint(map.begin(), std::make_pair("baz", std::make_shared<std::string>("quux")));
	ENSURE_EQUAL(res->first, "baz");
	ENSURE_EQUAL(*res->second, "quux");
	check_contents(map, {{"foo","bar"},{"baz","quux"}});
	res=map.emplace_hint(map.begin(), std::piecewise_construct, 
	                     std::make_tuple("xen"), 
	                     std::make_tuple(std::make_shared<std::string>("hom")));
	ENSURE_EQUAL(res->first, "xen");
	ENSURE_EQUAL(*res->second, "hom");
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
	
	//atempting to emplace with an existing key should do nothing
	res=map.emplace_hint(map.begin(), "foo", std::make_shared<std::string>("drel"));
	ENSURE_EQUAL(res->first, "foo");
	ENSURE_EQUAL(*res->second, "bar");
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
}

TEST(InsertRef){
	TestMap map;
	
	TestValueType p0("foo", std::make_shared<std::string>("bar"));
	auto res=map.insert((const TestValueType&)p0);
	ENSURE_EQUAL(res.first->first, "foo");
	ENSURE_EQUAL(*res.first->second, "bar");
	ENSURE_EQUAL(res.second, true);
	check_contents(map, {{"foo","bar"}});
	TestValueType p1("baz", std::make_shared<std::string>("quux"));
	res=map.insert((const TestValueType&)p1);
	ENSURE_EQUAL(res.first->first, "baz");
	ENSURE_EQUAL(*res.first->second, "quux");
	ENSURE_EQUAL(res.second, true);
	check_contents(map, {{"foo","bar"},{"baz","quux"}});
	TestValueType p2("xen", std::make_shared<std::string>("hom"));
	res=map.insert((const TestValueType&)p2);
	ENSURE_EQUAL(res.first->first, "xen");
	ENSURE_EQUAL(*res.first->second, "hom");
	ENSURE_EQUAL(res.second, true);
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
	
	//atempting to emplace with an existing key should do nothing
	TestValueType p3("foo", std::make_shared<std::string>("drel"));
	res=map.insert((const TestValueType&)p3);
	ENSURE_EQUAL(res.first->first, "foo");
	ENSURE_EQUAL(*res.first->second, "bar");
	ENSURE_EQUAL(res.second, false);
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
}

TEST(InsertRRef){
	TestMap map;
	
	TestValueType p0("foo", std::make_shared<std::string>("bar"));
	auto res=map.insert(std::move(p0));
	ENSURE_EQUAL(res.first->first, "foo");
	ENSURE_EQUAL(*res.first->second, "bar");
	ENSURE_EQUAL(res.second, true);
	check_contents(map, {{"foo","bar"}});
	TestValueType p1("baz", std::make_shared<std::string>("quux"));
	res=map.insert(std::move(p1));
	ENSURE_EQUAL(res.first->first, "baz");
	ENSURE_EQUAL(*res.first->second, "quux");
	ENSURE_EQUAL(res.second, true);
	check_contents(map, {{"foo","bar"},{"baz","quux"}});
	TestValueType p2("xen", std::make_shared<std::string>("hom"));
	res=map.insert(std::move(p2));
	ENSURE_EQUAL(res.first->first, "xen");
	ENSURE_EQUAL(*res.first->second, "hom");
	ENSURE_EQUAL(res.second, true);
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
	
	//atempting to emplace with an existing key should do nothing
	TestValueType p3("foo", std::make_shared<std::string>("drel"));
	res=map.insert(std::move(p3));
	ENSURE_EQUAL(res.first->first, "foo");
	ENSURE_EQUAL(*res.first->second, "bar");
	ENSURE_EQUAL(res.second, false);
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
}

TEST(InsertConvertible){
	struct Thing{
		std::string k;
		std::shared_ptr<std::string> v;
		Thing(const std::string& k, const std::string& v):
		k(k),v(std::make_shared<std::string>(v)){}
		operator TestValueType() const{
			return TestValueType(k,v);
		}
	};
	TestMap map;
	
	Thing p0("foo", "bar");
	auto res=map.insert(p0);
	ENSURE_EQUAL(res.first->first, "foo");
	ENSURE_EQUAL(*res.first->second, "bar");
	ENSURE_EQUAL(res.second, true);
	check_contents(map, {{"foo","bar"}});
	Thing p1("baz", "quux");
	res=map.insert(p1);
	ENSURE_EQUAL(res.first->first, "baz");
	ENSURE_EQUAL(*res.first->second, "quux");
	ENSURE_EQUAL(res.second, true);
	check_contents(map, {{"foo","bar"},{"baz","quux"}});
	Thing p2("xen", "hom");
	res=map.insert(p2);
	ENSURE_EQUAL(res.first->first, "xen");
	ENSURE_EQUAL(*res.first->second, "hom");
	ENSURE_EQUAL(res.second, true);
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
	
	//atempting to emplace with an existing key should do nothing
	Thing p3("foo", "drel");
	res=map.insert(p3);
	ENSURE_EQUAL(res.first->first, "foo");
	ENSURE_EQUAL(*res.first->second, "bar");
	ENSURE_EQUAL(res.second, false);
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
}

TEST(InsertHintRef){
	TestMap map;
	
	TestValueType p0("foo", std::make_shared<std::string>("bar"));
	auto res=map.insert(map.begin(), (const TestValueType&)p0);
	ENSURE_EQUAL(res->first, "foo");
	ENSURE_EQUAL(*res->second, "bar");
	check_contents(map, {{"foo","bar"}});
	TestValueType p1("baz", std::make_shared<std::string>("quux"));
	res=map.insert(map.begin(), (const TestValueType&)p1);
	ENSURE_EQUAL(res->first, "baz");
	ENSURE_EQUAL(*res->second, "quux");
	check_contents(map, {{"foo","bar"},{"baz","quux"}});
	TestValueType p2("xen", std::make_shared<std::string>("hom"));
	res=map.insert(map.begin(), (const TestValueType&)p2);
	ENSURE_EQUAL(res->first, "xen");
	ENSURE_EQUAL(*res->second, "hom");
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
	
	//atempting to emplace with an existing key should do nothing
	TestValueType p3("foo", std::make_shared<std::string>("drel"));
	res=map.insert(map.begin(), (const TestValueType&)p3);
	ENSURE_EQUAL(res->first, "foo");
	ENSURE_EQUAL(*res->second, "bar");
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
}

TEST(InsertHintRRef){
	TestMap map;
	
	TestValueType p0("foo", std::make_shared<std::string>("bar"));
	auto res=map.insert(map.begin(), std::move(p0));
	ENSURE_EQUAL(res->first, "foo");
	ENSURE_EQUAL(*res->second, "bar");
	check_contents(map, {{"foo","bar"}});
	TestValueType p1("baz", std::make_shared<std::string>("quux"));
	res=map.insert(map.begin(), std::move(p1));
	ENSURE_EQUAL(res->first, "baz");
	ENSURE_EQUAL(*res->second, "quux");
	check_contents(map, {{"foo","bar"},{"baz","quux"}});
	TestValueType p2("xen", std::make_shared<std::string>("hom"));
	res=map.insert(map.begin(), std::move(p2));
	ENSURE_EQUAL(res->first, "xen");
	ENSURE_EQUAL(*res->second, "hom");
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
	
	//atempting to emplace with an existing key should do nothing
	TestValueType p3("foo", std::make_shared<std::string>("drel"));
	res=map.insert(map.begin(), std::move(p3));
	ENSURE_EQUAL(res->first, "foo");
	ENSURE_EQUAL(*res->second, "bar");
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
}

TEST(InsertHintConvertible){
	struct Thing{
		std::string k;
		std::shared_ptr<std::string> v;
		Thing(const std::string& k, const std::string& v):
		k(k),v(std::make_shared<std::string>(v)){}
		operator TestValueType() const{
			return TestValueType(k,v);
		}
	};
	TestMap map;
	
	Thing p0("foo", "bar");
	auto res=map.insert(map.begin(), p0);
	ENSURE_EQUAL(res->first, "foo");
	ENSURE_EQUAL(*res->second, "bar");
	check_contents(map, {{"foo","bar"}});
	Thing p1("baz", "quux");
	res=map.insert(map.begin(), p1);
	ENSURE_EQUAL(res->first, "baz");
	ENSURE_EQUAL(*res->second, "quux");
	check_contents(map, {{"foo","bar"},{"baz","quux"}});
	Thing p2("xen", "hom");
	res=map.insert(map.begin(), p2);
	ENSURE_EQUAL(res->first, "xen");
	ENSURE_EQUAL(*res->second, "hom");
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
	
	//atempting to emplace with an existing key should do nothing
	Thing p3("foo", "drel");
	res=map.insert(map.begin(), p3);
	ENSURE_EQUAL(res->first, "foo");
	ENSURE_EQUAL(*res->second, "bar");
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
}

TEST(InsertIteratorRange){
	std::list<TestPairType> raw={{"foo",std::make_shared<std::string>("bar")},
	                             {"baz",std::make_shared<std::string>("quux")},
	                             {"xen",std::make_shared<std::string>("hom")}};
	TestMap map;
	map.insert(raw.begin(), raw.end());
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
}

TEST(InsertInitializerList){
	std::list<TestValueType> raw={{"foo",std::make_shared<std::string>("bar")},
	                              {"baz",std::make_shared<std::string>("quux")},
	                              {"xen",std::make_shared<std::string>("hom")}};
	TestMap map;
	map.insert({{"foo",std::make_shared<std::string>("bar")},
	            {"baz",std::make_shared<std::string>("quux")},
	            {"xen",std::make_shared<std::string>("hom")}});
	check_contents(map, {{"foo","bar"},{"baz","quux"},{"xen","hom"}});
}

TEST(EraseIterator){
	TestMap map({{"foo",std::make_shared<std::string>("bar")},
	             {"baz",std::make_shared<std::string>("quux")},
	             {"xen",std::make_shared<std::string>("hom")}});
	auto res=map.erase(map.find("baz"));
	ENSURE_EQUAL(map.size(), 2U);
	ENSURE_EQUAL(res->first, "xen");
	ENSURE_EQUAL(*res->second, "hom");
	check_contents(map, {{"foo","bar"},{"xen","hom"}});
	
	res=map.erase(map.find("xen"));
	ENSURE_EQUAL(map.size(), 1U);
	ENSURE(res==map.end());
	check_contents(map, {{"foo","bar"}});
	
	res=map.erase(map.find("foo"));
	ENSURE_EQUAL(map.size(), 0U);
	ENSURE(res==map.end());
}

TEST(EraseConstIterator){
	TestMap map({{"foo",std::make_shared<std::string>("bar")},
	             {"baz",std::make_shared<std::string>("quux")},
	             {"xen",std::make_shared<std::string>("hom")}});
	auto res=map.erase((TestMap::const_iterator)map.find("baz"));
	ENSURE_EQUAL(map.size(), 2U);
	ENSURE_EQUAL(res->first, "xen");
	ENSURE_EQUAL(*res->second, "hom");
	check_contents(map, {{"foo","bar"},{"xen","hom"}});
	
	res=map.erase((TestMap::const_iterator)map.find("xen"));
	ENSURE_EQUAL(map.size(), 1U);
	ENSURE(res==map.end());
	check_contents(map, {{"foo","bar"}});
	
	res=map.erase((TestMap::const_iterator)map.find("foo"));
	ENSURE_EQUAL(map.size(), 0U);
	ENSURE(res==map.end());
}

TEST(EraseKey){
	TestMap map({{"foo",std::make_shared<std::string>("bar")},
	             {"baz",std::make_shared<std::string>("quux")},
	             {"xen",std::make_shared<std::string>("hom")}});
	auto res=map.erase("baz");
	ENSURE_EQUAL(map.size(), 2U);
	ENSURE_EQUAL(res, 1);
	check_contents(map, {{"foo","bar"},{"xen","hom"}});
	
	res=map.erase("drel");
	ENSURE_EQUAL(map.size(), 2U);
	ENSURE_EQUAL(res, 0);
	check_contents(map, {{"foo","bar"},{"xen","hom"}});
	
	res=map.erase("xen");
	ENSURE_EQUAL(map.size(), 1U);
	ENSURE_EQUAL(res, 1);
	check_contents(map, {{"foo","bar"}});
	
	res=map.erase("foo");
	ENSURE_EQUAL(map.size(), 0U);
	ENSURE_EQUAL(res, 1);
}

TEST(EraseIteratorRange){
	TestMap map({{"foo",std::make_shared<std::string>("bar")},
	             {"baz",std::make_shared<std::string>("quux")},
	             {"xen",std::make_shared<std::string>("hom")}});
	auto res=map.erase(map.begin(), map.find("xen"));
	ENSURE_EQUAL(res->first, "xen");
	ENSURE_EQUAL(*res->second, "hom");
	check_contents(map, {{"xen","hom"}});
	res=map.erase(map.begin(), map.end());
	ENSURE(res==map.end());
	ENSURE(map.empty());
}

TEST(Swap){
	TestMap m1({{"foo",std::make_shared<std::string>("bar")},
	            {"baz",std::make_shared<std::string>("quux")}});
	TestMap m2({{"xen",std::make_shared<std::string>("hom")},
	            {"drel",std::make_shared<std::string>("plugh")}});
	m1.swap(m2);
	check_contents(m2, {{"foo","bar"},{"baz","quux"}});
	check_contents(m1, {{"xen","hom"},{"drel","plugh"}});
}

TEST(Find){
	TestMap map({{"foo",std::make_shared<std::string>("bar")},
	             {"baz",std::make_shared<std::string>("quux")},
	             {"xen",std::make_shared<std::string>("hom")}});
	auto check_found=[&](std::string k, std::string v){
		auto it=map.find(k);
		ENSURE(it!=map.end());
		ENSURE_EQUAL(it->first, k);
		ENSURE_EQUAL(*it->second, v);
	};
	check_found("foo","bar");
	check_found("baz","quux");
	check_found("xen","hom");
	ENSURE(map.find("drel")==map.end());
}

TEST(FindConst){
	const TestMap map({{"foo",std::make_shared<std::string>("bar")},
	                   {"baz",std::make_shared<std::string>("quux")},
	                   {"xen",std::make_shared<std::string>("hom")}});
	auto check_found=[&](std::string k, std::string v){
		auto it=map.find(k);
		ENSURE(it!=map.end());
		ENSURE_EQUAL(it->first, k);
		ENSURE_EQUAL(*it->second, v);
	};
	check_found("foo","bar");
	check_found("baz","quux");
	check_found("xen","hom");
	ENSURE(map.find("drel")==map.end());
}

namespace convertible_key{
	struct Key{
		std::string s;
		Key(const char* d):s(d){}
		operator std::string() const{ return s; }
	};
	bool operator==(const std::string& s, const Key& t){ return s==t.s; }
	std::ostream& operator<<(std::ostream& os, const Key& t){
		return os << "Key(" << t.s << ')';
	}
}

TEST(FindConvertible){
	using namespace convertible_key;
	TestMap map({{"foo",std::make_shared<std::string>("bar")},
	             {"baz",std::make_shared<std::string>("quux")},
	             {"xen",std::make_shared<std::string>("hom")}});
	auto check_found=[&](const Key& k, std::string v){
		auto it=map.find(k);
		ENSURE(it!=map.end());
		ENSURE_EQUAL(it->first, k);
		ENSURE_EQUAL(*it->second, v);
	};
	check_found("foo","bar");
	check_found("baz","quux");
	check_found("xen","hom");
	ENSURE(map.find("drel")==map.end());
}

TEST(FindConstConvertible){
	using namespace convertible_key;
	const TestMap map({{"foo",std::make_shared<std::string>("bar")},
	                   {"baz",std::make_shared<std::string>("quux")},
	                   {"xen",std::make_shared<std::string>("hom")}});
	auto check_found=[&](const Key& k, std::string v){
		auto it=map.find(k);
		ENSURE(it!=map.end());
		ENSURE_EQUAL(it->first, k);
		ENSURE_EQUAL(*it->second, v);
	};
	check_found("foo","bar");
	check_found("baz","quux");
	check_found("xen","hom");
	ENSURE(map.find("drel")==map.end());
}

TEST(EqualRange){
	TestMap map({{"foo",std::make_shared<std::string>("bar")},
	             {"baz",std::make_shared<std::string>("quux")},
	             {"xen",std::make_shared<std::string>("hom")}});
	auto check_found=[&](std::string k, std::string v){
		auto r=map.equal_range(k);
		ENSURE(r.first!=map.end());
		ENSURE_EQUAL(std::distance(r.first,r.second),1,
		             "With unique keys equal_range should never return ranges longer than 1");
		ENSURE_EQUAL(r.first->first, k);
		ENSURE_EQUAL(*r.first->second, v);
	};
	check_found("foo","bar");
	check_found("baz","quux");
	check_found("xen","hom");
	{
		auto r=map.equal_range(std::string("drel"));
		ENSURE(r.first==map.end());
		ENSURE(r.second==map.end());
	}
}

TEST(EqualRangeConst){
	const TestMap map({{"foo",std::make_shared<std::string>("bar")},
	                   {"baz",std::make_shared<std::string>("quux")},
	                   {"xen",std::make_shared<std::string>("hom")}});
	auto check_found=[&](std::string k, std::string v){
		auto r=map.equal_range(k);
		ENSURE(r.first!=map.end());
		ENSURE_EQUAL(std::distance(r.first,r.second),1,
		             "With unique keys equal_range should never return ranges longer than 1");
		ENSURE_EQUAL(r.first->first, k);
		ENSURE_EQUAL(*r.first->second, v);
	};
	check_found("foo","bar");
	check_found("baz","quux");
	check_found("xen","hom");
	{
		auto r=map.equal_range(std::string("drel"));
		ENSURE(r.first==map.end());
		ENSURE(r.second==map.end());
	}
}

TEST(EqualRangeConvertible){
	using namespace convertible_key;
	TestMap map({{"foo",std::make_shared<std::string>("bar")},
	             {"baz",std::make_shared<std::string>("quux")},
	             {"xen",std::make_shared<std::string>("hom")}});
	auto check_found=[&](Key k, std::string v){
		auto r=map.equal_range(k);
		ENSURE(r.first!=map.end());
		ENSURE_EQUAL(std::distance(r.first,r.second),1,
		             "With unique keys equal_range should never return ranges longer than 1");
		ENSURE_EQUAL(r.first->first, k);
		ENSURE_EQUAL(*r.first->second, v);
	};
	check_found("foo","bar");
	check_found("baz","quux");
	check_found("xen","hom");
	{
		auto r=map.equal_range(Key("drel"));
		ENSURE(r.first==map.end());
		ENSURE(r.second==map.end());
	}
}

TEST(EqualRangeConvertibleConst){
	using namespace convertible_key;
	const TestMap map({{"foo",std::make_shared<std::string>("bar")},
	                   {"baz",std::make_shared<std::string>("quux")},
	                   {"xen",std::make_shared<std::string>("hom")}});
	auto check_found=[&](Key k, std::string v){
		auto r=map.equal_range(k);
		ENSURE(r.first!=map.end());
		ENSURE_EQUAL(std::distance(r.first,r.second),1,
		             "With unique keys equal_range should never return ranges longer than 1");
		ENSURE_EQUAL(r.first->first, k);
		ENSURE_EQUAL(*r.first->second, v);
	};
	check_found("foo","bar");
	check_found("baz","quux");
	check_found("xen","hom");
	{
		auto r=map.equal_range(Key("drel"));
		ENSURE(r.first==map.end());
		ENSURE(r.second==map.end());
	}
}

TEST(Indexing){
	TestMap map;
	map["foo"]=std::make_shared<std::string>("bar");
	check_contents(map, {{"foo","bar"}});
	ENSURE_EQUAL(*map["foo"], "bar");
	map["baz"]=std::make_shared<std::string>("quux");
	check_contents(map, {{"foo","bar"},{"baz","quux"}});
	ENSURE_EQUAL(*map["baz"], "quux");
	map["foo"]=std::make_shared<std::string>("xen");
	check_contents(map, {{"foo","xen"},{"baz","quux"}});
	ENSURE_EQUAL(*map["foo"], "xen");
}

TEST(At){
	TestMap map({{"foo",std::make_shared<std::string>("bar")},
	             {"baz",std::make_shared<std::string>("quux")},
	             {"xen",std::make_shared<std::string>("hom")}});
	ENSURE_EQUAL(*map.at("foo"), "bar");
	ENSURE_EQUAL(*map.at("baz"), "quux");
	ENSURE_EQUAL(*map.at("xen"), "hom");
	try{
		map.at("drel");
		FAIL("std::out_of_range should be thrown");
	}
	catch(std::out_of_range& ex){
		//expected
	}
}

TEST(AtConst){
	const TestMap map({{"foo",std::make_shared<std::string>("bar")},
	                   {"baz",std::make_shared<std::string>("quux")},
	                   {"xen",std::make_shared<std::string>("hom")}});
	ENSURE_EQUAL(*map.at("foo"), "bar");
	ENSURE_EQUAL(*map.at("baz"), "quux");
	ENSURE_EQUAL(*map.at("xen"), "hom");
	try{
		map.at("drel");
		FAIL("std::out_of_range should be thrown");
	}
	catch(std::out_of_range& ex){
		//expected
	}
}