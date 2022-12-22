#ifndef G3TEST_H_INCLUDED
#define G3TEST_H_INCLUDED

#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <boost/preprocessor/stringize.hpp>

namespace G3Test{

class TestFailure{
public:
	TestFailure(const std::string& m, const std::string& f, unsigned int l):
	message(m),file(f),line(l){}
private:
	std::string message;
	std::string file;
	unsigned int line;
	
	friend std::ostream& operator<<(std::ostream&, const TestFailure& tf);
};

class TestReport{
public:
	TestReport():successes(0),failures(0){}
	void addResult(bool success){ (success?successes:failures)++; }
	explicit operator bool() const{ return failures==0; }
	TestReport& operator&=(const TestReport& r){
		if(&r!=this){
			successes+=r.successes;
			failures+=r.failures;
		}
		return *this;
	}
private:
	unsigned int successes, failures;
	friend std::ostream& operator<<(std::ostream&, const TestReport& r);
};

class TestSuite{
public:
	void registerTest(std::string groupName, std::string testName, std::function<void()> testImpl);
	
	void listTests();
	TestReport runTests(bool verbose);
	TestReport runTests(const std::vector<std::string>& testNames, bool verbose);
private:
	std::set<std::string> groups;
	std::map<std::string, std::function<void()>> tests;
	
	bool runSingleTest(const std::pair<std::string, std::function<void()>>& test, bool verbose);
};

TestSuite& getTestSuite();

class TestRegisterer{
public:
	TestRegisterer(const std::string& group, const std::string& name, std::function<void()> impl);
};

inline void testAssertion(const std::string& file, unsigned int line, 
                          const std::string& predicate, bool predicateValue,
                          const std::string message=""){
	if(!predicateValue)
		throw TestFailure(message.empty()? predicate : predicate+": "+message, file, line);
}

template<typename T1, typename T2>
inline void testEquivalence(const std::string& file, unsigned int line, 
                            const T1& v1, const T2& v2, const std::string& v1s, const std::string& v2s,
                            const std::string message=""){
	if(!(v1==v2)){
		std::ostringstream ss;
		ss << "ENSURE_EQUAL(" << v1s << ", " << v2s << "): " << v1 << " != " << v2;
		if(!message.empty())
		ss << ": " << message;
		throw TestFailure(ss.str(), file, line);
	}
}

} //G3Test

/// Define a test group. use like:
/// TEST_GROUP(MyTests)
#define TEST_GROUP(GROUPNAME) \
namespace{ \
static const char* GetTestGroup(){ return BOOST_PP_STRINGIZE(GROUPNAME); } \
}

/// Define a test. Use like:
/// TEST(MyTest){
///     test code. . .
/// }
/// A test group should always be declared with TEST_GROUP before any tests.
#define TEST(TESTNAME) \
static void test_implementation_ ## TESTNAME(); \
namespace{ \
static G3Test::TestRegisterer register_ ## TESTNAME (GetTestGroup(), \
	BOOST_PP_STRINGIZE(TESTNAME), \
	test_implementation_ ## TESTNAME); \
} \
static void test_implementation_ ## TESTNAME()

#define ENSURE(predicate,...) \
G3Test::testAssertion(__FILE__, __LINE__, BOOST_PP_STRINGIZE(predicate), (predicate), ##__VA_ARGS__);

#define FAIL(...) \
G3Test::testAssertion(__FILE__, __LINE__, "FAILURE", false, ##__VA_ARGS__);

#define ENSURE_EQUAL(left,right,...) \
G3Test::testEquivalence(__FILE__, __LINE__, left, right, BOOST_PP_STRINGIZE(left), BOOST_PP_STRINGIZE(right), ##__VA_ARGS__);

#endif //G3TEST_H_INCLUDED