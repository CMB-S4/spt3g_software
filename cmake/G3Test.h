#ifndef G3TEST_H_INCLUDED
#define G3TEST_H_INCLUDED

#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#ifdef USE_PYTHON
#include <pybindings.h>
#endif

#define STRINGIZE(s) STRINGIZE_DIRECT(s)
#define STRINGIZE_DIRECT(s) #s

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

/// Define a test group. Use like:
/// TEST_GROUP(MyTests)
#define TEST_GROUP(GROUPNAME) \
namespace{ \
static const char* GetTestGroup(){ return STRINGIZE(GROUPNAME); } \
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
	STRINGIZE(TESTNAME), \
	test_implementation_ ## TESTNAME); \
} \
static void test_implementation_ ## TESTNAME()

/// Require that a predicate must be true, or mark the containing test as
/// failed. An optional message may be passed after the predicate.
#define ENSURE(predicate,...) \
G3Test::testAssertion(__FILE__, __LINE__, STRINGIZE(predicate), (predicate), ##__VA_ARGS__);

/// Unconditionally mark the containing test as failed. A message may
/// optionally be passed.
/// This is useful in the presence of other conditional logic, e.g. exception
/// handling, wherein one may wish to treat an undesired control flow path as
/// a test failure.
#define FAIL(...) \
G3Test::testAssertion(__FILE__, __LINE__, "FAILURE", false, ##__VA_ARGS__);

/// Require that two value are equal at runtime, or mark the containing test
/// as failed. An optional message may be passed after the two expressions to
/// be compared.
/// It can be useful to use this variant, rather than writing a simple ENSURE()
/// with a predicate which performs an equality comparison because both the
/// left and right operands and their values are printed in the event of a
/// failure.
#define ENSURE_EQUAL(left,right,...) \
G3Test::testEquivalence(__FILE__, __LINE__, left, right, STRINGIZE(left), STRINGIZE(right), ##__VA_ARGS__);


//=============================================================================
// Implementation details

/// A helper macro to allow condiftional use of python without conditional
/// directives needing to appear inside G3TEST_MAIN_IMPL.
/// This expands to either a catch block or nothing, so to ensure that code is
/// always well-formed, a try block using this should always have at least one
/// other catch block.
#ifdef USE_PYTHON
#define G3TEST_CATCH_PYTHON_ERROR \
catch(const py::error_already_set& ex){ \
	PyErr_Print(); \
}
#else
#define G3TEST_CATCH_PYTHON_ERROR /*no python error handling*/
#endif

/// Users should generally not use this macro directly.
/// It expands to the full implemtnation which needs to be compiled for the
/// test infrastructure, and is present here so that the entire
/// infrastsructure can be supplied b y this one header.
/// To form a complete test executable, this macro should be expanded in
/// exactly one translation unit, which is often handles by a cmake macro.
#define G3TEST_MAIN_IMPL \
namespace G3Test{ \
std::ostream& operator<<(std::ostream& os, const TestFailure& tf){ \
	return os << tf.file << ':' << tf.line << ": " << tf.message; \
} \
 \
std::ostream& operator<<(std::ostream& os, const TestReport& r){ \
	return os << r.successes << " pass" << (r.successes==1?" ":"es ") \
		<< r.failures << " failur" << (r.failures==1?"e":"es"); \
} \
 \
static TestSuite testSuite; \
 \
TestSuite& getTestSuite(){ return testSuite; } \
 \
TestRegisterer::TestRegisterer(const std::string& group, const std::string& name, std::function<void()> impl){ \
	getTestSuite().registerTest(group,name,impl); \
} \
 \
void TestSuite::registerTest(std::string groupName, std::string testName, std::function<void()> testImpl){ \
	groups.insert(groupName); \
	tests.insert(std::make_pair(groupName+"::"+testName, testImpl)); \
} \
 \
void TestSuite::listTests(){ \
	for(const auto& t : tests) \
		std::cout << t.first << '\n'; \
} \
 \
TestReport TestSuite::runTests(bool verbose){ \
	TestReport report; \
	for(const auto& t : tests) \
		report.addResult(runSingleTest(t, verbose)); \
	return report; \
} \
 \
TestReport TestSuite::runTests(const std::vector<std::string>& testNames, bool verbose){ \
	TestReport report; \
	for(const auto& name : testNames){ \
		if(groups.count(name)){ \
			/*run all tests in the group*/ \
			std::string prefix=name+"::"; \
			for(const auto& t : tests){ \
				if(t.first.find(prefix)==0) \
					report.addResult(runSingleTest(t, verbose)); \
			} \
		} \
		else{ \
			auto it=tests.find(name); \
			if(it!=tests.end()) \
				report.addResult(runSingleTest(*it, verbose)); \
			else \
				std::cerr << "Test not found: " << name << std::endl; \
		} \
	} \
	return report; \
} \
 \
bool TestSuite::runSingleTest(const std::pair<std::string, std::function<void()>>& test, bool verbose){ \
	struct SilenceOutput{ \
		SilenceOutput():bufOut(std::cout.rdbuf()),bufErr(std::cerr.rdbuf()){ \
			std::cout.rdbuf(s.rdbuf()); \
			std::cerr.rdbuf(s.rdbuf()); \
		} \
		~SilenceOutput(){ \
			std::cout.rdbuf(bufOut); \
			std::cerr.rdbuf(bufErr); \
		} \
		std::ostringstream s; \
		std::streambuf* bufOut; \
		std::streambuf* bufErr; \
	}; \
	bool pass=false; \
	try{ \
		std::cout << test.first << ": "; \
		std::cout.flush(); \
		std::unique_ptr<SilenceOutput> silence; \
		if(!verbose) \
			silence.reset(new SilenceOutput); \
		test.second(); \
		pass=true; \
	}catch(const TestFailure& tf){ \
		std::cout << tf << std::endl; \
	}catch(const std::exception& ex){ \
		std::cout << "Uncaught std::exception: " << ex.what() << std::endl; \
	} \
	G3TEST_CATCH_PYTHON_ERROR \
	catch(...){ \
		std::cout << "Unidentified flying object" << std::endl; \
	} \
	if(pass) \
		std::cout << "PASS" << std::endl; \
	else \
		std::cout << "FAIL" << std::endl; \
	return pass; \
} \
 \
} /*namespace G3Test*/ \
 \
int main(int argc, char* argv[]){ \
	bool quiet=true; \
	bool list=false; \
	bool runAll=true; \
	std::vector<std::string> testsToRun; \
 \
	for(int i=1; i<argc; i++){ \
		std::string opt=argv[i]; \
		if(opt=="--help" || opt=="-h"){ \
			std::cout << "Usage: \n" \
			<< " --all, -a                     Run all tests [default]\n" \
			<< " --run-tests, -r [Test Names]  Run only the named tests\n" \
			<< "                               Specifying a group name will\n" \
			<< "                               run all tests in that group.\n" \
			<< " --list, -l                    List test names\n" \
			<< " --quiet, -q                   Suppress output printed by tests\n" \
			<< "                               [default]\n" \
			<< " --verbose, -v                 Show output printed by tests\n"; \
			return 0; \
		} \
		else if(opt=="--quiet" || opt=="-q"){ \
			quiet=true; \
		} \
		else if(opt=="--verbose" || opt=="-v"){ \
			quiet=false; \
		} \
		else if(opt=="--list" || opt=="-l"){ \
			list=true; \
			runAll=false; \
		} \
		else if(opt=="--all" || opt=="-a"){ \
			runAll=true; \
		} \
		else if(opt=="--run-tests" || opt=="-r"){ \
			runAll=false; \
			i++; \
			while(i<argc && argv[i][0]!='-') \
				testsToRun.push_back(argv[i++]); \
			i--; \
		} \
		else{ \
			std::cerr << "Unknown option: " << opt << std::endl; \
			return 1; \
		} \
	} \
	 \
	G3Test::TestSuite& tests=G3Test::getTestSuite(); \
	G3Test::TestReport report; \
	if(list) \
		tests.listTests(); \
	if(!testsToRun.empty()) \
		report&=tests.runTests(testsToRun,!quiet); \
	if(runAll) \
		report&=tests.runTests(!quiet); \
	if(runAll || !testsToRun.empty()) \
		std::cout << report << std::endl; \
	return ((bool)report?0:1); \
}

#endif //G3TEST_H_INCLUDED
