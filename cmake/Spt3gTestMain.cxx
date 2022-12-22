#include <G3Test.h>

#include <memory>

namespace G3Test{
std::ostream& operator<<(std::ostream& os, const TestFailure& tf){
	return os << tf.file << ':' << tf.line << ": " << tf.message;
}

std::ostream& operator<<(std::ostream& os, const TestReport& r){
	return os << r.successes << " pass" << (r.successes==1?" ":"es ")
		<< r.failures << " failur" << (r.failures==1?"e":"es");
}

static TestSuite testSuite;

TestSuite& getTestSuite(){ return testSuite; }

TestRegisterer::TestRegisterer(const std::string& group, const std::string& name, std::function<void()> impl){
	getTestSuite().registerTest(group,name,impl);
}

void TestSuite::registerTest(std::string groupName, std::string testName, std::function<void()> testImpl){
	groups.insert(groupName);
	tests.insert(std::make_pair(groupName+"::"+testName, testImpl));
}

void TestSuite::listTests(){
	for(const auto& t : tests)
		std::cout << t.first << '\n';
}

TestReport TestSuite::runTests(bool verbose){
	TestReport report;
	for(const auto& t : tests)
		report.addResult(runSingleTest(t, verbose));
	return report;
}

TestReport TestSuite::runTests(const std::vector<std::string>& testNames, bool verbose){
	TestReport report;
	for(const auto& name : testNames){
		if(groups.count(name)){
			//run all tests in the group
			std::string prefix=name+"::";
			for(const auto& t : tests){
				if(t.first.find(prefix)==0)
					report.addResult(runSingleTest(t, verbose));
			}
		}
		else{
			auto it=tests.find(name);
			if(it!=tests.end())
				report.addResult(runSingleTest(*it, verbose));
			else
				std::cerr << "Test not found: " << name << std::endl;
		}
	}
	return report;
}

bool TestSuite::runSingleTest(const std::pair<std::string, std::function<void()>>& test, bool verbose){
	struct SilenceOutput{
		SilenceOutput():bufOut(std::cout.rdbuf()),bufErr(std::cerr.rdbuf()){
			std::cout.rdbuf(s.rdbuf());
			std::cerr.rdbuf(s.rdbuf());
		}
		~SilenceOutput(){
			std::cout.rdbuf(bufOut);
			std::cerr.rdbuf(bufErr);
		}
		std::ostringstream s;
		std::streambuf* bufOut;
		std::streambuf* bufErr;
	};
	bool pass=false;
	try{
		std::cout << test.first << ": ";
		std::cout.flush();
		std::unique_ptr<SilenceOutput> silence;
		if(!verbose)
			silence.reset(new SilenceOutput);
		test.second();
		pass=true;
	}catch(const TestFailure& tf){
		std::cout << tf << std::endl;
	}catch(const std::exception& ex){
		std::cout << "Uncaught std::exception: " << ex.what() << std::endl;
	}
	#ifdef USE_PYTHON
	catch(const boost::python::error_already_set& ex){
		PyErr_Print();
	}
	#endif
	catch(...){
		std::cout << "Unidentified flying object" << std::endl;
	}
	if(pass)
		std::cout << "PASS" << std::endl;
	else
		std::cout << "FAIL" << std::endl;
	return pass;
}

} //namespace G3Test

int main(int argc, char* argv[]){
	bool quiet=true;
	bool list=false;
	bool runAll=true;
	std::vector<std::string> testsToRun;

	for(unsigned int i=1; i<argc; i++){
		std::string opt=argv[i];
		if(opt=="--help" || opt=="-h"){
			std::cout << "Usage: \n"
			<< " --all, -a                     Run all tests [default]\n"
			<< " --run-tests, -r [Test Names]  Run only the named tests\n"
			<< "                               Specifying a group name will\n"
			<< "                               run all tests inthat group.\n"
			<< " --list, -l                    List test names\n"
			<< " --quiet, -q                   Suppress output printed by tests\n"
			<< "                               [default]\n"
			<< " --verbose, -v                 Show output printed by tests\n";
			return 0;
		}
		else if(opt=="--quiet" || opt=="-q"){
			quiet=true;
		}
		else if(opt=="--verbose" || opt=="-v"){
			quiet=false;
		}
		else if(opt=="--list" || opt=="-l"){
			list=true;
			runAll=false;
		}
		else if(opt=="--all" || opt=="-a"){
			runAll=true;
		}
		else if(opt=="--run-tests" || opt=="-r"){
			runAll=false;
			i++;
			while(i<argc && argv[i][0]!='-')
				testsToRun.push_back(argv[i++]);
			i--;
		}
		else{
			std::cerr << "Unknown option: " << opt << std::endl;
			return 1;
		}
	}
	
	G3Test::TestSuite& tests=G3Test::getTestSuite();
	G3Test::TestReport report;
	if(list)
		tests.listTests();
	if(!testsToRun.empty())
		report&=tests.runTests(testsToRun,!quiet);
	if(runAll)
		report&=tests.runTests(!quiet);
	if(runAll || !testsToRun.empty())
		std::cout << report << std::endl;
	return ((bool)report?0:1);
}