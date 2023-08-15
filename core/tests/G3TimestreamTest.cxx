#include <G3Test.h>

#include <algorithm>
#include <deque>
#include <vector>

#include <core/G3Timestream.h>

TEST_GROUP(G3Timestream)

TEST(SizeValueConstructor){
	for(std::size_t size=0; size<10; size++){
		for(double value : {0, 1}){
			G3Timestream ts(size, value);
			ENSURE(ts.size()==size, "Timestream should be the correct length");
			for(std::size_t i=0; i<size; i++){
				ENSURE(ts[i]==value, "Each timestream sample should match the initializer");
			}
		}
	}
}

TEST(IteratorConstructor){
	std::vector<double> d1={1.7,3.9,13.2};
	std::deque<double> d2={8.8,-6.5,4.9};
	
	G3Timestream t1(d1.begin(), d1.end());
	G3Timestream t2(d2.begin(), d2.end());
	
	ENSURE(std::equal(d1.begin(), d1.end(), t1.begin()),
	       "Timestream should contain the same data as the container from which it was initialized");
	
	ENSURE(std::equal(d2.begin(), d2.end(), t2.begin()),
	       "Timestream should contain the same data as the container from which it was initialized");
}