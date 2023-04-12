#include <G3Test.h>

#include <core/G3Timestream.h>

#include <boost/make_shared.hpp>

TEST_GROUP(G3TimestreamMap)

template<typename SampleType>
void testMakeCompact(){
	std::vector<std::string> keys={"a","b","c","d"};
	G3Time t1(0), t2(3);
	const std::size_t nSamples=4;
	G3TimestreamMap tsm=G3TimestreamMap::MakeCompact<SampleType>(keys, nSamples, t1, t2);
	
	ENSURE_EQUAL(tsm.size(),keys.size(), "MakeCompact should produce the correct number of timestreams");
	ENSURE_EQUAL(tsm.GetStartTime(),t1, "MakeCompact should produce a map with the given start time");
	ENSURE_EQUAL(tsm.GetStopTime(),t2, "MakeCompact should produce a map with the given stop time");
	ENSURE_EQUAL(tsm.NSamples(),nSamples, "MakeCompact should produce a map with the given number of samples");
	ENSURE(tsm.CheckAlignment(), "MakeCompact should produce aligned timestreams");
}

TEST(MakeCompactDouble){ testMakeCompact<double>(); }
TEST(MakeCompactFloat){ testMakeCompact<float>(); }
TEST(MakeCompactInt32){ testMakeCompact<int32_t>(); }
TEST(MakeCompactInt64){ testMakeCompact<int64_t>(); }

TEST(MakeCompactUnits){
	std::vector<std::string> keys={"a","b","c","d"};
	G3Time t1(0), t2(3);
	const std::size_t nSamples=4;
	G3TimestreamMap tsm=G3TimestreamMap::MakeCompact<double>(keys, nSamples, t1, t2, G3Timestream::Voltage);
	
	for(const auto& item : tsm)
		ENSURE_EQUAL(item.second->units, G3Timestream::Voltage, "All timestreams' units should be set");
}

TEST(MakeCompactCompressed){
	std::vector<std::string> keys={"a","b","c","d"};
	G3Time t1(0), t2(3);
	const std::size_t nSamples=4;
	G3TimestreamMap tsm=G3TimestreamMap::MakeCompact<int32_t>(keys, nSamples, t1, t2, G3Timestream::Counts, 6);
	
	for(const auto& item : tsm)
		ENSURE_EQUAL(item.second->GetCompressionLevel(), 6, "All timestreams should be set to use compression");
}

template<typename SampleType>
void testMakeCompactExisting(){
	std::vector<std::string> keys={"a","b","c","d"};
	G3Time t1(0), t2(3);
	const std::size_t nSamples=4;
	boost::shared_ptr<SampleType[]> data=boost::make_shared<SampleType[]>(keys.size()*nSamples);
	//construct a recognizable data pattern
	for(std::size_t i=0; i<keys.size(); i++){
		for(std::size_t j=0; j<nSamples; j++){
			data[i*nSamples + j] = (SampleType)(20*i*nSamples + 2*j)*(j%2?1:-1);
		}
	}
	
	G3TimestreamMap tsm=G3TimestreamMap::MakeCompact<SampleType>(keys, nSamples, data, t1, t2);
	
	ENSURE_EQUAL(tsm.size(),keys.size(), "MakeCompact should produce the correct number of timestreams");
	ENSURE_EQUAL(tsm.GetStartTime(),t1, "MakeCompact should produce a map with the given start time");
	ENSURE_EQUAL(tsm.GetStopTime(),t2, "MakeCompact should produce a map with the given stop time");
	ENSURE_EQUAL(tsm.NSamples(),nSamples, "MakeCompact should produce a map with the given number of samples");
	ENSURE(tsm.CheckAlignment(), "MakeCompact should produce aligned timestreams");
	ENSURE_EQUAL(data.use_count(), keys.size()+1, "Each timestream should hold a reference to the data block");
	
	for(std::size_t i=0; i<keys.size(); i++){
		const auto& key=keys[i];
		auto it=tsm.find(key);
		ENSURE(it!=tsm.end(), "Each expected key shjould be present in the map");
		ENSURE((bool)it->second, "Each timestream pointer should be valid");
		const G3Timestream& ts=*it->second;
		ENSURE(ts.size()==nSamples);
		for(std::size_t j=0; j<nSamples; j++){
			ENSURE_EQUAL(ts[j],(SampleType)(20*i*nSamples + 2*j)*(j%2?1:-1));
		}
	}
}

TEST(MakeCompactExistingDouble){ testMakeCompactExisting<double>(); }
TEST(MakeCompactExistingFloat){ testMakeCompactExisting<float>(); }
TEST(MakeCompactExistingInt32){ testMakeCompactExisting<int32_t>(); }
TEST(MakeCompactExistingInt64){ testMakeCompactExisting<int64_t>(); }

TEST(SetStartTime){
	std::vector<std::string> keys={"a","b","c","d"};
	G3Time t1(0), t2(3), t3(47);
	const std::size_t nSamples=4;
	G3TimestreamMap tsm=G3TimestreamMap::MakeCompact<double>(keys, nSamples, t1, t2);
	
	ENSURE_EQUAL(tsm.GetStartTime(),t1, "MakeCompact should produce a map with the given start time");
	tsm.SetStartTime(t3);
	ENSURE_EQUAL(tsm.GetStartTime(),t3, "SetStartTime should change the timestream start time");
	ENSURE(tsm.CheckAlignment(), "All start times should match after SetStartTime");
}

TEST(SetStopTime){
	std::vector<std::string> keys={"a","b","c","d"};
	G3Time t1(0), t2(3), t3(47);
	const std::size_t nSamples=4;
	G3TimestreamMap tsm=G3TimestreamMap::MakeCompact<double>(keys, nSamples, t1, t2);
	
	ENSURE_EQUAL(tsm.GetStopTime(),t2, "MakeCompact should produce a map with the given stop time");
	tsm.SetStopTime(t3);
	ENSURE_EQUAL(tsm.GetStopTime(),t3, "SetStopTime should change the timestream start time");
	ENSURE(tsm.CheckAlignment(), "All stop times should match after SetStopTime");
}

TEST(SetUnits){
	std::vector<std::string> keys={"a","b","c","d"};
	G3Time t1(0), t2(3);
	const std::size_t nSamples=4;
	G3TimestreamMap tsm=G3TimestreamMap::MakeCompact<double>(keys, nSamples, t1, t2);
	
	tsm.SetUnits(G3Timestream::Voltage);
	for(const auto& key : keys){
		ENSURE_EQUAL(tsm[key]->units, G3Timestream::Voltage, "SetUnits should set the units for each timestream");
	}
}

TEST(SetCompression){
	std::vector<std::string> keys={"a","b","c","d"};
	G3Time t1(0), t2(3);
	const std::size_t nSamples=4;
	G3TimestreamMap tsm=G3TimestreamMap::MakeCompact<double>(keys, nSamples, t1, t2);
	
	unsigned int compressionLevel=5;
	tsm.SetFLACCompression(5);
	for(const auto& key : keys){
		ENSURE_EQUAL(tsm[key]->GetCompressionLevel(), compressionLevel, "SetFLACCompression should set the compression level for each timestream");
	}
}