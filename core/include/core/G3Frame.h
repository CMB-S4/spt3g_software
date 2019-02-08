#ifndef _G3_FRAME_H
#define _G3_FRAME_H

#include <unordered_map>
#include <string>
#include <iostream>
#include <G3.h>
#include <G3Logging.h>

/*
 * The frame is a collection of data indexed by name (basically an
 * std::map<std::string, G3FrameObjectPtr>). Each object stored in the
 * frame must subclass G3FrameObject so that it can be serialized to disk and
 * so that RTTI can be used to discover type information.
 */

class G3FrameObject {
public:
	G3FrameObject() {}
	virtual ~G3FrameObject() {}

	// Return a verbose human-readable description of the data. By default,
	// returns a string with the object's type.
	virtual std::string Description() const;

	// Return a one-line summary of the data shown when frames are printed
	// or for __str__() in python. By default, calls Description().
	virtual std::string Summary() const;

	// Serialization hook for boost. Must be reimplemented for child
	// classes.
	template <class A> void serialize(A &ar, const unsigned v);
};

G3_POINTERS(G3FrameObject);
G3_SERIALIZABLE(G3FrameObject, 1);

class G3Frame {
public:
	// Frame type codes. Indicates type of data in frame.
	enum FrameType {
		Timepoint = 'T',
		Housekeeping = 'H',
		Observation = 'O',
		Scan = 'S',
		Map = 'M',
		InstrumentStatus = 'I',
		Wiring = 'W',
		Calibration = 'C',
		GcpSlow = 'G',
		PipelineInfo = 'P',
		EndProcessing = 'Z',
		None = 'N',
	};
		
	G3Frame(FrameType type = G3Frame::None);

	FrameType type;
	
	// Add and remove data from the frame. Note that retrieved data is
	// const. Will return a null pointer if the requested object is not
	// in the frame. Get<> does RTTI for you.
	G3FrameObjectConstPtr operator [](const std::string &) const;
	void Put(const std::string &name, G3FrameObjectConstPtr);
	template<typename T> boost::shared_ptr<const T> Get(const std::string &) const;
	void Delete(const std::string &);

	// Check if an object is in the frame. Second version makes sure it
	// exists and is the right type.
	bool Has(const std::string &) const;
	template<typename T> bool Has(const std::string &name) const;

	// Get list of keys. This should be an iterator eventually.
	std::vector<std::string> Keys() const;

	// Number of keys in frame
	size_t size() const { return map_.size(); }

	G3Frame &operator = (const G3Frame &);

	// Serialize (or deserialize) from an IO stream.
	template <typename T> void save(T &os) const;
	template <typename T> void load(T &is);

private:
	struct blob_container {
		G3FrameObjectConstPtr frameobject;
		boost::shared_ptr<std::vector<char> > blob;
	};
	static void blob_decode(struct blob_container &blob);
	typedef std::unordered_map<std::string, blob_container> G3MapType;
	friend std::ostream& operator<<(std::ostream& os, const G3Frame &);
	mutable G3MapType map_;

	SET_LOGGER("G3Frame");
};

G3_POINTERS(G3Frame);

template <typename T>
boost::shared_ptr<const T> G3Frame::Get(const std::string &name) const
{
	boost::shared_ptr<const T> ret =
	    boost::dynamic_pointer_cast<const T>((*this)[name]);

	if (!ret)
		log_fatal("Requesting key %s %s", name.c_str(),
		    Has(name) ? "of the wrong type" : "not in frame");

	return ret;
}

template <typename T>
bool G3Frame::Has(const std::string &name) const
{
	return !!boost::dynamic_pointer_cast<const T>((*this)[name]);
}

std::ostream& operator<<(std::ostream& os, const G3FrameObject &);
std::ostream& operator<<(std::ostream& os, const G3Frame &);
std::ostream& operator<<(std::ostream& os, const G3Frame::FrameType &);

#endif

