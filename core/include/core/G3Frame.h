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
		Ephemeris = 'E',
		LightCurve = 'L',
		Statistics = 'R',
		EndProcessing = 'Z',
		None = 'N',
	};
		
	G3Frame(FrameType type = G3Frame::None);

	FrameType type;
	
	// Source filename (if output from G3Reader)
	std::string _filename;

	// Add and remove data from the frame. Note that retrieved data is
	// const. Get<> is like operator [], except that it does a dynamic
	// cast to a type for you.
	//
	// If the data cannot be retrieved:
	// - operator [] will return a null pointer if the element is not
	//   present
	// - Get<> will either log_fatal (default) or return a null pointer if
	//   the requested object is not in the frame or cannot be cast to the
	//   given type. Whether it does log_fatal() or returns a NULL pointer
	//   is set by the value of exception_on_error -- if true, it will
	//   log_fatal(), printing an informative error and throwing an
	//   exception.
	G3FrameObjectConstPtr operator [](const std::string &) const;
	template<typename T> std::shared_ptr<const T> Get(
	    const std::string &, bool exception_on_error = true) const;
	void Put(const std::string &name, G3FrameObjectConstPtr);
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
	template <typename T> void saves(T &os) const;
	template <typename T> void loads(T &is);

	// Serialize (or deserialize) from a cereal archive.
	template <class A> void save(A &ar, unsigned) const;
	template <class A> void load(A &ar, unsigned);

	// JSON representation of a frame (if JSON output is enabled, otherwise
	// outputs valid JSON telling you that JSON output isn't enabled)
	std::string asJSON() const;

	// Routines for handling the stored serialized copies of data.
	// These are all const because they only manipulate caches and
	// so change only performance rather than data contents of the frame.

	// Drop all serialized data either for already-decoded objects
	// (default) or all objects after decoding them (if decode_all is
	// true). Saves memory at the expense of CPU time if reserialized.
	void DropBlobs(bool decode_all=false) const;

	// Force immediate serialization of all objects. Will save some
	// CPU time later during serialization of the frame in exchange for
	// spending the exact same amount of CPU time right now.
	void GenerateBlobs(bool drop_objects=false) const;

	// Drop all decoded objects in favor of their serialized copies, where
	// those serialized copies already exist. Saves memory for frames
	// about to be written at the expense of CPU time to re-decode them
	// if they are accessed again later.
	void DropObjects() const;

private:
	struct blob_container {
		G3FrameObjectConstPtr frameobject;
		std::shared_ptr<std::vector<char> > blob;
	};
	static void blob_decode(struct blob_container &blob);
	static void blob_encode(struct blob_container &blob);
	typedef std::unordered_map<std::string, blob_container> G3MapType;
	friend std::ostream& operator<<(std::ostream& os, const G3Frame &);
	mutable G3MapType map_;

	SET_LOGGER("G3Frame");
};

G3_POINTERS(G3Frame);
CEREAL_CLASS_VERSION(G3Frame, 1);
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(G3Frame, cereal::specialization::member_load_save);

template <typename T>
std::shared_ptr<const T> G3Frame::Get(const std::string &name,
  bool exception_on_error) const
{
	std::shared_ptr<const T> ret =
	    std::dynamic_pointer_cast<const T>((*this)[name]);

	if (exception_on_error && !ret)
		log_fatal("Requesting key %s %s", name.c_str(),
		    Has(name) ? "of the wrong type" : "not in frame");

	return ret;
}

template <typename T>
bool G3Frame::Has(const std::string &name) const
{
	return !!std::dynamic_pointer_cast<const T>((*this)[name]);
}

std::ostream& operator<<(std::ostream& os, const G3FrameObject &);
std::ostream& operator<<(std::ostream& os, const G3Frame &);
std::ostream& operator<<(std::ostream& os, const G3Frame::FrameType &);

#endif

