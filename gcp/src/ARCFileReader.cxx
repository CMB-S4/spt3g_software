#include <pybindings.h>

#include <G3Module.h>
#include <G3TimeStamp.h>
#include <G3Map.h>
#include <G3Vector.h>
#include <G3Data.h>
#include <dataio.h>
#include <gcp/Experiments.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/filesystem.hpp>

#include <string.h>
#include <arpa/inet.h>

#include <deque>
#include <vector>

#ifdef __FreeBSD__
#include <sys/endian.h>
#endif

#ifdef __APPLE__

#include <libkern/OSByteOrder.h>

#define htobe16(x) OSSwapHostToBigInt16(x)
#define htole16(x) OSSwapHostToLittleInt16(x)
#define be16toh(x) OSSwapBigToHostInt16(x)
#define le16toh(x) OSSwapLittleToHostInt16(x)

#define htobe32(x) OSSwapHostToBigInt32(x)
#define htole32(x) OSSwapHostToLittleInt32(x)
#define be32toh(x) OSSwapBigToHostInt32(x)
#define le32toh(x) OSSwapLittleToHostInt32(x)

#define htobe64(x) OSSwapHostToBigInt64(x)
#define htole64(x) OSSwapHostToLittleInt64(x)
#define be64toh(x) OSSwapBigToHostInt64(x)
#define le64toh(x) OSSwapLittleToHostInt64(x)

#endif

enum {
	ARC_SIZE_RECORD = 0,
	ARC_ARRAYMAP_RECORD = 1,
	ARC_FRAME_RECORD = 2
};

enum {
	/* Flags */
	REG_EXC = 0x100,

	/* Types */
	REG_UTC = 0x200, /* MJD + ms */
	REG_BOOL = 0x800,
	REG_CHAR = 0x1000,
	REG_UCHAR = 0x2000,
	REG_SHORT = 0x4000,
	REG_USHORT = 0x8000,
	REG_INT = 0x10000,
	REG_UINT = 0x20000,
	REG_FLOAT = 0x40000,
	REG_DOUBLE = 0x80000,
};

#define REG_TYPEMASK	0xffe01

class ARCFileReader : public G3Module {
public:
	ARCFileReader(const std::string &path,
	    Experiment experiment=Experiment::SPT);
	ARCFileReader(const std::vector<std::string> & filename,
	    Experiment experiment=Experiment::SPT);
	virtual ~ARCFileReader() {}

	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);

private:
	void StartFile(const std::string & path);
	
	boost::iostreams::filtering_istream stream_;

	struct block_stats {
		int flags;
		int mode;
		int base;
		int addr;
		int dim[3];
		std::string axes_help;
		int offset;
		size_t width;
	};
	typedef std::map<std::string, block_stats> board_params;
	typedef std::map<std::string, board_params> arr_template;
	std::map<std::string, arr_template> array_map_;

	off_t frame_length_;
	uint64_t ms_jiffie_base_;

	void ParseArrayMap(uint8_t *buffer, size_t size);
	G3FrameObjectPtr GCPToFrameObject(uint8_t *buffer,
	    const struct block_stats &block, int depth = 0,
	    int base_offset = 0);
	std::deque<std::string> filename_;
	std::string cur_file_;

	Experiment experiment;
	G3TimePtr GCPToTime(uint8_t *buffer, off_t offset);

	SET_LOGGER("ARCFileReader");
};


ARCFileReader::ARCFileReader(const std::string &path,
    Experiment experiment) : experiment(experiment)
{
	if (experiment == Experiment::SPT || experiment == Experiment::BK) {
		ms_jiffie_base_ = G3Units::ms;
	} else if (experiment == Experiment::PB) {
		ms_jiffie_base_ = 86400/INT_MAX;
	} else {
		log_fatal("Unrecognized Experiment");
	}

	boost::filesystem::path fpath(path);
	if (!boost::filesystem::exists(fpath) ||
	    !boost::filesystem::is_regular_file(fpath))
		log_fatal("Could not find file %s", path.c_str());
	StartFile(path);
}



ARCFileReader::ARCFileReader(const std::vector<std::string> &filename,
    Experiment experiment) : experiment(experiment)
{
	if (experiment == Experiment::SPT || experiment == Experiment::BK) {
		ms_jiffie_base_ = G3Units::ms;
	} else if (experiment == Experiment::PB) {
		ms_jiffie_base_ = 86400/INT_MAX;
	} else {
		log_fatal("Unrecognized Experiment");
	}

	if (filename.size() == 0)
		log_fatal("Empty file list provided to G3Reader");

	for (auto i = filename.begin(); i != filename.end(); i++){
		boost::filesystem::path fpath(*i);
		if (!boost::filesystem::exists(fpath) ||
		    !boost::filesystem::is_regular_file(fpath)) 
			log_fatal("Could not find file %s", i->c_str());
		filename_.push_back(*i);
	}

	const std::string path = filename_.front();
	filename_.pop_front();
	StartFile(path);

}


void ARCFileReader::StartFile(const std::string & path)
{
	int32_t size, opcode;
	uint8_t *buffer;

	// Open file, including whatever decompression/network access/etc.
	// may be required
	g3_istream_from_path(stream_, path);
	cur_file_ = path;

	/*
	 * GCP ARC file header:
	 *  ARC_SIZE_RECORD message (we ignore this)
	 *  ARC_ARRAYMAP_RECORD: definition of fields
	 *  ARC_FRAME_RECORD: data
	 */

	// First size record
	stream_.read((char *)&size, sizeof(size));
	size = ntohl(size);
	stream_.read((char *)&opcode, sizeof(opcode));
	opcode = ntohl(opcode);

	if (opcode != ARC_SIZE_RECORD)
		log_fatal("No ARC_SIZE_RECORD at beginning of %s",
		    path.c_str());
	if (size != 12)
		log_fatal("Incorrectly sized ARC_SIZE_RECORD (%d)", size);
	stream_.read((char *)&size, sizeof(size)); /* Skip size field */

	// Get array map
	stream_.read((char *)&size, sizeof(size));
	size = ntohl(size) - 8;
	stream_.read((char *)&opcode, sizeof(opcode));
	opcode = ntohl(opcode);
	if (opcode != ARC_ARRAYMAP_RECORD)
		log_fatal("No ARC_ARRAYMAP_RECORD at beginning of %s",
		    path.c_str());

	buffer = new uint8_t[size];
	stream_.read((char *)buffer, size);
	if (stream_.eof()) {
		delete [] buffer;
		log_fatal("%s truncated; unable to read register map",
		    path.c_str());
	}
	if (!stream_.good()) {
		delete [] buffer;
		log_fatal("Read error on %s while reading register map",
		    path.c_str());
	}

	ParseArrayMap(buffer, size);
	delete [] buffer;
}


void ARCFileReader::ParseArrayMap(uint8_t *buf, size_t size)
{
	int32_t version;
	uint16_t ntemplates;
	off_t off = 0, block_offset = 0;

	// Check serialization version
	memcpy(&version, buf + off, sizeof(version));
	version = ntohl(version); off += sizeof(version);
	if (version != 1)
		log_fatal("Unknown ARC array map version %d in %s\n", version,
		    cur_file_.c_str());

	memcpy(&ntemplates, buf + off, sizeof(ntemplates));
	ntemplates = ntohs(ntemplates); off += 2;

	for (uint16_t i = 0; i < ntemplates; i++) {
		uint16_t namelen, nboard;
		char template_name[255];
		arr_template templ;

		/* Template name */
		memcpy(&namelen, buf + off, sizeof(namelen));
		namelen = ntohs(namelen); off += 2;

		memcpy(template_name, buf + off, namelen);
		template_name[namelen] = '\0';
		off += namelen;

		/* Boards in this template */
		memcpy(&nboard, buf + off, sizeof(nboard));
		nboard = ntohs(nboard); off+= 2;

		/* Add implicit frame board */
		{
			board_params board;
			struct block_stats block_info;

			block_info.mode = 0;
			block_info.base = 0;
			block_info.addr = 0;
			block_info.dim[0] = 1;
			block_info.dim[1] = 0;
			block_info.dim[2] = 0;
			block_info.axes_help = "";

			block_info.offset = block_offset;
			block_info.flags = REG_UINT; block_offset += 4;
			board["status"] = block_info;

			block_info.offset = block_offset;
			block_info.flags = REG_UCHAR; block_offset += 1;
			board["received"] = block_info;

			block_info.offset = block_offset;
			block_info.flags = REG_UINT; block_offset += 4;
			board["nsnap"] = block_info;

			block_info.offset = block_offset;
			block_info.flags = REG_UINT; block_offset += 4;
			board["record"] = block_info;

			block_info.offset = block_offset;
			block_info.flags = REG_UTC; block_offset += 8;
			board["utc"] = block_info;

			if (experiment == Experiment::BK) {
				block_info.offset = block_offset;
				block_info.flags = REG_UINT; block_offset += 4;
				board["lst"] = block_info;
			}

			block_info.offset = block_offset;
			block_info.flags = REG_UINT; block_offset += 4;
			board["features"] = block_info;

			block_info.offset = block_offset;
			block_info.flags = REG_UINT; block_offset += 4;
			board["markSeq"] = block_info;

			templ["frame"] = board;
		}

		for (int j = 0; j < nboard; j++) {
			uint16_t nblock, block;
			uint32_t bases[4];
			char board_name[255];
			board_params board;

			/* Board name */
			memcpy(&namelen, buf + off, sizeof(namelen));
			namelen = ntohs(namelen); off += 2;

			memcpy(board_name, buf + off, namelen);
			board_name[namelen] = '\0';
			off += namelen;

			/* Register blocks */
			memcpy(&nblock, buf + off, sizeof(nblock));
			nblock = ntohs(nblock); off += 2;

			/* Add implicit board status block */
			{
				struct block_stats block_info;

				block_info.mode = 0;
				block_info.base = 0;
				block_info.addr = 0;
				block_info.dim[0] = 1;
				block_info.dim[1] = 0;
				block_info.dim[2] = 0;
				block_info.axes_help = "";
				block_info.offset = block_offset;
				block_info.flags = REG_UINT; block_offset += 4;
				board["status"] = block_info;
			}

			for (block = 0; block < nblock; block++) {
				int ltmp[7];
				char axes_help[255], block_name[255];
				struct block_stats block_info;

				/* Block name */
				memcpy(&namelen, buf + off, sizeof(namelen));
				namelen = ntohs(namelen); off += 2;
				memcpy(block_name, buf + off, namelen);
				block_name[namelen] = '\0';
				off += namelen;
				
				/* Block parameters */
				if (experiment == Experiment::BK) {
					// BK has one fewer dim registers
					memcpy(ltmp, buf + off, sizeof(ltmp) - sizeof(ltmp[0]));
					off += sizeof(ltmp) - sizeof(ltmp[0]);
				} else {
					memcpy(ltmp, buf + off, sizeof(ltmp));
					off += sizeof(ltmp);
				}
				block_info.flags = ntohl(ltmp[0]);
				block_info.mode = ntohl(ltmp[1]);
				block_info.base = ntohl(ltmp[2]);
				block_info.addr = ntohl(ltmp[3]);
				block_info.dim[0] = ntohl(ltmp[4]);
				block_info.dim[1] = ntohl(ltmp[5]);
				if ((experiment == Experiment::SPT) ||
					(experiment == Experiment::PB)) {
					block_info.dim[2] = ntohl(ltmp[6]);
				} else {
					block_info.dim[2] = 0;
				}
				block_info.offset = block_offset;

				/* Get type width */
				switch (block_info.flags & REG_TYPEMASK) {
				case REG_UTC:
					block_info.width = 8;
					break;
				case REG_BOOL: 
				case REG_CHAR: 
				case REG_UCHAR: 
					block_info.width = 1;
					break;
				case REG_SHORT:
				case REG_USHORT:
					block_info.width = 2;
					break;
				case REG_INT:
				case REG_UINT:
				case REG_FLOAT:
					block_info.width = 4;
					break;
				case REG_DOUBLE:
					block_info.width = 8;
					break;
				default:
					log_fatal("Unparseable flags %#x for "
					    "register %s.%s.%s in file %s\n",
					    block_info.flags, template_name,
					    board_name, block_name,
					    cur_file_.c_str());
				}

				/* Skip unarchived registers */
				if (!(block_info.flags & REG_EXC)) {
					if (experiment == Experiment::BK) {
						block_offset += block_info.width *
							((block_info.dim[0] == 0) ? 1 :
							  block_info.dim[0]) *
							((block_info.dim[1] == 0) ? 1 :
							  block_info.dim[1]);
					} else {
						block_offset += block_info.width *
							((block_info.dim[0] == 0) ? 1 :
							  block_info.dim[0]) *
							((block_info.dim[1] == 0) ? 1 :
							  block_info.dim[1]) *
							((block_info.dim[2] == 0) ? 1 :
							  block_info.dim[2]);
					}
				}

				if ((experiment == Experiment::SPT) ||
					(experiment == Experiment::PB)) {
					/* Human-readable type info */
					memcpy(&namelen, buf + off, sizeof(namelen));
					namelen = ntohs(namelen); off += 2;
					memcpy(axes_help, buf + off, namelen);
					axes_help[namelen] = '\0';
					off += namelen;
					block_info.axes_help = axes_help;
				} else {
					block_info.axes_help = "";
				}

				/* Add to map iff actually archived */
				if (!(block_info.flags & REG_EXC))
					board[block_name] = block_info;
			}

			if ((experiment == Experiment::SPT) ||
				(experiment == Experiment::PB)) {
				// Ignored VME addresses bases
				memcpy(bases, buf + off, sizeof(bases));
				off += sizeof(bases);
			} else if (experiment == Experiment::BK) {
				// Inter-frame whitespace
				off += 16;
			}

			templ[board_name] = board;
		}

		array_map_[template_name] = templ;
	}

	frame_length_ = block_offset;
}

G3TimePtr ARCFileReader::GCPToTime(uint8_t *buffer, off_t offset)
{
	uint32_t mjd, ms;

	memcpy(&mjd, buffer + offset, sizeof(mjd));
	memcpy(&ms, buffer + offset + 4, sizeof(ms));
	mjd = le32toh(mjd);
	ms = le32toh(ms);

	mjd -= 40587;	/* MJD to UNIX time */
	if (uint64_t(ms)*ms_jiffie_base_ > 86400ULL*100000000ULL)
		log_error("Fast time value %d longer than 1 day (%lf seconds)",
		    ms, uint64_t(ms)*ms_jiffie_base_/G3Units::s);
	return G3TimePtr(new G3Time(
	   uint64_t(mjd)*86400ULL*100000000ULL + uint64_t(ms)*ms_jiffie_base_));
}

static int GCP8ToInt(uint8_t *buffer, off_t offset)
{
	return int(buffer[offset]);
}

static int GCP16ToInt(uint8_t *buffer, off_t offset)
{
	uint16_t val;

	memcpy(&val, buffer + offset, sizeof(val));
	val = le16toh(val);

	return int(val);
}

static int GCP32ToInt(uint8_t *buffer, off_t offset)
{
	uint32_t val;

	memcpy(&val, buffer + offset, sizeof(val));
	val = le32toh(val);

	return int(val);
}

static double GCP32ToFloat(uint8_t *buffer, off_t offset)
{
	uint32_t val;
	float fval;

	memcpy(&val, buffer + offset, sizeof(val));
	val = le32toh(val);
	memcpy(&fval, &val, sizeof(val));

	return double(fval);
}

static double GCP64ToFloat(uint8_t *buffer, off_t offset)
{
	uint64_t val;
	double fval;

	memcpy(&val, buffer + offset, sizeof(val));
	val = le64toh(val);
	memcpy(&fval, &val, sizeof(val));

	return fval;
}

G3FrameObjectPtr ARCFileReader::GCPToFrameObject(uint8_t *buffer,
    const struct ARCFileReader::block_stats &block, int depth, int base_offset)
{
	if (base_offset == 0)
		base_offset = block.offset;

	int max_depth = 2;
	if (experiment == Experiment::BK)
		max_depth = 1;

	/* Check for two- and three-dimensional case (use vector of vectors) */
	if (depth < max_depth && block.dim[depth+1] != 0) {
		size_t stride = 1;
		G3VectorFrameObjectPtr root(new G3VectorFrameObject);

		for (int i = depth+1; i < max_depth+1; i++)
			if (block.dim[i] != 0)
				stride *= block.dim[i];

		for (int i = 0; i < block.dim[depth]; i++) {
			root->push_back(GCPToFrameObject(buffer, block, depth+1,
			    base_offset));
			base_offset += block.width*stride;
		}

		return root;
	}

	/* 1D case */
	if (block.dim[depth] > 0 && !(depth == 0 && block.dim[depth] == 1)) {
		switch (block.flags & REG_TYPEMASK) {
		case REG_UTC: {
			G3VectorFrameObjectPtr root(new G3VectorFrameObject);
			for (int i = 0; i < block.dim[depth]; i++) {
				root->push_back(GCPToTime(buffer, base_offset));
				base_offset += block.width;
			}
			return root; }
		case REG_CHAR: 
		case REG_UCHAR:
			if (buffer[base_offset + block.dim[depth] - 1] == '\0')
				return G3StringPtr(new G3String(
				    (char *)&buffer[base_offset]));
			// else fallthrough
		case REG_BOOL: {
			G3VectorIntPtr vec(new G3VectorInt);
			for (int i = 0; i < block.dim[depth]; i++) {
				vec->push_back(GCP8ToInt(buffer, base_offset));
				base_offset += block.width;
			}
			return vec; }
		case REG_SHORT:
		case REG_USHORT: {
			G3VectorIntPtr vec(new G3VectorInt);
			for (int i = 0; i < block.dim[depth]; i++) {
				vec->push_back(GCP16ToInt(buffer, base_offset));
				base_offset += block.width;
			}
			return vec; }
		case REG_INT:
		case REG_UINT: {
			G3VectorIntPtr vec(new G3VectorInt);
			for (int i = 0; i < block.dim[depth]; i++) {
				vec->push_back(GCP32ToInt(buffer, base_offset));
				base_offset += block.width;
			}
			return vec; }
		case REG_FLOAT: {
			G3VectorDoublePtr vec(new G3VectorDouble);
			for (int i = 0; i < block.dim[depth]; i++) {
				vec->push_back(GCP32ToFloat(buffer,
				    base_offset));
				base_offset += block.width;
			}
			return vec; }
		case REG_DOUBLE: {
			G3VectorDoublePtr vec(new G3VectorDouble);
			for (int i = 0; i < block.dim[depth]; i++) {
				vec->push_back(GCP64ToFloat(buffer,
				    base_offset));
				base_offset += block.width;
			}
			return vec; }
		default:
			log_fatal("Unknown type");
		}
	}

	switch (block.flags & REG_TYPEMASK) {
	case REG_UTC:
		return GCPToTime(buffer, base_offset);
		break;
	case REG_BOOL: 
	case REG_CHAR: 
	case REG_UCHAR: 
		return G3IntPtr(new G3Int(GCP8ToInt(buffer, base_offset)));
		break;
	case REG_SHORT:
	case REG_USHORT:
		return G3IntPtr(new G3Int(GCP16ToInt(buffer, base_offset)));
		break;
	case REG_INT:
	case REG_UINT:
		return G3IntPtr(new G3Int(GCP32ToInt(buffer, base_offset)));
		break;
	case REG_FLOAT:
		return G3DoublePtr(new G3Double(GCP32ToFloat(buffer,
		    base_offset)));
		break;
	case REG_DOUBLE:
		return G3DoublePtr(new G3Double(GCP64ToFloat(buffer,
		    base_offset)));
		break;
	default:
		log_fatal("Unknown type");
	}

	return (G3FrameObjectPtr());
}

void ARCFileReader::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
{
	G3FramePtr outframe(new G3Frame(G3Frame::GcpSlow));
	int32_t size, opcode;
	uint8_t *buffer;
	off_t off;

	while (stream_.peek() == EOF) {
		if (filename_.size() > 0) {
			const std::string path = filename_.front();
			filename_.pop_front();
			StartFile(path);
		} else {
			return;
		}
	}

	stream_.read((char *)&size, sizeof(size));
	size = ntohl(size) - 8;
	stream_.read((char *)&opcode, sizeof(opcode));
	opcode = ntohl(opcode);
	if (opcode != ARC_FRAME_RECORD)
		log_fatal("Message not an ARC_FRAME_RECORD mid-file");
	if (size != frame_length_)
		log_fatal("%zd-byte frame does not match expected length (%zd)",
		    (size_t)size, (size_t)frame_length_);
	buffer = new uint8_t[size];
	stream_.read((char *)buffer, size);

	for (auto temp = array_map_.begin(); temp != array_map_.end(); temp++) {
		G3MapFrameObjectPtr templ(new G3MapFrameObject);
		for (auto board = temp->second.begin();
		    board != temp->second.end(); board++) {
			G3MapFrameObjectPtr boarddat(new G3MapFrameObject);
			for (auto reg = board->second.begin();
			    reg != board->second.end(); reg++) {
				(*boarddat)[reg->first] = GCPToFrameObject(
				    buffer, reg->second);
			}
			(*templ)[board->first] = boarddat;
		}
		outframe->Put(temp->first, templ);
	}

	out.push_back(outframe);

	delete [] buffer;
}

PYBINDINGS("gcp") {
	using namespace boost::python;

	// Instead of EXPORT_G3MODULE since there are two constructors
	class_<ARCFileReader, bases<G3Module>, boost::shared_ptr<ARCFileReader>,
	    boost::noncopyable>("ARCFileReader",
	    "Read GCP archive file (or files if you pass an iterable of paths)."
	    " For non-SPT ARC file formats, please set Experiment to the "
	    " appropriate value.",
		init<std::string, optional<Experiment> >(
		        args("filename", "experiment")))
		.def(init<std::vector<std::string>, optional<Experiment> >(
		        args("filename", "experiment")))
		.def_readonly("__g3module__", true)
	;
}

