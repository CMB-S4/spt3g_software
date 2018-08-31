#include <arpa/inet.h>
#include <sstream>
#include <pybindings.h>
#include <dfmux/HardwareMap.h>

#include <container_pybindings.h>

std::string DfMuxChannelMapping::Description() const
{
	std::ostringstream s;
	s <<
	    "IP: " << ((ntohl(board_ip) >> 24) & 0xff) << "." <<
	     ((ntohl(board_ip) >> 16) & 0xff) << "." <<
	     ((ntohl(board_ip) >> 8) & 0xff) << "." <<
	     (ntohl(board_ip) & 0xff) << ", " <<
	    "Board: " << board_serial << " (slot " << board_slot << " crate "<< crate_serial<< "), " <<
	    "Module (1-indexed): " << (module + 1) <<
	       ", Channel (1-indexed): " << (channel + 1);
	return s.str();
}

std::string DfMuxChannelMapping::Summary() const
{
	std::ostringstream s;
	if (crate_serial != -1)
		s << crate_serial << "_" << board_slot << "/" << (module + 1) << "/" << (channel + 1);
	else
		s << board_serial << "/" << (module + 1) << "/" << (channel + 1);
	return s.str();
}

template <class A> void DfMuxChannelMapping::serialize(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("board_ip", board_ip);
	ar & make_nvp("board_serial", board_serial);
	ar & make_nvp("board_slot", board_slot);
	if (v > 1)
		ar & make_nvp("crate_serial", crate_serial);
	else
		crate_serial = 0;
	ar & make_nvp("module", module);
	ar & make_nvp("channel", channel);
}

G3_SERIALIZABLE_CODE(DfMuxChannelMapping);
G3_SERIALIZABLE_CODE(DfMuxWiringMap);

PYBINDINGS("dfmux") {
	EXPORT_FRAMEOBJECT(DfMuxChannelMapping, init<>(),
	  "Bolometer wiring information. Module and channel IDs are stored "
	  "zero-indexed, but be aware that they often printed one-indexed "
	  "for compatibility with pydfmux.")
	    .def_readwrite("board_ip", &DfMuxChannelMapping::board_ip,
	     "IP Address of the board, encoded as an int using struct")
	    .def_readwrite("board_serial", &DfMuxChannelMapping::board_serial,
	     "Serial number of the readout board to which this channel is "
	     "attached.")
	    .def_readwrite("board_slot", &DfMuxChannelMapping::board_slot,
	     "Crate slot of the board to which this channel is attached or -1 "
	     "if the board is not in a crate.")
	    .def_readwrite("crate_serial", &DfMuxChannelMapping::crate_serial,
	     "Serial number of the crate in which the readout board is housed "
	     "or -1 if the board is not in a crate.")
	    .def_readwrite("module", &DfMuxChannelMapping::module,
	     "0-indexed module/SQUID ID of the channel")
	    .def_readwrite("channel", &DfMuxChannelMapping::channel,
	     "0-indexed channel number on the parent module/SQUID")
	;

	register_g3map<DfMuxWiringMap>("DfMuxWiringMap", "Mapping from "
	   "logical detector ID string (same as used in timestreams) to wiring "
	   "information (the board, module, and channel to which a given "
	   "detector is connected)");
}

