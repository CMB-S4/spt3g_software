#ifndef _G3_DATAIO_H
#define _G3_DATAIO_H

#include <string>
#include <vector>
#include <boost/iostreams/filtering_stream.hpp>

typedef boost::iostreams::filtering_istream g3_istream;
typedef boost::iostreams::filtering_ostream g3_ostream;

// Function to obtain a filtering istream from a path/URL/someplace
// Knows how to decompress files if compressed (based on filename extension),
// read files from disk, and from network sockets

int g3_istream_from_path(g3_istream &stream, const std::string &path,
    float timeout=-1.0);

off_t g3_istream_seek(g3_istream &stream, off_t offset);
off_t g3_istream_tell(g3_istream &stream);

// Function to stream data from a char buffer

void g3_istream_from_buffer(g3_istream &stream, const char *buf, size_t len);

// Function to obtain a filtering ostream from a path
// Knows how to compress files (based on filename extension)

void g3_ostream_to_path(g3_ostream &stream, const std::string &path,
    bool append=false, bool counter=false);

// Function to return a count of streamed characters

size_t g3_ostream_count(g3_ostream &stream);

// Function to stream data to a char buffer

void g3_ostream_to_buffer(g3_ostream &stream, std::vector<char> &buf);

// Function to check whether the input or output file path is valid

void g3_check_input_path(const std::string &path);
void g3_check_output_path(const std::string &path);

#endif

