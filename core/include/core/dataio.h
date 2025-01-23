#ifndef _G3_DATAIO_H
#define _G3_DATAIO_H

#include <string>
#include <vector>
#include <boost/iostreams/filtering_stream.hpp>

typedef boost::iostreams::filtering_istream g3_istream;
typedef boost::iostreams::filtering_ostream g3_ostream;

/**
 * Configure a filtering stream for G3Frame decompression from a local or remote
 * file source.
 *
 * @param  stream   A reference to the filtering istream that will be configured
 *                  by this function.  Must be instantiated prior to this
 *                  function.
 * @param  path     A valid filename on disk, or a TCP socket address.  If a
 *                  filename, the compression scheme is determined from the file
 *                  extension.  Supported compression schemes are gzip or bzip2.
 *                  If a socket address may be in one of two forms:
 *                  "tcp://host:port" to connect to a host on a specific port
 *                  and read until EOF, or use an asterisk instead of a hostname
 *                  to listen on a specific port for the first connection and
 *                  read until EOF.
 * @param  timeout  Timeout in seconds for socket connections.
 * @param  buffersize Advisory buffer size in bytes for aggregating reads
 * @param  counter  If true, add a counter to the stream configuration,
 *                  for use by the g3_istream_tell function.
 * @return File descriptor for socket connections, or -1 for file input.
 */
int g3_istream_from_path(g3_istream &stream, const std::string &path,
    float timeout=-1.0, size_t buffersize=1024*1024, bool counter=false);

/**
 * Seek to a byte offset in an open input file stream.
 *
 * @param  stream   A reference to the filtering istream, e.g. as configured by
 *                  g3_istream_from_path.
 * @param  offset   File offset to seek to, in bytes, relative to the beginning
 *                  of the file.
 * @return New read head position, or -1 on error.
 */
off_t g3_istream_seek(g3_istream &stream, off_t offset);

/**
 * Return the current read head position in an open input file stream.
 *
 * @param  stream   A reference to the filtering istream, e.g. as configured by
 *                  g3_istream_from_path.
 * @return Current read head position, or -1 on error.
 */
off_t g3_istream_tell(g3_istream &stream);

/**
 * Configure a filtering stream for G3Frame decompression from a memory buffer.
 *
 * @param  stream   A reference to the filtering istream that will be configured
 *                  by this function.  Must be instantiated prior to this
 *                  function.
 * @param  buffer   A pointer to a char buffer in memory.
 * @param  len      Size of the buffer in bytes.
 */
void g3_istream_from_buffer(g3_istream &stream, const char *buf, size_t len);

/**
 * Configure a filtering stream for G3Frame compression to a local file.
 *
 * @param  stream   A reference to the filtering ostream that will be configured
 *                  by this function.  Must be instantiated prior to this
 *                  function.
 * @param  path     A valid filename on disk.  If a filename, the compression
 *                  scheme is determined from the file extension.  Supported
 *                  compression schemes are gzip or bzip2.
 * @param  append   If true, append to an existing file on disk.  Otherwise,
 *                  Create a new file or overwrite an existing file.
 * @param  counter  If true, add a counter filter to the stream configuration,
 *                  for use by the g3_ostream_count function.
 */
void g3_ostream_to_path(g3_ostream &stream, const std::string &path,
    bool append=false, bool counter=false);

/**
 * Count the number of bytes written to the output file stream.
 *
 * @param  stream   A reference to the filtering ostream, as configured by
 *                  g3_ostream_to_path with the counter argument set to true.
 * @return Number of bytes written to disk.
 */
size_t g3_ostream_count(g3_ostream &stream);

/**
 * Configure a filtering stream for G3Frame compression to a memory buffer.
 *
 * @param  stream   A reference to the filtering ostream that will be configured
 *                  by this function.  Must be instantiated prior to this
 *                  function.
 * @param  buffer   A reference a char buffer in memory.
 */
void g3_ostream_to_buffer(g3_ostream &stream, std::vector<char> &buf);

/**
 * Check that the input filename is a valid filename on disk.
 *
 * @throws runtime_error  If filename is invalid or missing.
 */
void g3_check_input_path(const std::string &path);

/**
 * Check that the output filename is a valid filename on disk.
 *
 * @throws runtime_error  If filename is empty, or its parent directory is
 *                        missing.
 */
void g3_check_output_path(const std::string &path);

#endif

