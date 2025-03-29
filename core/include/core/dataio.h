#ifndef _G3_DATAIO_H
#define _G3_DATAIO_H

#include <string>
#include <iostream>

/**
 * Configure a filtering stream for G3Frame decompression from a local or remote
 * file source.
 *
 * @param  stream   A reference to the input stream to be configured by this
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
 * @param  ext      Required file extension (excluding compression suffixes)
 */
void
g3_istream_from_path(std::istream &stream, const std::string &path,
    float timeout=-1.0, size_t buffersize=1024*1024, const std::string &ext=".g3");

/**
 * Return the file descriptor handle for socket connections.
 *
 * @param  stream   A reference to the input stream, as configured by
 *                  g3_istream_from_path.
 * @return fd       The socket file descriptor.
 */
int g3_istream_handle(std::istream &stream);

/**
 * Configure a filtering stream for G3Frame compression to a local file.
 *
 * @param  stream   A reference to the output stream to be configured by this
 *                  function.
 * @param  path     A valid filename on disk.  If a filename, the compression
 *                  scheme is determined from the file extension.  Supported
 *                  compression schemes are gzip or bzip2.
 * @param  append   If true, append to an existing file on disk.  Otherwise,
 *                  Create a new file or overwrite an existing file.
 * @param  ext      Required file extension (excluding compression suffixes)
 */
void
g3_ostream_to_path(std::ostream &stream, const std::string &path, bool append=false,
    size_t buffersize=1024*1024, const std::string &ext=".g3");

#endif

