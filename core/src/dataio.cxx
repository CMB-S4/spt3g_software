#include <G3Logging.h>
#include <dataio.h>

#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#ifdef Boost_BZIP2_FOUND
#include <boost/iostreams/filter/bzip2.hpp>
#endif
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/filesystem.hpp>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdlib.h>

int
g3_istream_from_path(boost::iostreams::filtering_istream &stream,
    const std::string &path, float timeout)
{
	stream.reset();
	if (boost::algorithm::ends_with(path, ".gz"))
		stream.push(boost::iostreams::gzip_decompressor());
	if (boost::algorithm::ends_with(path, ".bz2")) {
#ifdef Boost_BZIP2_FOUND
		stream.push(boost::iostreams::bzip2_decompressor());
#else
		log_fatal("Boost not compiled with bzip2 support.");
#endif
	}

	int fd = -1;

	// Figure out what kind of ultimate data source this is
	if (path.find("tcp://") == 0) {
		// TCP Socket. Two syntaxes:
		// - tcp://host:port -> connect to "host" on "port" and read
		//   until EOF
		// - tcp://*:port -> listen on "port" for the first connection
		//   and read until EOF

		std::string host = path.substr(path.find("://") + 3);
		if (host.find(":") == -1)
			log_fatal("Could not open URL %s: unspecified port",
			    path.c_str());
		std::string port = host.substr(host.find(":") + 1);
		host = host.substr(0, host.find(":"));

		log_debug("Opening connection to %s, port %s", host.c_str(),
		    port.c_str());

		if (strcmp(host.c_str(), "*") == 0) {
			// Listen for incoming connections
			struct sockaddr_in6 sin;
			int no = 0, yes = 1;
			int lfd;

			bzero(&sin, sizeof(sin));
			sin.sin6_family = AF_INET6;
			sin.sin6_port = htons(strtol(port.c_str(), NULL, 10));
		#ifdef SIN6_LEN
			sin.sin6_len = sizeof(sin);
		#endif
			
			lfd = socket(PF_INET6, SOCK_STREAM, 0);
			if (lfd <= 0)
				log_fatal("Could not listen on %s (%s)",
				    path.c_str(), strerror(errno));
			setsockopt(lfd, IPPROTO_IPV6, IPV6_V6ONLY, &no,
			    sizeof(no));
			setsockopt(lfd, SOL_SOCKET, SO_REUSEADDR, &yes,
			    sizeof(yes));

			if (bind(lfd, (struct sockaddr *)&sin, sizeof(sin)) < 0)
				log_fatal("Could not bind on port %s (%s)",
				    port.c_str(), strerror(errno));
			if (listen(lfd, 1) < 0)
				log_fatal("Could not listen on port %s (%s)",
				    port.c_str(), strerror(errno));

			log_debug("Waiting for connection on port %s",
			    port.c_str());
			fd = accept(lfd, NULL, NULL);
			log_debug("Accepted connection on port %s",
			    port.c_str());
			close(lfd);
		} else {
			// Connect to a listening host elsewhere

			struct addrinfo hints, *info, *r;
			int err;

			bzero(&hints, sizeof(hints));
			hints.ai_family = AF_UNSPEC;
			hints.ai_socktype = SOCK_STREAM;

			err = getaddrinfo(host.c_str(), port.c_str(), &hints,
			    &info);
			if (err != 0)
				log_fatal("Could not find host %s (%s)",
				    host.c_str(), gai_strerror(err));

			// Loop through possible addresses until we find one
			// that works.
			fd = -1;
			for (r = info; r != NULL; r = r->ai_next) {
				fd = socket(r->ai_family, r->ai_socktype,
				    r->ai_protocol);
				if (fd == -1)
					continue;

				if (connect(fd, r->ai_addr, r->ai_addrlen) ==
				    -1) {
					close(fd);
					fd = -1;
					continue;
				}

				break;
			}

			if (fd == -1)
				log_fatal("Could not connect to %s (%s)",
				    path.c_str(), strerror(errno));

                        if (timeout >= 0) {
                                struct timeval tv;
                                tv.tv_sec = (int)timeout;
                                tv.tv_usec = (int)(1e6 * (timeout - tv.tv_sec));
                                if (setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO,
                                               (char *)&tv, sizeof(tv)) < 0)
                                        log_fatal("Failed to set timeout on socket; errno=%i",
                                                  errno);
                        }

			if (info != NULL)
				freeaddrinfo(info);
		}

		boost::iostreams::file_descriptor_source fs(fd,
		    boost::iostreams::close_handle);
		stream.push(fs);
	} else {
		// Simple file case
		stream.push(boost::iostreams::file_source(path,
		    std::ios::binary));
	}

	return fd;
}
