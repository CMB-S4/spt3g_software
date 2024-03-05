=====
serve
=====

Example web server for .g3 files

The spt3g-json-serve program will serves .g3 files within a directory hierarchy in json format, providing
a crude way to expose a files to a web browser.


Running the server
==================

Command line usage::

  Usage: spt3g-json-serve [OPTION]
  	-h,--help		display this message
  	-P,--Python		initialize python interpreter
  	-p,--port		port to use (default 2726)
  	-b,--bind-address	bind address (default 0.0.0.0)
  	-d,--doc-root		base dir to serve (default is .)
  	-t,--threads		number of threads to use (default 4)
  	-v,--verbose		verbose mode

The example server doesn't support authentication and has not rigorously been
tested for security (though it does try to avoid the most obvious directory
traversal attacks). If it is to be exposed to the internet, you probably should
probably at least reverse proxy it behind a basic authentication screendoor.


Implemented endpoints
====================

*/[?json]
---

List a directory within basedir. By default, this uses a clunky HTML browser,
appending the ``?json`` parameter, the directory listing will be returned
in a json format more suited to be programtically used by a less clunky browser.


*.g3[.gz|.bz2][?N=0&nskip=0]
---

List a .g3 file as JSON. Optional parameters ``N`` and ``nskip`` control the
number of frames (0 means all) and the number of frames from the begining to
skip, respectively.  Note that this just affects what gets transferred over the
wire, the server must still deserialize the file up to nskip + N frames each
time (in the future, caching may be implemented, though the OS file cache probably helps a little bit).

/version
---
List a version number (in JSON format)

/favicon.ico
---
What else, but a favicon?
