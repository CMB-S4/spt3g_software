/**
 * First pass at a .3g json server. Using httplib for now.
 *
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu>
 *
 */


/* enable compression, hopefully? */
#define CPPHTTPLIB_ZLIB_SUPPORT
#define CPPHTTPLIB_USE_POLL
#define VERSION "0.0.0"

#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/python.hpp>
#include <mutex> 
#include <vector>
#include <streambuf>
#include <csignal> 

#include "httplib/httplib.h"
#include "core/G3.h"
#include "core/G3Frame.h"
#include "core/dataio.h"
#include "core/pybindings.h"
#include "serve/favicon.h" 
namespace fs = boost::filesystem; // from <boost/filesystem.hpp>

static unsigned short port = 9000;
static std::string bind_address = "0.0.0.0";
static std::string basedir = "." ;
static unsigned short nthreads = 4;

//mutex for cout synchronization on worker threads
static std::mutex cout_mutex; 

static void
usage()
{
	std::cout << "Usage: spt3g-json-serve [OPTION]" << std::endl;
	std::cout << "  -h,--help  display this message" << std::endl;
	std::cout << "  -p,--port  port to use (default 8080)" << std::endl;
	std::cout << "  -b,--bind-address  bind address (default 0.0.0.0)" << std::endl;
	std::cout << "  -d,--doc-root base dir to serve (default is .)" << std::endl;
	std::cout << "  -t,--threads number of threads to use (default 4) " << std::endl;
}


static int
parse(int nargs, char ** args)
{
	for (int i = 1; i < nargs; i++)
	{
		if (!strcmp(args[i],"-p") || !strcmp(args[i],"--port"))
		{
			port = atoi(args[++i]);
			if (port <= 0 || port > 65536)
			{
				std::cerr << "Invalid port: " << port << std::endl;
				return 1;
			}
		}
		else if (!strcmp(args[i],"-b") || !strcmp(args[i],"--bind-address"))
		{
			bind_address = args[++i];
		}
		else if (!strcmp(args[i],"-d") || !strcmp(args[i],"--doc-root"))
		{
			basedir = args[++i];
		}
		else if (!strcmp(args[i],"-t") || !strcmp(args[i],"--threads"))
		{
			int req_threads = atoi(args[++i]);
			if (req_threads > 0)
				nthreads = req_threads; 
 			else
				std::cerr << "Ignoring invalid number of threads" << std::endl;
		}
		else if (!strcmp(args[i],"-h") || !strcmp(args[i],"--help"))
		{
			usage();
			return 1;
		}
		else
		{
			std::cerr << "Unknown argument: " << args[i] << std::endl;
			usage();
			return 1;
		}
	}

	return 0;
}

static void
invalid(const httplib::Request & req, httplib::Response & resp)
{
	resp.status = 500;
  	std::string ans = "{ \"error\" : \"" + req.path + " is invalid\" }";
  	resp.set_content(ans,"application/json");
}



static void notfound(const httplib::Request & req, httplib::Response & resp)
{
	resp.status = 404;
	std::string ans = "{ \"error\" : \"" + req.path + " not found\" }";
	resp.set_content(ans,"application/json");
}


/* Something is funny about sink.os (because it doesn't override overflow?),
 * so let's just make our own buffered streamer */
class SinkStream : public std::streambuf 
{
	public:
		SinkStream(httplib::DataSink & sink, int bufsize = 1024)
		: sink_(sink), buf_(bufsize), filled_(0)
		{
			;
		}

		bool flush()
		{
			bool success = sink_.write(&buf_[0], filled_);
			filled_ = 0;
			return success;
		}

		int overflow(int c) override
		{
			if (c == EOF)
			{
				flush();
				return EOF;
			}
			else
			{
				buf_[filled_++] = (char) c;
				if (filled_ == buf_.size())
				{
					return flush();
				}
				return 0;
			}
		}
	private:
		httplib::DataSink & sink_;
		std::vector<char> buf_;
		size_t filled_;
};


// handle a g3 file 
static void handle_g3(const httplib::Request & req, httplib::Response & resp)
{
	auto t0= std::chrono::high_resolution_clock::now();
	if (!httplib::detail::is_valid_path(req.path))
	{
		invalid(req,resp);
		return;
	}

	//check if file is there
	fs::path p(basedir + req.path);
	if (!fs::exists(p))
	{
		notfound(req,resp);
		return;
	}

	int n_frames_to_read = -1;
	int n_skip = 0;

	//check if we limit the number of frames
	if (req.has_param("N"))
	{
		try
		{
			n_frames_to_read = std::stoi(req.get_param_value("N"));
		}
		catch(...)
		{
			resp.status = 500;
			resp.set_content("{ \"error\" : \"invalid N\" }", "application/json");
			return;
		}
	}

	//check if we skip some frames
	if (req.has_param("skip"))
	{
		try
		{
			n_skip = std::stoi(req.get_param_value("skip"));
		}
		catch(...)
		{
			resp.status = 500;
			resp.set_content("{ \"error\" : \"invalid skip\" }", "application/json");
			return;
		}
	}

	//TODO investigate chunked content here
	resp.set_content_provider("application/json",
		[n_frames_to_read, n_skip, &req] ( size_t offset, httplib::DataSink & sink )
		{
			SinkStream sinkstr(sink);
			std::ostream os (&sinkstr);
			try
			{
				//write a preamble:
				std::string preamble = "{\n  \"filename\" : \"";
				preamble += req.path;
				preamble += "\",\n  \"frames\": [";
				os << preamble;

				//set up the frame 
				G3Frame frame;
				std::string fullpath = basedir + req.path;
				boost::iostreams::filtering_istream ifs;
				g3_istream_from_path(ifs, fullpath, 0.01);

				bool first = true;
				int i = 0;
				while(ifs.peek() != EOF)
				{
					if (!first)
						os << "  ,  "  << std::endl;

					frame.load(ifs);
					if (i++ >= n_skip) 
					{
						frame.saveJSON(os);
						first = false;
					}
					if (n_frames_to_read >= 0 && i >= n_frames_to_read + n_skip) 
						break;
				}
				std::string end = "  ]\n}\n" ;
				os << end;
				sinkstr.flush();
				sink.done();
				return true;
			}
			catch(...)
			{
				sinkstr.flush(); 
				std::string err = "],\n\"error\" : \"reading error\" }";
				sink.write(err.c_str(),err.length());
				sink.done();
				return false;
			}
		}
		,
		[t0,&req](bool success) // prints out request and how long it took 
		{
			std::chrono::duration<double> diff = std::chrono::high_resolution_clock::now()-t0;
      std::lock_guard<std::mutex> cout_guard(cout_mutex); 
			std::cout<< "["<< req.remote_addr << "]" <<  req.path << "::" <<  
				(success ? "success: " : "fail: ")  << diff.count() << " s" << std::endl;
		}
	);

}


// handle a directory
static void
handle_dir(const httplib::Request & req, httplib::Response & resp)
{
	if (!httplib::detail::is_valid_path(req.path))
	{
		invalid(req,resp);
		return;
	}

	fs::path dir(basedir + req.path);
	if (!fs::is_directory(dir))
	{
		notfound(req,resp);
		return;
	}

	std::stringstream ss_dirs; ss_dirs << "[";
	std::stringstream ss_files;  ss_files << "[";
	std::stringstream ss_sizes;  ss_sizes << "[";

	int ndirs = 0,nfiles = 0;

	for (auto x : fs::directory_iterator(dir))
    {
		//skip hidden stuff
		if (x.path().string()[0]  == '.' && (x.path().string().length() < 2 || x.path().string()[1] != '/'))
		  continue;

		if (fs::is_directory(x.path()))
		{
			if (ndirs)
				ss_dirs << ",";
			ss_dirs << "\"";
			ss_dirs << x.path().string().substr(basedir.length());
			ss_dirs << "/\"";
			ndirs++;
		}
		else if (x.path().extension() == ".g3" || ( (x.path().extension() == ".gz") && (x.path().stem().extension() == ".g3")))
		{
			if (nfiles)
			{
				ss_files << ",";
				ss_sizes << ",";
			}
			ss_files << "\"";
			ss_files << x.path().string().substr(basedir.length());
			ss_sizes << fs::file_size(x.path());
			ss_files << "\"";
			nfiles++;
		}
	}
	ss_dirs << "]";
	ss_files << "]";
	ss_sizes << "]";
	resp.set_content("{\n\"files\" : " + ss_files.str() + ",\n" + "\"dirs\" :" + ss_dirs.str() +
	                 ",\n\"file_sizes\" : " + ss_sizes.str() + "\n}" , "application/json");
}

int main(int nargs, char **argv)
{
	if (parse(nargs,argv))
	{
		return 1;
	}

	//This is needed for "random python things" stored in frames to work
	G3PyContext::initializePython(); 
  G3PyContext::enable(); 
  signal(SIGPIPE, SIG_IGN);

	httplib::Server server;

	//set number of threads
	server.new_task_queue = [] { return new httplib::ThreadPool(nthreads); }  ;

	//this might be useful in the future 
	server.Get("/version", 
		[](const httplib::Request & req, httplib::Response & resp)
		{
			resp.set_content("{ \"spt3g-json-serve-version\" : \"" VERSION "\" } ", "application/json");
		}
	);

	//favicon! 
	server.Get("/favicon.ico", 
		[](const httplib::Request & req, httplib::Response & resp)
		{
			resp.set_content( (const char*) favicon, sizeof(favicon), "image/ico");
		}
	);

	//handler for something ending with .g3
	server.Get("/.*\\.g3", handle_g3);
	server.Get("/.*\\.g3.gz", handle_g3);

	//handler for directory listing
	server.Get(".*/", handle_dir);

    // keep alive seems to delay sending response with content provider
    // TODO: figure out what the right thing to do is here
    server.set_keep_alive_timeout(0);

	if (!server.listen(bind_address.c_str(), port))
	{
		std::cerr << "Could not bind to " << bind_address << ":" << port << std::endl;
		return 1;
	}

  G3PyContext::deinitializePython(); 

	return 0;
}

