/**
 * First pass at a .g3 json server. Using httplib for now.
 *
 * This server is intended to serve a file hierarchy of .g3 files over the web.
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
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/python.hpp>
#include <mutex>
#include <vector>
#include <streambuf>
#include <csignal>

#include "core/G3.h"
#include "core/G3Frame.h"
#include "core/dataio.h"
#include "core/pybindings.h"

#include "httplib/httplib.h"
#include "serve/favicon.h"
namespace fs = boost::filesystem; // from <boost/filesystem.hpp>

static unsigned short port = 2726;
static std::string bind_address = "0.0.0.0";
static std::string basedir = "." ;
static unsigned short nthreads = 4;
static bool use_python = false;
static bool verbose = false;

//mutex for cout synchronization on worker threads
static std::mutex cout_mutex;

static const char * VALID_G3_EXTENSIONS[] = { ".g3", ".g3.gz", ".g3.bz2" };
const int N_VALID_G3_EXTENSIONS = sizeof(VALID_G3_EXTENSIONS)/sizeof(*VALID_G3_EXTENSIONS); 

static void
usage()
{
	std::cout << "Usage: spt3g-json-serve [OPTION]" << std::endl;
	std::cout << "\t-h,--help\t\tdisplay this message" << std::endl;
	std::cout << "\t-P,--Python\t\tinitialize python interpreter" << std::endl;
	std::cout << "\t-p,--port\t\tport to use (default 2726)" << std::endl;
	std::cout << "\t-b,--bind-address\tbind address (default 0.0.0.0)" << std::endl;
	std::cout << "\t-d,--doc-root\t\tbase dir to serve (default is .)" << std::endl;
	std::cout << "\t-t,--threads\t\tnumber of threads to use (default 4) " << std::endl;
	std::cout << "\t-v,--verbose\t\tverbose mode" << std::endl;
}

// pointer to server for stopping
httplib::Server * the_server = nullptr;

static void sig_handler(int signum)
{
	if (the_server) the_server->stop();
}


static int
parse(int nargs, char ** args)
{
	for (int i = 1; i < nargs; i++)
	{
		if (!strcmp(args[i],"-p") || !strcmp(args[i],"--port")) {
			port = atoi(args[++i]);
			if (port <= 0 || port > 65536)
			{
				std::cerr << "Invalid port: " << port << std::endl;
				return 1;
			}
		}
		else if (!strcmp(args[i],"-b") || !strcmp(args[i],"--bind-address")) {
    bind_address = args[++i];
		}
		else if (!strcmp(args[i],"-d") || !strcmp(args[i],"--doc-root")) {
			basedir = args[++i];
		}
		else if (!strcmp(args[i],"-t") || !strcmp(args[i],"--threads")) {
			int req_threads = atoi(args[++i]);
			if (req_threads > 0)
				nthreads = req_threads;
 			else
				std::cerr << "Ignoring invalid number of threads" << std::endl;
		}
		else if (!strcmp(args[i],"-h") || !strcmp(args[i],"--help")) {
			usage();
			return 1;
		}
		else if (!strcmp(args[i],"-P") || !strcmp(args[i],"--Python")) {
			use_python = true;
		}
		else if (!strcmp(args[i],"-v") || !strcmp(args[i],"--version")) {
			verbose = true;
		}
		else {
			std::cerr << "Unknown argument: " << args[i] << std::endl;
			usage();
			return 1;
		}
	}

	return 0;
}


static void
invalid(const httplib::Request & req, httplib::Response & resp, bool html = false)
{
	resp.status = 500;
	if (html) {
		std::string ans = " <html><head><title>Invalid Path</title></head><body>";
		ans+= req.path;
		ans+= " is invalid.</body></html>";
		resp.set_content(ans, "text/html");
	}
 	else {
		std::string ans = "{ \"error\" : \"" + req.path + " is invalid\" }";
		resp.set_content(ans,"application/json");
  }
}



static void notfound(const httplib::Request & req, httplib::Response & resp, bool html = false)
{
	resp.status = 404;
	if (html) {
		std::string ans = " <html><head><title>Path not found</title></head><body>";
		ans+= req.path;
		ans+= " not oound.</body></html>";
		resp.set_content(ans, "text/html");
	}
	else {
		std::string ans = "{ \"error\" : \"" + req.path + " not found\" }";
		resp.set_content(ans,"application/json");
	}
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
			if (c == EOF) {
				flush();
				return EOF;
			}
			else {
				buf_[filled_++] = (char) c;
				if (filled_ == buf_.size()) {
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


// handle a g3 file, producing json
// two ways to modify behavior
// ?N=nframes  to limit number of frames prinited (default is all, equivalent to negative)
// ?skip=nskip  to skip the first nskip frames (default is 0) 
static void handle_g3(const httplib::Request & req, httplib::Response & resp)
{
	auto t0= std::chrono::high_resolution_clock::now();
	if (!httplib::detail::is_valid_path(req.path)) {
		if (verbose) {
			std::cout << "[ " << req.path << "::" << "file invalid ]"  << std::endl;
    }
		invalid(req,resp);
		return;
	}

	//check if file is there
	fs::path p(basedir + req.path);
	if (!fs::exists(p)) {
		if (verbose) {
			std::cout << "[ " << p << "::" << "file not found ]"  << std::endl;
    }
		notfound(req,resp);
		return;
	}

	int n_frames_to_read = 0;
	int n_skip = 0;

	//check if we limit the number of frames
	if (req.has_param("N")) {
		try {
			n_frames_to_read = std::stoi(req.get_param_value("N"));
		}
		catch(...) {
			resp.status = 500;
			resp.set_content("{ \"error\" : \"invalid N\" }", "application/json");
			return;
		}
	}

	//check if we skip some frames
	if (req.has_param("skip")) {
		try {
			n_skip = std::stoi(req.get_param_value("skip"));
		}
		catch(...) {
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
			try {
				//write a preamble:
				std::string preamble = "{\n\"filename\" : \"";
				preamble += req.path;
				preamble += "\",\n\"frames\": [";
				os << preamble;

				//set up the frame
				G3Frame frame;
				std::string fullpath = basedir + req.path;
				boost::iostreams::filtering_istream ifs;
				g3_istream_from_path(ifs, fullpath, 0.01);

				bool first = true;
				int i = 0;
				while(ifs.peek() != EOF) {
					if (!first)
						os << ","  << std::endl;

					frame.load(ifs);
					if (i++ >= n_skip) {
						frame.saveJSON(os);
						first = false;
					}
					if (n_frames_to_read > 0 && i >= n_frames_to_read + n_skip)
						break;
				}
				std::string end = "]\n}\n" ;
				os << end;
				sinkstr.flush();
				sink.done();
				return true;
			}
			catch(...) {
				sinkstr.flush();
				std::string err = "],\n\"error\" : \"reading error\" }";
				sink.write(err.c_str(),err.length());
				sink.done();
				return false;
			}
		}
		,
		[t0,&req, n_frames_to_read, n_skip](bool success) // prints out request and how long it took
		{
			if (verbose) {
				std::chrono::duration<double> diff = std::chrono::high_resolution_clock::now()-t0;
				std::lock_guard<std::mutex> cout_guard(cout_mutex);
				std::cout<< "["<< req.remote_addr << "]" <<req.path << "(N=" << n_frames_to_read << ",n_skip="
          << n_skip << ")::" << (success ? "success: " : "fail: ")<< diff.count() << " s" << std::endl;
			}
		}
	);

}


template<typename T> static void
append_str_array(std::string & s, std::vector<T> & v)
{
	std::stringstream ss("[");
	for (size_t i = 0; i < v.size(); i++) {
		if (!std::is_arithmetic<T>::value) ss << "\"";
		ss << v[i];
		if (!std::is_arithmetic<T>::value) ss << "\"";
		if (i < v.size()-1) ss << ",";
	}
	ss << "]";
	s += ss.str();
}

static bool
valid_g3_ext(const std::string & filename)
{
	bool valid = false;
	for (const char * ext : VALID_G3_EXTENSIONS) {
		if (boost::algorithm::ends_with(filename, ext))
			valid = true;
	}
	return valid;
}



// handle a directory
// a directory can elicit either an html or json response
// we'll make the html response default, and the json response
// available with the ?json parameter
static void
handle_dir(const httplib::Request & req, httplib::Response & resp)
{

  //check if we should output json or html, default to html

  bool output_html = true;
  if (req.has_param("json"))
    output_html = false;

	if (!httplib::detail::is_valid_path(req.path)) {
		if (verbose) {
			std::cout << "[ " << req.path << "::" << "dir invalid ]"  << std::endl;
		}
		invalid(req,resp, output_html);
		return;
	}

	fs::path dir(basedir + req.path);
	if (!fs::is_directory(dir)) {
		std::cout << "[ " << req.path << "::" << "dir not found ]"  << std::endl;
		notfound(req,resp, output_html);
		return;
	}

	if (verbose) {
		std::cout << "[ " << req.path << "::" << "listing ]"  << std::endl;
	}
	std::vector<std::string> dirs;
	std::vector<std::string> files;
	std::vector<size_t> file_sizes;

	for (auto x : fs::directory_iterator(dir)) {
		//skip hidden files
		if (x.path().string()[0]== '.' && (x.path().string().length() < 2 || x.path().string()[1] != '/'))
			continue;

		if (fs::is_directory(x.path()))
		{
			dirs.push_back(x.path().string().substr(basedir.length()));
		}

		else if (valid_g3_ext(x.path().string())) {
			files.push_back(x.path().string().substr(basedir.length()));
			file_sizes.push_back( fs::file_size(x.path()));
		}
	}

	//json response
  if (!output_html)
	{
		std::string json_response = "{\n\"files\" : ";
		append_str_array(json_response,files);
		json_response += ",\n\"dirs\" : ";
		append_str_array(json_response,dirs);
		json_response += ",\n\"file_sizes\" : ";
		append_str_array(json_response,file_sizes);
		json_response += "\n}";
		resp.set_content(json_response,"application/json");
		return;
	}

	//html response
	std::stringstream html_response;
	html_response << "<html>\n<head>\n<title>Listing of ";
	html_response << req.path << "</title>";
	html_response << R"js(
		<script type="text/javascript">

	  function show_g3(what)
    {
			var req = new XMLHttpRequest();
			req.onreadystatechange = function() {
				if (this.readyState == 4 && this.status == 200) {
					document.getElementById("result").innerHTML = this.responseText;
				}
			}

			var skip = document.getElementById('skip').value;
			var N = document.getElementById('N').value;
			req.open("GET", what+'?N='+N+'&skip='+skip);
      req.send();
    }

    </script>
  )js";

	html_response << R"css(
  <style type="text/css">

    .float-container {
			border: 3px;
			padding: 10 px;
		}

		.ls {
			width: 25%;
			float: left;
			padding: 10px;
		}

		.content {
			width: 50%;
			float: left;
			padding: 10px;
		}


	</style>
	)css";

	html_response << "</head>\n<body><h2>Listing of ";
	html_response << req.path;
	html_response << "</h2><hr>\n";

	html_response << "Number of requested frames (0 for all): "
									 "<input type='number' value='0' id='N' min='0'><br>\n";
	html_response << "Number of frames to skip : "
									 "<input type='number' value='0' id='skip' min='0'><br>\n";
	html_response << "<hr><div class='float-container'>\n;";
	html_response << "<div class='ls'>";

	// list dirs first
	// are we top? if not, put a link to ..
	if (req.path!="/")
		html_response << "<b><a href='..'>..</a></b><br>\n";

	for (std::string dir : dirs) {
		html_response << "<b><a href='" << dir  << "/'>"  << dir << "/</a></b><br>\n";
	}

  html_response.precision(4);
	//then loop over files, listing the size also
	for (size_t i = 0; i < files.size(); i++) {
		double sz = file_sizes[i];
		const char  size_suffixes[] = { 'B','K','M','G','T' };
		int suffix_index = 0;
		while (sz > 1024 && suffix_index < sizeof(size_suffixes) - 1) {
			sz /= 1024;
			suffix_index++;
		}

		html_response << "<a href='javascript:show_g3(\"" << files[i]  << "\");'>"
									<< files[i] << "</a> [ " << sz << size_suffixes[suffix_index]
									<< "] [<a href='"  << files[i] <<"'>json</a>]<br>\n";
	}

	html_response << "</div>";

  html_response << R"content(
  <div class='content'>
	<pre id='result'>Result will appear here</pre>
	</div>
  </div>
  </body></html>
  )content";

	resp.set_content(html_response.str() , "text/html");
}


int main(int nargs, char **argv)
{
	if (parse(nargs,argv)) {
		return 1;
	}

	//This USED to be needed for "random python things" stored in frames to work
  // (though it was buggy and didn't play well with multiple threads). Since Sasha's changes in PR 147,
  // you shouldn't need the python interpreter anymore.
	std::unique_ptr<G3PythonInterpreter> interp = use_python ?
		std::make_unique<G3PythonInterpreter>(false) : nullptr;

	signal(SIGPIPE, SIG_IGN);
	signal(SIGINT, sig_handler);

	httplib::Server server;

	//set number of threads
	server.new_task_queue = [] { return new httplib::ThreadPool(nthreads); };

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

	//handler for directory listing
	server.Get(".*/", handle_dir);


	//handler for something ending with .g3
	for (const char * ext : VALID_G3_EXTENSIONS) {
		// regex escape the extension
		std::string escaped_ext = ext;
		boost::replace_all(escaped_ext,".", "\\.");
 	  std::string target = "/.*" + escaped_ext;
		server.Get(target, handle_g3);
	}

	// keep alive seems to delay sending response with content provider
	// TODO: figure out what the right thing to do is here
	server.set_keep_alive_timeout(0);

	the_server = &server;

	std::cout << "spt3g-json-serve is serving " << basedir << " on " << bind_address << ":" << port << std::endl;

	if (!server.listen(bind_address.c_str(), port)) {
		std::cerr << "Could not bind to " << bind_address << ":" << port << std::endl;
		return 1;
	}


	return 0;
}

