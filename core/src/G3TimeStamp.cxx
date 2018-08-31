#include <pybindings.h>
#include <serialization.h>

#include <stdint.h>
#include <time.h>
#include <sstream>
#include <iomanip>

#include <G3TimeStamp.h>
#include <G3Units.h>


std::string G3Time::Description() const
{
	struct tm tm;
	std::ostringstream desc;
	char ftime[255];
	time_t t;

	t = time_t(time / G3Units::s);
	gmtime_r(&t, &tm);

	strftime(ftime, sizeof(ftime), "%d-%b-%Y:%H:%M:%S", &tm);

	// Funky casts are to avoid FP errors with non-representability of .1
	desc << ftime << "." << std::setfill('0') << std::setw(9) <<
		uint64_t(time % uint64_t(G3Units::s))*uint64_t(1./G3Units::ns);

	return desc.str();
}

std::string G3Time::isoformat() const
{
	struct tm tm;
	std::ostringstream desc;
	char ftime[255];
	time_t t;

	t = time_t(time / G3Units::s);
	gmtime_r(&t, &tm);

	strftime(ftime, sizeof(ftime), "%Y-%m-%dT%H:%M:%S", &tm);

	// Funky casts are to avoid FP errors with non-representability of .1
	desc << ftime << "." << std::setfill('0') << std::setw(9) <<
		uint64_t(time % uint64_t(G3Units::s))*uint64_t(1./G3Units::ns);

	return desc.str();
}

template <class A> void G3Time::serialize(A &ar, const unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("timestamp", time);
}

G3_SERIALIZABLE_CODE(G3Time);

G3Time::G3Time(int y, int d, int h, int m, int s, int ss)
{
	struct tm tm;
	tm.tm_year = y + 100 /* tm_year starts in 1900 */;
	tm.tm_yday = d;
	tm.tm_hour = h;
	tm.tm_min = m;
	tm.tm_sec = s;

	tm.tm_mon = 0;
	tm.tm_mday = tm.tm_yday;
	int64_t last_code = int64_t(timegm(&tm)) * G3Units::second;
	last_code += (uint64_t)ss;
	time = last_code;
}

G3Time::G3Time(std::string t)
{
	const char *rv;
	struct tm tm;
	G3TimeStamp subsecond = 0;

	// Try a number of conversions in order. If one fails, strptime()
	// return NULL, so try again.

	// 1. 26-Jan-2012:23:29:36
	rv = strptime(t.c_str(), "%d-%b-%Y:%H:%M:%S", &tm);

	// 2. 20120126_232936
	if (rv == NULL)
		rv = strptime(t.c_str(), "%Y%m%d_%H%M%S", &tm);

	// 3. 120126_232936
	if (rv == NULL)
		rv = strptime(t.c_str(), "%y%m%d_%H%M%S", &tm);
	
	// 4. 120126 23:29:36
	if (rv == NULL)
		rv = strptime(t.c_str(), "%y%m%d %H:%M:%S", &tm);

	// 5. 2012-01-26T23:29:36+00:00
	if (rv == NULL) {
		rv = strptime(t.c_str(), "%Y-%m-%dT%H:%M:%S%z", &tm);
		tm.tm_sec -= tm.tm_gmtoff; // timegm() doesn't respect TZ
	}

	// 6. 2012-01-26T23:29:36
	if (rv == NULL)
		rv = strptime(t.c_str(), "%Y-%m-%dT%H:%M:%S", &tm);

	// 7. 2012-01-26 23:29:36+00:00
	if (rv == NULL) {
		rv = strptime(t.c_str(), "%Y-%m-%d %H:%M:%S%z", &tm);
		tm.tm_sec -= tm.tm_gmtoff; // timegm() doesn't respect TZ
	}

	// If there is a string of digits after the field after a period,
	// interpret that as a decimal subsecond
	if (rv != NULL && *rv == '.') {
		char *endptr;
		int i;
		G3TimeStamp prefactor = G3Units::s;
		subsecond = strtol(rv+1, &endptr, 10);

		// Figure out what to divide by to get the units right for the
		// number of digits given. Note that, because the internal
		// storage is as an integer, diving the prefactor by more
		// than the unit time interval results in it being zero.
		// Instead, drop precision on the subsecond field if there
		// were too many digits.
		for (i = 0; i < endptr - (rv+1) && prefactor >= 10; i++)
			prefactor /= 10;
		for (; i < endptr - (rv+1); i++)
			subsecond /= 10;
		subsecond *= prefactor;
	}

	if (rv == NULL)
		log_fatal("Could not convert time string \"%s\"", t.c_str());

	time = G3TimeStamp(timegm(&tm)*G3Units::s) + subsecond;
}

bool G3Time::operator==(const G3Time & other) const
{
	return time == other.time;
}

bool G3Time::operator>(const G3Time & other) const
{
	return time > other.time;
}

bool G3Time::operator<(const G3Time & other) const
{
	return time < other.time;
}

bool G3Time::operator>=(const G3Time & other) const
{
	return time >= other.time;
}

bool G3Time::operator<=(const G3Time & other) const
{
	return time <= other.time;
}

bool G3Time::operator!=(const G3Time & other) const
{
	return time != other.time;
}

std::string G3Time::GetFileFormatString() const
{
	time_t sec = time/G3Units::second;
	char string_buffer[20]; // 16 + padding
	strftime (string_buffer, 19, "%Y%m%d_%H%M%S",
		  gmtime(&sec));
	return std::string(string_buffer);
}

G3Time G3Time::Now()
{
	struct timeval tv;

	gettimeofday(&tv, NULL);

	return G3Time(tv.tv_sec*G3TimeStamp(G3Units::s) +
	    tv.tv_usec*G3TimeStamp(G3Units::us));
}

#define MJD_TO_UNIX_DAYS 40587

double G3Time::GetMJD()
{
	return time / G3Units::day + MJD_TO_UNIX_DAYS;
}

void G3Time::SetMJD(double mjd)
{
	time = (mjd - MJD_TO_UNIX_DAYS) * G3Units::day;
}

G3Time::operator double() const
{
	return double(time);
}

G3Time::operator long() const
{
	return long(time);
}

G3Time
G3Time::operator +(G3TimeStamp delta) const
{
	return G3Time(time + delta);
}

G3Time
G3Time::operator -(G3TimeStamp delta) const
{
	return G3Time(time - delta);
}

G3Time &
G3Time::operator +=(G3TimeStamp delta)
{
	time += delta;
	return *this;
}

G3Time &
G3Time::operator -=(G3TimeStamp delta)
{
	time -= delta;
	return *this;
}

static G3Time
g3time_fadd(G3Time &t, double delta)
{
	return (t + G3TimeStamp(delta));
}

static G3Time
g3time_fsub(G3Time &t, double delta)
{
	return (t - G3TimeStamp(delta));
}

static G3TimePtr
g3time_from_timestamp(boost::python::object obj)
{
	namespace bp = boost::python;
	G3TimeStamp t;

	bp::extract<G3Time> is_a_time(obj);
	if (is_a_time.check())
		return G3TimePtr(new G3Time(is_a_time()));

	// This ends up shadowing the string constructor, so add dispatch
	bp::extract<std::string> is_a_string(obj);
	if (is_a_string.check())
		return G3TimePtr(new G3Time(is_a_string()));

	if (PyFloat_Check(obj.ptr()))
		return G3TimePtr(new G3Time(
		    PyFloat_AsDouble(obj.ptr())));

#if PY_MAJOR_VERSION < 3
	t = PyInt_AsUnsignedLongLongMask(obj.ptr());
#else
	t = PyLong_AsLongLong(obj.ptr());
#endif

	if (PyErr_Occurred() != NULL)
		bp::throw_error_already_set();

	return G3TimePtr(new G3Time(t));
}

PYBINDINGS("core") {
	namespace bp = boost::python;

	EXPORT_FRAMEOBJECT(G3Time, init<>(), "UTC Time")
	    .def(bp::init<int, int , int , int, int, int>("Create a timestamp object from IRIG B code", bp::args("y", "d", "h", "m", "s", "ss")))
	    .def(bp::init<std::string>("Create a time object from a string representation. Supported formats are: YYYYMMDD_HHMMSS, YYMMDD_HHMMSS, YYMMDD HH:MM:SS, DD-Mon-YYYY:HH:MM:SS, YYYY-MM-DDTHH:MM:SS[+TZ] (ISO 8601). All can have a fraction of second field after a dot."))
	    .def("__init__", bp::make_constructor(g3time_from_timestamp, bp::default_call_policies(), (bp::arg("timestamp"))), "Create a G3Time from a numeric timestamp")
	    .def("GetFileFormatString", &G3Time::GetFileFormatString, "Get a string corresponding to how SPTpol and GCP name files for this time")
	    .def("isoformat", &G3Time::isoformat, "Return the ISO 8601 formatted timestamp string")
	    .def("Now", &G3Time::Now, "Return a G3Time object corresponding to the current system time")
	    .staticmethod("Now")
	    .def_readwrite("time", &G3Time::time, "Time relative to the UNIX epoch")
	    .add_property("mjd", &G3Time::GetMJD, &G3Time::SetMJD, "Time in MJD")
	    .def(bp::self == bp::self)
	    .def(bp::self != bp::self)
	    .def(bp::self < bp::self)
	    .def(bp::self <= bp::self)
	    .def(bp::self > bp::self)
	    .def(bp::self >= bp::self)
	    .def("__add__", &G3Time::operator +)
	    .def("__add__", &g3time_fadd)
	    .def("__sub__", &G3Time::operator -)
	    .def("__sub__", &g3time_fsub)
	    .def("__float__", &G3Time::operator double)
	    .def("__int__", &G3Time::operator long)
	;
	register_pointer_conversions<G3Time>();
}

