#include <pybindings.h>
#include <serialization.h>
#include <cmath>
#include <typeinfo>
#include <coordinateutils/CutSkyHealpixMap.h>

#ifdef OPENMP_FOUND
#include <omp.h>
#endif

G3_SERIALIZABLE_CODE(HealpixHitPix);
G3_SERIALIZABLE_CODE(CutSkyHealpixMap);

static void
get_flipped_pixel_and_u_scaling(size_t pix, size_t nside, int is_nested,
    int sym_group, long &out_pix, double &out_u_scale)
{
	g3_assert(sym_group >= 0 && sym_group <= 4);

	double theta_scaling[] = {1, -1, 1, -1};
	double theta_off[] = {0, M_PI, 0, M_PI};
	double phi_off[] = {0, 0, M_PI, M_PI};

	double u_scalings[] = {1,-1,1,-1};

	auto p2a_ptr = pix2ang_ring;
	auto a2p_ptr = ang2pix_ring;

	if (is_nested) {
		p2a_ptr = pix2ang_nest;
		a2p_ptr = ang2pix_nest;
	}

	double theta, phi;
	p2a_ptr(nside, pix, &theta, &phi);

	theta = theta * theta_scaling[sym_group] + theta_off[sym_group];
	phi = fmod(phi + phi_off[sym_group], 2*M_PI);

	if (theta >= M_PI || phi >= 2 * M_PI) log_warn("Theta %lf, Phi %lf\n", theta, phi);
	a2p_ptr(nside, theta, phi, &out_pix);


	if (out_pix > nside * nside * 12) log_warn("Shit son %ld\n", out_pix);

	out_u_scale = u_scalings[sym_group];
}

template <class A>
void HealpixHitPix::serialize(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("is_nested_", is_nested_);
	ar & make_nvp("nside_", nside_);
	ar & make_nvp("ipixmin_", ipixmin_);
	ar & make_nvp("ipixmax_", ipixmax_);
	ar & make_nvp("total_", total_);
	ar & make_nvp("pixinds_", pixinds_);
	ar & make_nvp("coord_ref", coord_ref);
	ar & make_nvp("full_to_cut_inds_", full_to_cut_inds_);
}

std::string
HealpixHitPix::Description() const
{
	std::ostringstream os;

	os << "Nside " << nside_ << " Healpix map with " <<
	    pixinds_.size() << " filled pixels";
	return os.str();
}

template <class A>
void CutSkyHealpixMap::serialize(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3SkyMap", base_class<G3SkyMap>(this));
	ar & make_nvp("hitpix", hitpix);
	ar & make_nvp("nside_", nside_);
	ar & make_nvp("is_nested_", is_nested_);
}

HealpixHitPix::HealpixHitPix(const std::vector<int64_t> & pixinds,
    size_t nside, bool is_nested, MapCoordReference coord_reference)
  : coord_ref(coord_reference), is_nested_(is_nested), nside_(nside),
    map_info_(init_map_info(nside, 1), free_map_info)
{
	pixinds_ = pixinds;
	total_ = pixinds.size()-1;

	//notes
	ipixmin_ = 12 * nside * nside;
	ipixmax_ = 0;
	for (size_t i=0; i < pixinds.size(); i++) {
		if (pixinds[i] < ipixmin_) ipixmin_ = pixinds[i];
		if (pixinds[i] > ipixmax_) ipixmax_ = pixinds[i];
	}

	//full_to_cut_inds_ is the inverse;
	full_to_cut_inds_ = std::vector<int64_t>(
	    ipixmax_ - ipixmin_ + 1, total_);
	for (size_t i=0; i < pixinds.size(); i++)
		full_to_cut_inds_[pixinds[i] -  ipixmin_] = i;
}

HealpixHitPix::HealpixHitPix(size_t nside, bool is_nested, MapCoordReference coord_ref)
  : coord_ref(coord_ref), is_nested_(is_nested), nside_(nside),
    ipixmin_(0), ipixmax_(0), total_(0),
    map_info_(init_map_info(nside, 1), free_map_info)
{
}

HealpixHitPix::HealpixHitPix(const FlatSkyMap &in_map, size_t nside,
    bool is_nested)
  : coord_ref(in_map.coord_ref), is_nested_(is_nested), nside_(nside),
    map_info_(init_map_info(nside, 1), free_map_info)
{
	pixels_from_map(in_map);
}

HealpixHitPix::HealpixHitPix(const CutSkyHealpixMap &in_map, size_t nside,
    bool is_nested)
  : coord_ref(in_map.coord_ref), is_nested_(is_nested), nside_(nside),
    map_info_(init_map_info(nside, 1), free_map_info)
{
	pixels_from_map(in_map);
}

void HealpixHitPix::pixels_from_map(const G3SkyMap & in_map)
{
        g3_assert(in_map.coord_ref != MapCoordReference::Local);

	std::vector<double> alpha(1, 0);
	std::vector<double> delta(1, 0);
	std::vector<int> out_inds(1, 0);
	size_t overflow_index = in_map.xdim() * in_map.ydim();
	size_t n_pix = 12 * nside_ * nside_;

	std::vector<bool> is_hit(n_pix, false);
	ipixmin_ = 0;
	ipixmax_ = 0;
	total_ = 0;

	for (size_t i = 0; i < n_pix; i++) {
		double theta, phi;

		if (is_nested_)
			pix2ang_nest(nside_, i, &theta, &phi);
		else
			pix2ang_ring(nside_, i, &theta, &phi);
		delta[0] = (M_PI/2.0 - theta) * G3Units::rad;
		alpha[0] = (phi) * G3Units::rad;

		out_inds = in_map.angles_to_pixels(alpha, delta);
		if (out_inds[0] != overflow_index) {
			is_hit[i] = true;
			if (total_ == 0)
				ipixmin_ = i;
			ipixmax_ = i;
			total_++;
		}
	}

	if (total_ == 0)
		log_fatal("No pixels defined in this map.");

	pixinds_ = std::vector<int64_t>(total_ + 1);
	pixinds_[total_] = -1;
	full_to_cut_inds_ = std::vector<int64_t>(
	    ipixmax_ - ipixmin_ + 1, total_);
	size_t index_for_pixinds = 0;
	for (size_t i = ipixmin_; i < ipixmax_ + 1; i++) {
		if (is_hit[i]) {
			pixinds_[index_for_pixinds] = i;
			full_to_cut_inds_[i -  ipixmin_] = index_for_pixinds;
			index_for_pixinds++;
		}
	}
}

long
HealpixHitPix::get_fullsky_index(long cutsky_index) const
{
	if (cutsky_index < 0 || cutsky_index > total_)
		return -1;
	return pixinds_[cutsky_index];
}

long
HealpixHitPix::get_cutsky_index(long fullsky_index) const
{
	if (fullsky_index < ipixmin_ || fullsky_index > ipixmax_)
		return total_;
	return full_to_cut_inds_[fullsky_index - ipixmin_];
}

//assumed map to be ring ordering
CutSkyHealpixMap::CutSkyHealpixMap(boost::python::object v,
    size_t full_sky_map_nside, HealpixHitPixPtr hitpix, bool is_weighted,
    G3SkyMap::MapPolType pol_type, G3Timestream::TimestreamUnits u,
    int init_sym_group) :
      G3SkyMap(hitpix->coord_ref, hitpix->total_, 1, is_weighted, u,
       pol_type, true),
      hitpix(hitpix),  nside_(hitpix->get_nside()),
      is_nested_(hitpix->is_nested())
{
	size_t nside = hitpix->get_nside();
	bool is_nested = hitpix->is_nested();
	int sym_group = init_sym_group;
	bool is_u = (pol_type == G3SkyMap::U);

	Py_buffer view;
	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) != -1) {
		if (strcmp(view.format, "d") == 0) {
			for (size_t i = hitpix->ipixmin_;
			    i <= hitpix->ipixmax_; i++) {
				long tmp_full_ind;
				double u_scale;

				get_flipped_pixel_and_u_scaling(i, nside,
				    is_nested, sym_group, tmp_full_ind,
				    u_scale);

				u_scale = is_u ? u_scale : 1.0;

				if (hitpix->is_nested() ) {
					long tmp_ind;
					nest2ring(hitpix->get_nside(), tmp_full_ind, &tmp_ind);
					tmp_full_ind = tmp_ind;
				}

				int64_t local_ind =
				    hitpix->full_to_cut_inds_[
				    i - hitpix->ipixmin_];
				if (local_ind < 0) {
					local_ind = data_.size()-1;
				}
				g3_assert(local_ind < data_.size());
				g3_assert(local_ind >= 0);

				data_[local_ind] =
				    ((double *)view.buf)[tmp_full_ind] * u_scale;
			}
		} else if (strcmp(view.format, "f") == 0) {
			for (size_t i = hitpix->ipixmin_; i <= hitpix->ipixmax_;
			    i++) {
				long tmp_full_ind;
				double u_scale;
				get_flipped_pixel_and_u_scaling(i, nside,
				    is_nested, sym_group, tmp_full_ind,
				    u_scale);
				u_scale = is_u ? u_scale : 1.0;


				if (hitpix->is_nested() ) {
					long tmp_ind;
					nest2ring(hitpix->get_nside(), tmp_full_ind, &tmp_ind);
					tmp_full_ind = tmp_ind;
				}

				int64_t local_ind =
				    hitpix->full_to_cut_inds_[
				    i - hitpix->ipixmin_];
				if (local_ind < 0) {
					local_ind = data_.size()-1;
				}
				g3_assert(local_ind < data_.size());
				g3_assert(local_ind >= 0);

				data_[local_ind] =
				    ((float *)view.buf)[tmp_full_ind] * u_scale;
			}
		} else {
			log_fatal("Unknown type code %s", view.format);
		}

		PyBuffer_Release(&view);
		return;
	}

	throw boost::python::error_already_set();
}

CutSkyHealpixMap::CutSkyHealpixMap(boost::python::object v,
    HealpixHitPixPtr hitpix, bool is_weighted, G3SkyMap::MapPolType pol_type,
    G3Timestream::TimestreamUnits u) :
      G3SkyMap(v, hitpix->coord_ref, is_weighted, u, pol_type),
      hitpix(hitpix),  nside_(hitpix->get_nside()),
      is_nested_(hitpix->is_nested())
{
}

CutSkyHealpixMap::CutSkyHealpixMap(HealpixHitPixPtr hitpix, bool is_weighted,
    G3SkyMap::MapPolType pol_type, G3Timestream::TimestreamUnits u) :
      G3SkyMap(hitpix->coord_ref, hitpix->total_, 1, is_weighted, u,
       pol_type, false),
      hitpix(hitpix),  nside_(hitpix->get_nside()),
      is_nested_(hitpix->is_nested())
{
}

G3SkyMapPtr
CutSkyHealpixMap::Clone(bool copy_data) const
{
	if (copy_data)
		return boost::make_shared<CutSkyHealpixMap>(*this);
	else
		return boost::make_shared<CutSkyHealpixMap>(hitpix,
		    is_weighted, pol_type, units);
}

CutSkyHealpixMap::CutSkyHealpixMap() :
    G3SkyMap(MapCoordReference::Equatorial, 0), nside_(0), is_nested_(0)
{
}

bool CutSkyHealpixMap::IsCompatible(const G3SkyMap & other) const {
	try {
		const CutSkyHealpixMap & hpx = dynamic_cast<const CutSkyHealpixMap&>(other);
		return (G3SkyMap::IsCompatible(other) &&
			(nside_ == hpx.nside_) &&
			(is_nested_ == hpx.is_nested_) &&
			(hitpix && hpx.hitpix) &&
			(hitpix->ipixmin_ == hpx.hitpix->ipixmin_) &&
			(hitpix->ipixmax_ == hpx.hitpix->ipixmax_));
	} catch(const std::bad_cast& e) {
		return false;
	}
}

std::vector<double>
CutSkyHealpixMap::get_fullsky_map()
{

	log_debug("%zu %zu\n", hitpix->pixinds_.size(), data_.size());
	g3_assert(hitpix->pixinds_.size() == data_.size());
	size_t nside = hitpix->get_nside();

	std::vector<double> full_sized_map(12 * nside * nside, 0);

	for (size_t i = 0; i < hitpix->total_; i++) {
		size_t index = hitpix->pixinds_[i];
		g3_assert(index < full_sized_map.size() + 1);
		full_sized_map[index] = data_[i];
	}
	return full_sized_map;
}


long
HealpixHitPix::angle_to_pixel(double alpha, double delta, bool cutsky) const
{
	alpha /= G3Units::rad;
	double theta = (90 * G3Units::deg - delta) / G3Units::rad;

	long outpix = -1;

	if ( std::isnan(theta) || std::isnan(alpha) ) {
		return -1;
	}

	if (is_nested_)
		ang2pix_nest(nside_, theta, alpha, &outpix);
	else
		ang2pix_ring(nside_, theta, alpha, &outpix);

	if (cutsky)
		outpix = get_cutsky_index(outpix);

	return outpix;
}

std::vector<double>
HealpixHitPix::pixel_to_angle(long pixel, bool cutsky) const
{
	if (cutsky)
		pixel = get_fullsky_index(pixel);

	double alpha, delta;
	if (is_nested_)
		pix2ang_nest(nside_, pixel, &delta, &alpha);
	else
		pix2ang_ring(nside_, pixel, &delta, &alpha);

	alpha *= G3Units::rad;
	delta = 90 * G3Units::deg - delta * G3Units::rad;

	std::vector<double> out = {alpha, delta};
	return out;
}

void
HealpixHitPix::get_rebin_angles(long pixel, size_t scale,
    std::vector<double> & alphas, std::vector<double> & deltas, bool cutsky) const
{
	if (nside_ % scale != 0)
		log_fatal("Nside must be a multiple of rebinning scale");

	if (cutsky)
		pixel = get_fullsky_index(pixel);

	if (!is_nested_)
		ring2nest(nside_, pixel, &pixel);

	alphas = std::vector<double>(scale * scale);
	deltas = std::vector<double>(scale * scale);

	size_t nside_rebin = nside_ * scale;
	long pixmin = pixel * scale * scale;
	for (size_t i = 0; i < (scale * scale); i++) {
		long p = pixmin + i;
		double theta, phi;
		pix2ang_nest(nside_rebin, p, &theta, &phi);
		alphas[i] = phi * G3Units::rad;
		deltas[i] = 90 * G3Units::deg - theta * G3Units::rad;
	}
}

void
HealpixHitPix::get_interp_pixels_weights(double alpha, double delta,
    std::vector<long> & pix, std::vector<double> & weight, bool cutsky) const
{
	alpha /= G3Units::rad;
	delta /= G3Units::rad;

	double theta = M_PI/2.0 - delta;
	double w[4];
	long fullpix[4];
	get_interp_weights(map_info_.get(), theta, alpha, fullpix, w);

	if (is_nested_) {
		for (size_t i = 0; i < 4; i++) {
			ring2nest(nside_, fullpix[i], fullpix + i);
		}
	}

	pix = std::vector<long>(4, -1);
	weight = std::vector<double>(4, 0);
	for (size_t i = 0; i < 4; i++) {
		pix[i] = cutsky ? get_cutsky_index(fullpix[i]) : fullpix[i];
		weight[i] = w[i];
	}
}

size_t
CutSkyHealpixMap::angle_to_pixel(double alpha, double delta) const
{
	return hitpix->angle_to_pixel(alpha, delta, true);
}

std::vector<double>
CutSkyHealpixMap::pixel_to_angle(size_t pixel) const
{
	return hitpix->pixel_to_angle(pixel, true);
}

void
CutSkyHealpixMap::get_rebin_angles(long pixel, size_t scale,
    std::vector<double> & alphas, std::vector<double> & deltas) const
{
	hitpix->get_rebin_angles(pixel, scale, alphas, deltas, true);
}

void
CutSkyHealpixMap::get_interp_pixels_weights(double alpha, double delta,
    std::vector<long> & pix, std::vector<double> & weight) const
{
	hitpix->get_interp_pixels_weights(alpha, delta, pix, weight, true);
}


G3SkyMapPtr CutSkyHealpixMap::rebin(size_t scale) const
{
	if (nside_ % scale != 0)
		log_fatal("Map nside must be a multiple of rebinning scale");

	if (scale == 1)
		return Clone(true);

	HealpixHitPixPtr hpx = boost::make_shared<HealpixHitPix>(*this,
	    nside_ / scale, is_nested_);
	CutSkyHealpixMap out(hpx, is_weighted, pol_type, units);
        out.EnsureAllocated();

	size_t scale2 = scale * scale;

#ifdef OPENMP_FOUND
#pragma omp parallel for
#endif
	for (long i = 0; i < out.size(); i++) {
		long ipmin = out.hitpix->get_fullsky_index(i);
		if (!is_nested_)
			ring2nest(out.nside_, ipmin, &ipmin);
		ipmin *= scale2;
		double norm = 0;
		for (size_t j = 0; j < scale2; j++) {
			long ip = ipmin + j;
			if (!is_nested_)
				nest2ring(nside_, ip, &ip);
			ip = hitpix->get_cutsky_index(ip);
			if (ip >= size()) continue;
			out[i] += data_[ip];
			norm += 1.;
		}
		out[i] /= norm;
	}

	return boost::make_shared<CutSkyHealpixMap>(out);
}

std::string
CutSkyHealpixMap::Description() const
{
	std::ostringstream os;

	os << "Nside " << nside_ << " Healpix map";
	return os.str();
}

#define CUT_SKY_MAP_DOCSTR \
	"CutSkyHealpixMap is a G3SkyMap with the extra meta information for a "\
	"Healpix-pixelized sky.\nIn practice it behaves (mostly) like a 1D "\
	"numpy array.\n\n"                                                   \
	"Meta information stored:\n\n"\
	"    get_nside() : returns the healpix resolution parameter\n"\
	"    is_nested() : returns True if the pixels are in nested ordering,\n"\
	"        false if in ring ordering\n"\
	"\n\n"\
	"Map array includes only the pixels that are within the sky patch covered\n"\
	"by the map, as defined by the HealpixHitPix hitpix attribute.  To obtain a\n"\
	"full-sized map for use with standard healpix tools, use get_fullsky_map()\n"\
	"utility function.\n\n"\
	"The other meta information is inheritted from G3SkyMap that lives in core. \n\n" \
	"For reasons (skymap __setitem__ has to handle both 1d and 2d \n"\
	" semantics) the CutSkyHealpixMap has a slightly unintuitive way of \n"\
	" setting values when using a slicing operator.  Instead of being\n"\
	" able to slice directly you need to cast it to be an array first: \n\n" \
	"    np.asarray(your_cut_sky_map)[:] = the_numpy_array_you_are_assigning\n\n\n"\
	"If you find that you need numpy functionality from a CutSkyHealpixMap,\n"\
	" using np.asarray will convert it to a numpy array without copying the data.\n" \
	" any changes to the resulting numpy array will affect the data stored in the\n" \
	" CutSkyHealpixMap."

PYBINDINGS("coordinateutils")
{
	namespace bp = boost::python;

	EXPORT_FRAMEOBJECT(HealpixHitPix, init<>(),
	  "The HealpixHitPix is a class that encapsulates the "
	  "transform between cut sky and full sky healpix maps.  "
	  "It provides the mapping between the indices.\n\n"
	  "Out of necessity this object also defines a patch of "
	  "sky, and therefore the mapping between pixel and sky angle.")
	    .def(bp::init<const FlatSkyMap &, size_t, int>(
		(bp::arg("in_map"), bp::arg("nside"),
		 bp::arg("is_nested"))))
	    .def(bp::init<const CutSkyHealpixMap &, size_t, int>(
		(bp::arg("in_map"), bp::arg("nside"),
		 bp::arg("is_nested"))))
            .def(bp::init<const std::vector<int64_t> &, size_t, int,
                 MapCoordReference>(
	         (bp::arg("pix_inds"), bp::arg("nside"), bp::arg("is_nested"),
		  bp::arg("coord_ref"))))
	    .def_pickle(g3frameobject_picklesuite<HealpixHitPix>())
	    .def("get_size", &HealpixHitPix::get_size, "Number of hit pixels")
	    .def("get_nside", &HealpixHitPix::get_nside, "Map resolution parameter")
	    .def("is_nested", &HealpixHitPix::is_nested,
	        "True if pixels are in nested ordering, False if in ring ordering")
	    .def("get_fullsky_index", &HealpixHitPix::get_fullsky_index,
	         bp::arg("cut_sky_index"))
	    .def("get_cutsky_index", &HealpixHitPix::get_cutsky_index,
		 bp::arg("full_sky_index"))
	    .def_readwrite("pixinds", &HealpixHitPix::pixinds_,
	        "Vector of full sky indices corresponding to each cut sky pixel")
	    ;

	bp::class_<CutSkyHealpixMap, bp::bases<G3SkyMap>, CutSkyHealpixMapPtr>
		("CutSkyHealpixMap", CUT_SKY_MAP_DOCSTR, bp::no_init)
	    .def(bp::init<HealpixHitPixPtr, bool,
		 G3SkyMap::MapPolType, G3Timestream::TimestreamUnits>(
			 (bp::arg("hitpix"),
			  bp::arg("is_weighted") = true,
			  bp::arg("pol_type") = G3SkyMap::None,
			  bp::arg("units") = G3Timestream::Tcmb),
			 "Initialize an empty sky map from sky patch "
			 "defined by the hitpix input"))
	    .def(bp::init<bp::object, size_t, HealpixHitPixPtr, bool,
		 G3SkyMap::MapPolType, G3Timestream::TimestreamUnits, int >(
			 (bp::arg("full_sky_map"),
			  bp::arg("full_sky_smap_nside"),
			  bp::arg("hitpix"),
			  bp::arg("is_weighted") = true,
			  bp::arg("pol_type") = G3SkyMap::None,
			  bp::arg("units") = G3Timestream::Tcmb,
			  bp::arg("init_sym_group") = 0),
			 "Initialize a sky map by copying data from the full sky map "
			 "into the sky patch defined by the hitpix input"))
	    .def(bp::init<bp::object, HealpixHitPixPtr, bool,
		 G3SkyMap::MapPolType, G3Timestream::TimestreamUnits>(
			 (bp::arg("cut_sky_map"), bp::arg("hitpix"),
			  bp::arg("is_weighted") = true,
			  bp::arg("pol_type") = G3SkyMap::None,
			  bp::arg("units") = G3Timestream::Tcmb),
			 "Create a sky map from a numpy array"))
	    .def(boost::python::init<>())
	    .def_pickle(g3frameobject_picklesuite<CutSkyHealpixMap>())
	    .def("get_fullsky_map",
		 &CutSkyHealpixMap::get_fullsky_map,
		 bp::arg("out_fullsky_map"),
		 "Return a full sky array for use with standard healpix tools")
	    .def("pixel_to_angle", &CutSkyHealpixMap::pixel_to_angle,
	         "Compute the sky coordinates of the given pixel")
	    .def("is_nested", &CutSkyHealpixMap::is_nested,
	         "True if pixels are in nested ordering, false if in ring ordering")
	    .def("get_nside", &CutSkyHealpixMap::get_nside,
	         "Map resolution parameter")

	    .def_readwrite("hitpix", &CutSkyHealpixMap::hitpix)
	;

        register_pointer_conversions<CutSkyHealpixMap>();
	bp::implicitly_convertible<CutSkyHealpixMapPtr, G3SkyMapPtr>();
	bp::implicitly_convertible<CutSkyHealpixMapConstPtr,
	    G3SkyMapConstPtr>();
	bp::implicitly_convertible<CutSkyHealpixMapPtr, G3SkyMapConstPtr>();
}

