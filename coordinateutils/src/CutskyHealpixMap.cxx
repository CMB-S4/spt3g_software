#include <pybindings.h>
#include <serialization.h>
#include <math.h>
#include <coordinateutils/coordinateutils.h>
#include <coordinateutils/CutSkyHealpixMap.h>

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
	phi = phi + phi_off[sym_group];
	a2p_ptr(nside, theta, phi, &out_pix);

	out_u_scale = u_scalings[sym_group];
}

template <class A>
void HealpixHitPix::serialize(A &ar, unsigned u)
{
	using namespace cereal;

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
void CutSkyHealpixMap::serialize(A &ar, unsigned u)
{
        using namespace cereal;

        ar & make_nvp("G3SkyMap", base_class<G3SkyMap>(this));
        ar & make_nvp("hitpix", hitpix);
        ar & make_nvp("nside_", nside_);
        ar & make_nvp("is_nested_", is_nested_);
}

HealpixHitPix::HealpixHitPix(const std::vector<int64_t> & pixinds,
    size_t nside, bool is_nested, MapCoordReference coord_reference)
{
	nside_ = nside;
	is_nested_ = is_nested;
	pixinds_ = pixinds;
	total_ = pixinds.size()-1;
	coord_ref = coord_reference;

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

HealpixHitPix::HealpixHitPix(const FlatSkyMap &flat_map, size_t nside,
    bool is_nested) 
  : coord_ref(flat_map.coord_ref), is_nested_(is_nested), nside_(nside)
{
        g3_assert(flat_map.coord_ref != MapCoordReference::Local);
	
	std::vector<double> alpha(1, 0);
	std::vector<double> delta(1, 0);
	std::vector<int> out_inds(1, 0);
	size_t overflow_index = flat_map.xdim() * flat_map.ydim();
	size_t n_pix = 12 * nside * nside;

	std::vector<bool> is_hit(n_pix, false);
	ipixmin_ = 0;
	ipixmax_ = 0;
	total_ = 0;

	for (size_t i = 0; i < n_pix; i++) {
		double theta, phi;

		if (is_nested)
			pix2ang_nest(nside, i, &theta, &phi);
		else
			pix2ang_ring(nside, i, &theta, &phi);
		delta[0] = (M_PI/2.0 - theta) * G3Units::rad;
		alpha[0] = (phi) * G3Units::rad;

		out_inds = flat_map.angles_to_pixels(alpha, delta);
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

//assumed to be ring ordering
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

CutSkyHealpixMap::CutSkyHealpixMap() :
    G3SkyMap(MapCoordReference::Equatorial, 0), nside_(0), is_nested_(0)
{
}

std::vector<double>
CutSkyHealpixMap::get_full_sized_healpix_map()
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


std::vector<int> CutSkyHealpixMap::angles_to_pixels(
    const std::vector<double> & alphas, 
    const std::vector<double> & deltas) const { 
	std::vector<int> v( alphas.size(), 0);
	for (size_t i=0; i < alphas.size(); i++){
		v[i] = ang_2_pix_(alphas[i], deltas[i]);
	}
	return v;
}


long
CutSkyHealpixMap::ang_2_pix_(double alpha, double delta) const {
	alpha /= G3Units::rad;
	delta /= G3Units::rad;

	double theta = M_PI/2.0 - delta;
	long outpix = 0;

	if ( isnan(theta) || isnan(alpha) ) {
		return -1;
	}

	if (is_nested_)
		ang2pix_nest(nside_, theta, alpha, &outpix);
	else
		ang2pix_ring(nside_, theta, alpha, &outpix);
	return hitpix->get_cutsky_index(outpix);
}

std::vector<double>
CutSkyHealpixMap::pixel_to_angle(size_t pix) const {
	double alpha=0;
	double delta=0;

	pix = hitpix->get_fullsky_index(pix);

	double theta, phi;
	if (is_nested_)
		pix2ang_nest(nside_, pix, &theta, &phi);
	else
		pix2ang_ring(nside_, pix, &theta, &phi);

	alpha = phi * G3Units::rad;
	delta = M_PI/2.0 * G3Units::rad - theta * G3Units::rad;
	std::vector<double> ov(2,0); 
	ov.push_back(alpha);
	ov.push_back(delta);
	return ov;
}

void
CutSkyHealpixMap::get_interpolated_weights(double alpha, double delta,
    long pix[4], double weight[4]) const
{
	WCSMapInfoPtr ugh(init_map_info(nside_, 0), free_map_info);
	alpha /= G3Units::rad;
	delta /= G3Units::rad;

	double theta = M_PI/2.0 - delta;
	get_interp_weights(ugh.get(), theta, alpha, pix, weight);
}

double
CutSkyHealpixMap::get_interp_precalc(long pix[4], double weight[4]) const
{
	double outval = 0;
	for (size_t i = 0; i < 4; i++) {
		size_t cutsky_pix = hitpix->get_cutsky_index(pix[i]);
		outval += data_[cutsky_pix] * weight[i];
	}
	return outval;
}

double
CutSkyHealpixMap::get_interpolated_value(double alpha, double delta) const
{
	WCSMapInfoPtr ugh(init_map_info(nside_, 0), free_map_info);

	alpha /= G3Units::rad;
	delta /= G3Units::rad;

	double theta = M_PI/2.0 - delta;

	long pix[4];
	double weight[4];
	get_interp_weights(ugh.get(), theta, alpha, pix, weight);

	double outval = 0;
	for (size_t i = 0; i < 4; i++) {
		size_t cutsky_pix = hitpix->get_cutsky_index(pix[i]);
		outval += data_[cutsky_pix] * weight[i];
	}

	return outval;
}

std::string
CutSkyHealpixMap::Description() const
{
	std::ostringstream os;

	os << "Nside " << nside_ << " Healpix map";
	return os.str();
}

PYBINDINGS("coordinateutils")
{
	namespace bp = boost::python;

	EXPORT_FRAMEOBJECT(HealpixHitPix, init<>(),
	  "The HealpixHitPix is a class that encapsulates the "
	  "transform between cut sky and full sky healpix maps.  "
	  "It provides the mapping between the indices.\n\n"
	  "Out of necessity this object also defines a patch of "
	  "sky.")
	    .def(bp::init<const FlatSkyMap &, size_t, int>(
		(bp::arg("flat_map"), bp::arg("nside"),
		 bp::arg("is_nested"))))
            .def(bp::init<const std::vector<int64_t> &, size_t, int,
                 MapCoordReference>(
	         (bp::arg("pix_inds"), bp::arg("nside"), bp::arg("is_nested"),
		  bp::arg("coord_ref"))))
	    .def_pickle(g3frameobject_picklesuite<HealpixHitPix>())
	    .def("get_nside", &HealpixHitPix::get_nside)
	    .def("is_nested", &HealpixHitPix::is_nested)
	    .def("get_fullsky_index", &HealpixHitPix::get_fullsky_index,
	         bp::arg("cut_sky_index"))
	    .def("get_cutsky_index", &HealpixHitPix::get_cutsky_index,
		 bp::arg("full_sky_index"))
	    .def_readwrite("pixinds", &HealpixHitPix::pixinds_)
	    ;

	// XXX docstrings
	bp::class_<CutSkyHealpixMap, bp::bases<G3SkyMap>, CutSkyHealpixMapPtr>
		("CutSkyHealpixMap", bp::no_init)
	    .def(bp::init<HealpixHitPixPtr, bool,
		 G3SkyMap::MapPolType, G3Timestream::TimestreamUnits>(
			 (bp::arg("hitpix"),
			  bp::arg("is_weighted"),
			  bp::arg("pol_type"),
			  bp::arg("units"))
			 ))

	    .def(bp::init<bp::object, size_t, HealpixHitPixPtr, bool,
		 G3SkyMap::MapPolType, G3Timestream::TimestreamUnits, int >(
			 (bp::arg("full_sky_map"),
			  bp::arg("full_sky_smap_nside"),
			  bp::arg("hitpix"),
			  bp::arg("is_weighted") = true,
			  bp::arg("pol_type") = G3SkyMap::None,
			  bp::arg("units") = G3Timestream::Kcmb,
			  bp::arg("init_sym_group") = 0)))
		.def(bp::init<bp::object, HealpixHitPixPtr, bool,
		     G3SkyMap::MapPolType, G3Timestream::TimestreamUnits>(
			     (bp::arg("cut_sky_map"), bp::arg("hitpix"),
			      bp::arg("is_weighted") = true,
			      bp::arg("pol_type") = G3SkyMap::None,
			      bp::arg("units") = G3Timestream::Kcmb))
			)
		.def(boost::python::init<>())
		.def_pickle(g3frameobject_picklesuite<CutSkyHealpixMap>())
		.def("get_full_sized_healpix_map",
		     &CutSkyHealpixMap::get_full_sized_healpix_map,
		     bp::arg("out_full_sized_map"))
		.def("pixel_to_angle", &CutSkyHealpixMap::pixel_to_angle)
		.def("is_nested", &CutSkyHealpixMap::is_nested)
		.def("get_nside", &CutSkyHealpixMap::get_nside)
		
	    .def_readwrite("hitpix", &CutSkyHealpixMap::hitpix)
	;

        register_pointer_conversions<CutSkyHealpixMap>();
	bp::implicitly_convertible<CutSkyHealpixMapPtr, G3SkyMapPtr>();
	bp::implicitly_convertible<CutSkyHealpixMapConstPtr,
	    G3SkyMapConstPtr>();
	bp::implicitly_convertible<CutSkyHealpixMapPtr, G3SkyMapConstPtr>();
}

