from spt3g.core import G3Module, MapCoordReference, G3FrameType, G3VectorDouble
from spt3g.core import G3Units as units
from spt3g.coordinateutils import (coord_eq_to_gal, coord_gal_to_eq,
                                   pol_angle_eq_to_gal, pol_angle_gal_to_eq)

class ChangeCoordSys(G3Module):
    def __init__(self, coord_ref_old, coord_ref_new, alt_key_in, alt_key_out,
                 az_key_in, az_key_out, pol_key_in=None, pol_key_out=None):
        """
        ChangeCoordSys converts the pointing and, optionally, polarization
        angle, between two coordinate systems.
        """
        super(ChangeCoordSys, self).__init__()
        self.coord_ref_old = coord_ref_old
        self.coord_ref_new = coord_ref_new
        self.alt_key_in = alt_key_in
        self.alt_key_out = alt_key_out
        self.az_key_in = az_key_in
        self.az_key_out = az_key_out
        self.pol_key_in = pol_key_in
        self.pol_key_out = pol_key_out

    def Process(self, frame):
        if (not frame or frame.type is not G3FrameType.Map):
            return
        if (self.coord_ref_old == MapCoordReference.Equatorial and
                self.coord_ref_new == MapCoordReference.Galactic):
            alt_old = frame[self.alt_key_in]
            az_old = frame[self.az_key_in]
            alt_new, az_new = coord_eq_to_gal(alt_old/units.rad, az_old/units.rad)
            frame[self.alt_key_out] = G3VectorDouble(alt_new*units.rad)
            frame[self.az_key_out] = G3VectorDouble(az_new*units.rad)
            if self.pol_key_in is not None and self.pol_key_out is not None:
                pol_old = frame[self.pol_key_in]
                pol_new = pol_angle_eq_to_gal(alt_old/units.rad, az_old/units.rad,
                                              pol_old/units.rad)
                frame[self.pol_key_out] = pol_new*units.rad

        elif (self.coord_ref_old == MapCoordReference.Galactic and
                self.coord_ref_new == MapCoordReference.Equatorial):
            alt_old = frame[self.alt_key_in]
            az_old = frame[self.az_key_in]
            alt_new, az_new = coord_gal_to_eq(alt_old/rad, az_old/rad)
            frame[self.alt_key_out] = alt_new*rad
            frame[self.az_key_out] = az_new*rad
            if self.pol_key_in is not None and self.pol_key_out is not None:
                pol_old = frame[self.pol_key_in]
                pol_new = pol_angle_gal_to_eq(alt_old/units.rad, az_old/units.rad,
                                              pol_old/units.rad)
                frame[self.pol_key_out] = pol_new*units.rad
        else:
            raise NotImplementedError(
                    "Don't know how to convert from {0} to {1}".format(
                    coord_ref_old, coord_ref_new))
