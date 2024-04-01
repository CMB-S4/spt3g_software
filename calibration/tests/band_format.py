from spt3g import core
from spt3g.calibration import BolometerProperties, set_band_format, band_to_value, get_band_units

bp = BolometerProperties()
bp.band = 95.47 * core.G3Units.GHz

assert bp.band_string == "95GHz"
assert bp.band_vstring == "95"
assert band_to_value(bp.band_string) == 95 * core.G3Units.GHz
assert bp.band_vstring + get_band_units() == bp.band_string

set_band_format(0, "MHz")

assert bp.band_string == "95470MHz"
assert bp.band_vstring == "95470"
assert band_to_value(bp.band_string) == bp.band
assert bp.band_vstring + get_band_units() == bp.band_string

set_band_format(2, "GHz")

assert bp.band_string == "95.47GHz"
assert bp.band_vstring == "95.47"
assert band_to_value(bp.band_string) == bp.band
assert bp.band_vstring + get_band_units() == bp.band_string
