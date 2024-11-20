from spt3g import core
from spt3g.calibration import BolometerProperties, set_band_format, band_to_value, get_band_units

bp = BolometerProperties()
bp.band = 94.67 * core.G3Units.GHz

assert bp.band_string == "95GHz"
assert bp.band_vstring == "95"
assert band_to_value(bp.band_string) == 95 * core.G3Units.GHz
assert bp.band_vstring + get_band_units() == bp.band_string

set_band_format(0, "MHz")

assert bp.band_string == "94670MHz"
assert bp.band_vstring == "94670"
assert band_to_value(bp.band_string) == bp.band
assert bp.band_vstring + get_band_units() == bp.band_string

set_band_format(2, "GHz")

assert bp.band_string == "94.67GHz"
assert bp.band_vstring == "94.67"
assert band_to_value(bp.band_string) == bp.band
assert bp.band_vstring + get_band_units() == bp.band_string

set_band_format(-1, "GHz")
assert bp.band_string == "90GHz"
assert bp.band_vstring == "90"
assert band_to_value(bp.band_string) == 90 * core.G3Units.GHz
assert bp.band_vstring + get_band_units() == bp.band_string
