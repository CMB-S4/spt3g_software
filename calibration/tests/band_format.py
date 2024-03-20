from spt3g import core
from spt3g.calibration import BolometerProperties, set_band_format, band_to_value

bp = BolometerProperties()
bp.band = 95.47 * core.G3Units.GHz

assert bp.band_string == "95GHz"
assert band_to_value(bp.band_string) == 95 * core.G3Units.GHz

set_band_format(0, "MHz")

assert bp.band_string == "95470MHz"
assert band_to_value(bp.band_string) == bp.band

set_band_format(2, "GHz")

assert bp.band_string == "95.47GHz"
assert band_to_value(bp.band_string) == bp.band
