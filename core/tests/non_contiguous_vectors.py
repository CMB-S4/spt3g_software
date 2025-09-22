from spt3g import core
import numpy as np

data = np.random.randn(10, 10)

# check that vectors created from non-contiguous input arrays
# store the correct elements

vec = core.G3VectorDouble(data[:, 0])
np.testing.assert_array_equal(np.asarray(vec), data[:, 0])

vec = core.G3VectorDouble(data[0, :])
np.testing.assert_array_equal(np.asarray(vec), data[0, :])
