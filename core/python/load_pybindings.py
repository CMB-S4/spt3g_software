import platform,sys

if platform.system().startswith('freebsd') or platform.system().startswith('FreeBSD'):
	# C++ modules are extremely fragile when loaded with RTLD_LOCAL,
	# which is what Python uses on FreeBSD by default, and maybe other
	# systems. Convince it to use RTLD_GLOBAL.

	# See thread by Abrahams et al:
	# http://mail.python.org/pipermail/python-dev/2002-May/024074.html
	sys.setdlopenflags(0x102)

lib_prefix = "libspt3g-"

if platform.system().startswith('Darwin'):
    # OSX compatibility requires .dylib suffix
    lib_suffix = ".dylib"
else:
    lib_suffix = ".so"

def load_pybindings(name, path):
	import os
	mod = sys.modules[name]
	p = os.path.split(path[0])
	from spt3g import dload
	m = dload.load_dynamic(name, p[1], p[0] + "/" + lib_prefix + p[1] + lib_suffix)
	sys.modules[name] = mod # Don't override Python mod with C++

	for (k,v) in m.__dict__.items():
		if not k.startswith("_"):
			mod.__dict__[k] = v
