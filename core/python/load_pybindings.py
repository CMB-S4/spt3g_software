import os, platform, sys

if platform.system().startswith('freebsd') or platform.system().startswith('FreeBSD'):
    # C++ modules are extremely fragile when loaded with RTLD_LOCAL,
    # which is what Python uses on FreeBSD by default, and maybe other
    # systems. Convince it to use RTLD_GLOBAL.

    # See thread by Abrahams et al:
    # http://mail.python.org/pipermail/python-dev/2002-May/024074.html
    sys.setdlopenflags(0x102)

try:
    from .. import dload

    have_dload = True
except ImportError:
    have_dload = False

if have_dload:
    lib_prefix = "libspt3g-"

    if platform.system().startswith('Darwin'):
        # OSX compatibility requires .dylib suffix
        lib_suffix = ".dylib"
    else:
        lib_suffix = ".so"
else:
    import importlib


def load_pybindings(name, path, deps=None):
    mod = sys.modules[name]

    if have_dload:
        p = os.path.split(path[0])
        m = dload.load_dynamic(name, p[1], p[0] + "/" + lib_prefix + p[1] + lib_suffix)
        # Don't override Python mod with C++
        sys.modules[name] = mod
    else:
        deps = deps or []
        if isinstance(deps, str):
            deps = [deps]
        if not deps and not name.endswith(".core"):
            deps = ["spt3g.core"]
        for dep in deps:
            # import dependencies
            importlib.import_module(dep, name)

        package, modname = name.rsplit(".", 1)
        m = importlib.import_module(f"{package}._lib{modname}")

    for (k,v) in m.__dict__.items():
        if not k.startswith("_"):
            mod.__dict__[k] = v
