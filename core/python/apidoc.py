from . import G3Module, G3FrameObject


__all__ = ["module_apidoc"]


# API categories for which to generate documentation, with a
# corresponding section title and sphinx-autodoc directive.
api_categories = {
    "object": ("Frame Objects", "autoclass"),
    "cmodule": ("Class-like Pipeline Modules", "autoclass"),
    "fmodule": ("Function-like Pipeline Modules", "autofunction"),
    "segment": ("Pipeline Segments", "autofunction"),
    "class": ("Useful Classes", "autoclass"),
    "function": ("Useful Functions", "autofunction"),
}


def find_objects(mod, cache, root):
    """
    Resursively find all objects in the input module that fit into one of the
    API categories listed in ``api_categories``.

    Arguments
    ---------
    mod : object
        Object to inspect
    cache : dict
        Dictionary of fully qualified names and their associated API category
    root : str
        Root name under which all objects in this cache should belong.  Any
        items outside of this namespace are not included.
    """
    for x, obj in mod.__dict__.items():
        if x.startswith("_") or x.endswith("_"):
            continue

        # construct fully-qualified object name
        if hasattr(obj, "__wrapped__"):
            obj = obj.__wrapped__
        try:
            if hasattr(obj, '__module__'):
                itemname = '%s.%s' % (obj.__module__, obj.__name__)
            else:
                itemname = obj.__name__
        except:
            continue

        # skip objects outside of namespace
        if not itemname.startswith(root):
            continue
        # skip already-classified objects
        if itemname in cache:
            continue

        # classify
        isclass = isinstance(obj, type)
        ismod = isclass and issubclass(obj, G3Module)
        isobj = isclass and issubclass(obj, G3FrameObject)
        iscxx = 'Boost.Python' in str(type(obj))

        # add to cache
        if ismod or hasattr(obj, "__g3module__"):
            cache[itemname] = "cmodule" if isclass else "fmodule"
        elif isobj or hasattr(obj, '__g3frameobject__'):
            cache[itemname] = "object"
        elif hasattr(obj, "__pipesegment__"):
            cache[itemname] = "segment"
        elif iscxx or hasattr(obj, '__g3usefulfunc__'):
            cache[itemname] = "class" if isclass else "function"
        else:
            cache[itemname] = "skip"
            if not isclass and hasattr(obj, '__dict__'):
                # recurse
                find_objects(obj, cache, root)


def rst_header(title, line, overline=False):
    """
    Add a valid RST underline (and optional matching overline) to create a
    section header.
    """
    out = title + "\n" + line * len(title) + "\n\n"
    if overline:
        out = line * len(title) + "\n" + out
    return out


def module_apidoc(module_path):
    """
    Create API documentation for a submodule of the spt3g library.

    The output separates each of the API types into separate sections, and lists
    all relevant objects with appropriate sphinx-autodoc directives.
    """

    mod = __import__(module_path)
    name = module_path.split(".")[-1]
    if '.' in module_path and hasattr(mod, name):
        mod = mod.__dict__[name]

    # find all relevant objects
    cache = {}
    find_objects(mod, cache, module_path)

    # drop all skipped objects
    for k in list(cache):
        if cache[k] == "skip":
            cache.pop(k)
    if not len(cache):
        return

    txt = rst_header("API Documentation", "=")

    # build API sections
    for kind, (title, directive) in api_categories.items():
        objs = sorted([k for k, v in cache.items() if v == kind])
        if not len(objs):
            continue

        txt += rst_header(title, "-")
        txt += "\n\n".join([f".. {directive}:: {k}" for k in objs]) + "\n\n"

    return txt
