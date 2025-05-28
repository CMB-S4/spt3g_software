from . import G3Module, G3FrameObject, usefulfunc
import textwrap


__all__ = ["G3Documenter", "module_apidoc"]


def rst_header(title, line, overline=False):
    """
    Add a valid RST underline (and optional matching overline) to create a
    section header.
    """
    out = title + "\n" + line * len(title) + "\n\n"
    if overline:
        out = line * len(title) + "\n" + out
    return out


@usefulfunc
class G3Documenter:
    """
    Class for inspecting sub-modules of the SPT-3G software package and
    generating valid RST for use with sphinx-autodoc.
    """

    # API categories for which to generate documentation, with a
    # corresponding section title, description and sphinx-autodoc directive.
    _categories = {
        "object": (
            "Frame Objects",
            "Serializable objects that may be stored as entries in "
            ":py:class:`~spt3g.core.G3Frame` objects.",
            "autoclass",
        ),
        "cmodule": (
            "Class-like Pipeline Modules",
            """
            Classes that can be added to :py:class:`~spt3g.core.G3Pipeline`
            instances, either written as callable objects in Python (see
            :ref:`class-modules`) or as :py:class:`~spt3g.core.G3Module`-derived
            classes in C++ (see :ref:`cxx-modules`).  Aside from adding these into
            pipelines using the standard ``Add()`` method, such objects may also be
            instantiated outside of pipelines to make use of additional features.

            One typical use case is to extract some kind of data from a series of
            frames::

                class ExtractData:
                    def __init__(self, key="RawTimestreams"):
                        self.key = key
                        self.data = []
                    def __call__(self, frame):
                        if self.key in frame:
                            self.data.append(frame[self.key])

                # instantiate caching class
                extractor = ExtractData()

                # use instance in a pipeline
                pipe = core.G3Pipeline()
                pipe.Add(core.G3Reader, filename="file.g3")
                pipe.Add(extractor)
                pipe.Run()

                # use extracted data in later processing
                print(extractor.data)

            Another use case may be to call additional methods of a class.  For
            example, the :py:class:`~spt3g.core.G3Writer` class can be used as a
            context manager, and also report the position of the file pointer after
            each frame is written::

                with core.G3Writer("file.g3") as writer:
                    for frame in frames:
                        writer(frame)
                        print("Bytes written:", writer.tell())

            Class-like modules permit a variety of additional features beyond
            their standard usage in pipelines.
            """,
            "autoclass",
        ),
        "fmodule": (
            "Function-like Pipeline Modules",
            "Python functions that can be added to :py:class:`~spt3g.core.G3Pipeline` "
            "instances. Such functions may also be called directly with a "
            ":py:class:`~spt3g.core.G3Frame` object as the first argument, and do not "
            "necessarily need to be used in a pipeline. See :ref:`function-modules` "
            "for more detail.",
            "autofunction",
        ),
        "segment": (
            "Pipeline Segments",
            "Python functions that combine a sequence of pipeline modules, and can be "
            "added to :py:class:`~spt3g.core.G3Pipeline` instances. "
            "See :ref:`pipesegments` for more detail.",
            "autofunction",
        ),
        "class": (
            "Useful Classes",
            "Various Python and C++ classes that are part of the public API.",
            "autoclass",
        ),
        "function": (
            "Useful Functions",
            "Various Python and C++ functions that are part of the public API.",
            "autofunction",
        ),
        "decorator": (
            "Decorators",
            "Decorator functions for indicating any of the above types and/or "
            "modifying function behavior.",
            "autodecorator",
        )
    }

    def __init__(self, root):
        self.root = root
        self.cache = {}

    @property
    def categories(self):
        """
        Dictionary of API categories for which to generate documentation,
        containing a section title, description, and sphinx-autodoc directive as
        a tuple for each category.
        """
        return self._categories

    def find_objects(self, mod=None):
        """
        Recursively find all objects in the input root that fit into one of the
        API categories listed in ``categories``.

        Arguments
        ---------
        mod : object
            Object to inspect

        Returns
        -------
        cache : dict
            Recursively populated dictionary of fully qualified names and their
            associated API category.
        """
        if mod is None:
            mod = __import__(self.root)
            name = self.root.split(".")[-1]
            if '.' in self.root and hasattr(mod, name):
                mod = mod.__dict__[name]

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
            if not itemname.startswith(self.root):
                continue
            # skip already-classified objects
            if itemname in self.cache:
                continue

            # classify
            isclass = isinstance(obj, type)
            ismod = isclass and issubclass(obj, G3Module)
            isobj = isclass and issubclass(obj, G3FrameObject)
            iscxx = 'Boost.Python' in str(type(obj))

            # add to cache
            if hasattr(obj, "__g3decorator__"):
                self.cache[itemname] = "decorator"
            elif ismod or hasattr(obj, "__g3module__"):
                self.cache[itemname] = "cmodule" if isclass else "fmodule"
            elif isobj or hasattr(obj, '__g3frameobject__'):
                self.cache[itemname] = "object"
            elif hasattr(obj, "__pipesegment__"):
                self.cache[itemname] = "segment"
            elif iscxx or hasattr(obj, '__g3usefulfunc__'):
                self.cache[itemname] = "class" if isclass else "function"
            else:
                self.cache[itemname] = "skip"
                if not isclass and hasattr(obj, '__dict__'):
                    # recurse
                    self.find_objects(obj)

        # drop all skipped objects
        if self.root == mod.__name__:
            for k in list(self.cache):
                if self.cache[k] == "skip":
                    self.cache.pop(k)

        return self.cache

    def generate(self):
        """
        Construct valid RST for each of the API types in separate sections, and
        list all relevant objects with appropriate sphinx-autodoc directives.

        Returns
        -------
        str :
            String containing valid RST to append to a document
        """
        # find all relevant objects
        cache = self.find_objects()
        if not len(cache):
            return

        txt = rst_header("API Documentation", "=")

        # build API sections
        body = ""
        for kind, (title, description, directive) in self.categories.items():
            objs = sorted([k for k, v in cache.items() if v == kind])
            if not len(objs):
                continue

            body += rst_header(title, "-")
            body += textwrap.dedent(description) + "\n\n"
            body += "\n\n".join([f".. {directive}:: {k}" for k in objs]) + "\n\n"

        if len(body):
            txt += ".. contents:: Contents\n   :depth: 1\n   :local:\n\n"
            txt += body

        return txt


@usefulfunc
def module_apidoc(module_path):
    """
    Create API documentation for a submodule of the spt3g library.

    The output separates each of the API types into separate sections, and lists
    all relevant objects with appropriate sphinx-autodoc directives.

    Arguments
    ---------
    module_path : str
        Python module path to inspect, e.g. "spt3g.core"

    Returns
    -------
    str :
        String containing valid RST to append to a document
    """

    return G3Documenter(module_path).generate()
