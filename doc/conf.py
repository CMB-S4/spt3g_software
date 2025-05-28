from sphinx.ext import autodoc
import pprint
import re
import textwrap
import spt3g

# sphinx configuration options

project = "SPT-3G Software"
version = spt3g.__version__
copyright = "2017-, SPT-3G Collaboration"
author = "SPT-3G Collaboration"
master_doc = "index"
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
}
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_logo = "spt_outlines.png"
html_favicon = "spt_outlines_white.png"
exclude_patterns = ["intro_*"]
extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
]

# docstring parser configuration
napoleon_custom_sections = [
    ("Emits", "returns_style"),
    ("Equivalent to", "Example")
]
maximum_signature_line_length = 85
python_use_unqualified_type_names = True

# autodoc configuration
autodoc_class_signature = "mixed"
autoclass_content = "both"
autodoc_member_order = "groupwise"
autodoc_typehints = "description"
autodoc_typehints_format = "short"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
    "class-doc-from": "both",
    "member-order": "groupwise",
    "exclude-members": "__new__",
}

# custom documenters for handling C++ signatures and overloaded functions

selfarg = re.compile(r"(\(\s*\([\w.]+\)arg1\s\[,\s)", re.VERBOSE)
selfarg1 = re.compile(r"(\(\s*\([\w.]+\)arg1(?:,\s*)?)", re.VERBOSE)

optarg1 = re.compile(r"\(\[([\(\)\w._]+=[\w.\(\)=,\s_\-\'\"]+)\]\)", re.VERBOSE)
optarg2 = re.compile(r"\[,(\s*[\(\)\w._]+=[\w.\(\)=,\s_\'\"\-\[\]]+)\]", re.VERBOSE)

typearg = re.compile(r"\(([\w.]+)\)([\w._]+)([=,\)]?)", re.VERBOSE)

arg2 = re.compile(r"\(arg2:(\s*[\w.]+)\)", re.VERBOSE)


class BindingDocsMixin:
    """
    Mixin class to handle signature cleanup in C++ bindings,
    including annotations and overloaded signatures
    """
    def get_doc(self):
        if self._new_docstrings is not None:
            return self._new_docstrings

        docstrings = super().get_doc()
        if not docstrings:
            return docstrings

        if 'Boost.Python' not in str(type(self.object)):
            if "->" not in str(docstrings):
                return docstrings

        for i, doclines in enumerate(docstrings):
            # drop empty lines
            while not doclines[0].strip():
                doclines.pop(0)

            # collect all signature lines
            sigs = []
            docs = []
            docs1 = []
            for j, line in enumerate(doclines):
                # CXX bindings typically have return annotations
                if not "->" in line:
                    docs1.append(line)
                    continue

                # remove unnecessary return annotations
                if line.strip().endswith(":"):
                    line = line.rstrip().rstrip(":").rstrip()
                if line.startswith("__init__"):
                    line = line.partition(" -> ")[0]
                if line.endswith("-> None"):
                    line = line.replace(" -> None", "")
                if isinstance(
                    self, (autodoc.ClassLevelDocumenter, autodoc.ClassDocumenter)
                ):
                    # remove self from method argument list
                    m = selfarg.search(line)
                    if m:
                        line = line.replace(m.group(0), "([")
                    m = selfarg1.search(line)
                    if m:
                        line = line.replace(m.group(0), "(")

                # mangle signatures to produce python-like type annotations
                m = optarg2.search(line)
                while m:
                    line = line.replace(m.group(0), "," + m.group(1))
                    m = optarg2.search(line)
                m = optarg1.search(line)
                if m:
                    line = line.replace(m.group(0), "(" + m.group(1) + ")")
                m = typearg.search(line)
                while m:
                    line = line.replace(
                        m.group(0), m.group(2) + ": " + m.group(1) + m.group(3)
                    )
                    m = typearg.search(line)
                m = arg2.search(line)
                if m:
                    line = line.replace(m.group(0), "(arg:" + m.group(1) + ")")
                line = line.replace("()", "(\u200b)")
                line = re.sub('<(.*) at (.*)>', '<\\1>', str(line))

                # collate overloaded signatures and corresponding docstrings
                if line.strip():
                    if "".join(docs1).strip():
                        if sigs:
                            n = len(docs1[0]) - len(docs1[0].lstrip())
                            docs1[0] = docs1[0][:n] + "Signature {}: {}".format(
                                len(sigs), docs1[0].lstrip()
                            )
                        docs.extend(docs1 + [""])
                    sigs.append(line)
                    docs1.clear()
                doclines[j] = line

            if len(sigs):
                # finish collation
                doclines.clear()
                doclines.extend(sigs)
                if "".join(docs1).strip():
                    if len(sigs) > 1 and (len(docs) or i > 0):
                        n = len(docs1[0]) - len(docs1[0].lstrip())
                        docs1[0] = docs1[0][:n] + "Signature {}: {}".format(
                            len(sigs), docs1[0].lstrip()
                        )
                    docs.extend(docs1 + [""])
                doclines.extend(docs)

        return docstrings


# Create mixin subclasses for all relevant documenters
class FunctionDocumenter(BindingDocsMixin, autodoc.FunctionDocumenter):
    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        # ensure that pipesegments are documented like functions
        if hasattr(member, "__pipesegment__"):
            return True
        return super().can_document_member(member, membername, isattr, parent)


class MethodDocumenter(BindingDocsMixin, autodoc.MethodDocumenter):
    pass


class PropertyDocumenter(autodoc.PropertyDocumenter):
    def get_doc(self):
        if self._new_docstrings is not None:
            return self._new_docstrings

        docstrings = super().get_doc()
        if not docstrings or "->" not in str(docstrings):
            func = self._get_property_getter()
            if func is None or getattr(func, "__doc__", None) is None:
                return docstrings
            gdoc = func.__doc__
            if any([gdoc.strip() in "\n".join(d).strip() for d in docstrings]):
                return docstrings
            docstrings.append(gdoc.split("\n"))

        for doclines in docstrings:
            for j, line in enumerate(doclines):
                if "->" in line:
                    self._cxx_type_ann = line.partition(" -> ")[-1].strip()
                    doclines[j] = ""
                    return docstrings

        return docstrings

    def add_directive_header(self, sig):
        # update property type annotation
        super().add_directive_header(sig)
        sourcename = self.get_sourcename()
        if hasattr(self, "_cxx_type_ann"):
            self.add_line("   :type: " + self._cxx_type_ann, sourcename)


# better formatting for enum name/value dicts
class AttributeDocumenter(autodoc.AttributeDocumenter):
    def import_object(self, raiseerror=False):
        self._is_enum_dict = False
        ret = super().import_object(raiseerror)
        if isinstance(self.object, dict):
            if self.objpath[-1] in ["names", "values"]:
                self._is_enum_dict = True
        return ret

    def should_suppress_value_header(self):
        if self._is_enum_dict:
            return True
        return super().should_suppress_value_header()

    def add_directive_header(self, sig):
        super().add_directive_header(sig)
        if self._is_enum_dict:
            sourcename = self.get_sourcename()
            self.add_line("   :type: dict", sourcename)

    def get_doc(self):
        if not self._is_enum_dict:
            return super().get_doc()

        if self.objpath[-1] == "names":
            txt = "Mapping of enum names to objects."
        else:
            txt = "Mapping of enum values to objects."

        v = pprint.pformat(self.object)
        v = textwrap.indent(v, "    ")
        return [[txt, "", ".. code-block::", "", ""] + v.split("\n")]


class ClassDocumenter(BindingDocsMixin, autodoc.ClassDocumenter):
    def add_line(self, line, source, *lineno):
        # remove empty bases directive
        if line.strip() == "Bases:":
            return
        super().add_line(line, source, *lineno)


def handle_bases(app, name, obj, options, bases: list):
    """
    Exclude generic base classes from inheritance lists
    """
    for b in list(bases):
        if 'Boost.Python' in str(b):
            bases.remove(b)
    if object in bases:
        bases.remove(object)


def setup(app):
    app.connect('autodoc-process-bases', handle_bases)
    app.add_autodocumenter(FunctionDocumenter, override=True)
    app.add_autodocumenter(MethodDocumenter, override=True)
    app.add_autodocumenter(PropertyDocumenter, override=True)
    app.add_autodocumenter(AttributeDocumenter, override=True)
    app.add_autodocumenter(ClassDocumenter, override=True)
