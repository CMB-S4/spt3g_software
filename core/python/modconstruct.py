from . import G3Module, G3Pipeline, G3PipelineInfo, G3Frame, G3FrameType, G3Time, G3ModuleConfig
import types
import re
import textwrap


def g3decorator(func):
    """
    Mark argument as a decorator that can be found by automated
    documentation tools.
    """
    func.__g3decorator__ = True
    return func


@g3decorator
def usefulfunc(func):
    '''
    Mark argument as a useful function that can be found by automated
    documentation tools.

    Example
    -------
    ::

        @core.usefulfunc
        def do_some_science(data):
            science(data)
    '''

    func.__g3usefulfunc__ = True

    return func

# document thyself
usefulfunc.__g3usefulfunc__ = True


class PipeSegmentDocstring(type):
    @property
    def __doc__(cls):
        return '''
    Use as a decorator for a pre-assembled set of pipeline modules. Makes a
    pseudo-module consisting of several inputs. Use this to introspect the
    segment to find out what it does, or use :func:`pipesegment_nodoc` if your
    module does anything complicated.

    Example
    -------
    ::

        @core.pipesegment
        def standardfiltering(
            pipe,
            PolyOrder=4,
            MaskedHighPassEll=6000,
            Input='CalTimestreams',
            Output='FilteredTimestreams',
        ):
            pipe.Add(analysis.PolyFilter, PolyOrder=PolyOrder, Input=Input,
                Output='__Temp' + Output)
            pipe.Add(analysis.MaskedHighPass, MaskedHighPassEll=MaskedHighPassEll,
                Input='__Temp' + Output, Output=Output)
            def cleanup(frame):
                del frame['__Temp' + Output]
            pipe.Add(cleanup)

        pipe.Add(standardfiltering, PolyOrder=3)
    '''


@g3decorator
class pipesegment(metaclass=PipeSegmentDocstring):
    def __init__(self, func):
        self.func = self.__wrapped__ = func
        self.__pipesegment__ = func.__pipesegment__ = True

    def __call__(self, pipe, *args, **kwargs):
        return self.func(pipe, *args, **kwargs)

    @property
    def __doc__(self):
        """
        Create a dummy pipeline for introspection.  Generates a docstring when
        the __doc__ attribute is accessed.
        """
        if hasattr(self, "_autodoc"):
            return self._autodoc

        introdoc = textwrap.dedent(getattr(self.func, "__doc__", None) or "") or ""
        if introdoc:
            introdoc += "\n\n"
        introdoc += '\nEquivalent to\n-------------\n\n::\n\n'
        doclines = []
        class PotemkinPipe(object):
            def Add(self, thing, *args, **kwargs):
                if hasattr(thing, '__wrapped__'):
                    modname = thing.__wrapped__.__module__
                else:
                    modname = thing.__module__
                doc = 'pipe.Add(%s.%s' % (modname, thing.__name__)
                for arg in args:
                    doc += ', %s' % repr(arg)
                for arg in kwargs:
                    # remove object hashes
                    s = re.sub('<(.*) at (.*)>', '<\\1>', repr(kwargs[arg]))
                    doc += ', %s=%s' % (arg, s)
                doc += ')'
                doclines.append(doc)
        fake = PotemkinPipe()
        try:
            self.func(fake)
        except Exception as e:
            doclines.append('Exception evaluating equivalence (%s)' % (str(e), ))
        introdoc += '    ' + '\n    '.join(doclines) + '\n'

        self._autodoc = introdoc
        return self._autodoc

    @property
    def __name__(self):
        return self.func.__name__


@g3decorator
def pipesegment_nodoc(func):
    """
    Use as a decorator for a pre-assembled set of pipeline modules. Makes a
    pseudo-module consisting of several inputs.  Use this variant instead of
    :class:`pipesegment` to avoid introspection if your pipeline does anything
    complicated.

    Example
    -------
    ::

        @core.pipesegment_nodoc
        def standardfiltering(
            pipe,
            PolyOrder=4,
            MaskedHighPassEll=6000,
            Input='CalTimestreams',
            Output='FilteredTimestreams',
        ):
            pipe.Add(analysis.PolyFilter, PolyOrder=PolyOrder, Input=Input,
                Output='__Temp' + Output)
            pipe.Add(analysis.MaskedHighPass, MaskedHighPassEll=MaskedHighPassEll,
                Input='__Temp' + Output, Output=Output)
            def cleanup(frame):
                del frame['__Temp' + Output]
            pipe.Add(cleanup)

        pipe.Add(standardfiltering, PolyOrder=3)
    """
    func.__pipesegment__ = True
    return func


@g3decorator
def indexmod(func):
    '''
    Mark argument as a processing module that can be found by automated
    documentation tools.

    Example
    -------
    ::

        @core.indexmod
        def dostuff(frame):
            dosomestuff()
    '''

    func.__g3module__ = True

    return func

def build_pymodule(pycallable, *args, **kwargs):
    '''Convert a python callable and arguments into a core.G3Module by hook or by crook'''
    from .g3logging import log_fatal

    if isinstance(pycallable, G3Module):
        return pycallable

    if not callable(pycallable):
        log_fatal('Argument not a python callable', unit = 'G3Pipeline')

    if type(pycallable) == types.FunctionType:
        class PyFuncModule(G3Module):
            def Process(self, fr):
                return pycallable(fr, *args, **kwargs)
        return PyFuncModule()

    # If this is not a function, and it is callable, it is a class. If it is
    # a non-instantiated class, instantiate it.

    try:
        # issubclass() throws TypeError if argument is instantiated
        issubclass(pycallable, G3Module)
        isclass = True
    except TypeError:
        isclass = False

    # Instantiate if necessary. Not in try block so we don't miss TypeError
    if isclass:
        pycallable = pycallable(*args, **kwargs)
    else:
        # No way to handle arguments in this call; assert there are none
        if len(args) != 0 or len(kwargs) != 0:
            log_fatal('Cannot pass through arguments when passed instantiated class', unit = 'G3Pipeline')

    # See if it was a Python G3Module subclass
    if isinstance(pycallable, G3Module):
        return pycallable

    # This is a python callable that is not a module, so wrap it
    class PyCallObjModule(G3Module):
        def Process(self, fr):
            return pycallable(fr)

    return PyCallObjModule()

class _add_pipeline_info(G3Module):
    def __init__(self):
        from .. import version
        import socket, getpass

        G3Module.__init__(self)
        self.infoemitted = False
        self.buffer = []

        self.pipelineinfo = G3PipelineInfo()
        self.pipelineinfo.vcs_url = version.upstream_url
        self.pipelineinfo.vcs_branch = version.upstream_branch
        self.pipelineinfo.vcs_revision = version.revision
        self.pipelineinfo.vcs_localdiffs = version.localdiffs
        self.pipelineinfo.vcs_versionname = version.versionname
        self.pipelineinfo.vcs_fullversion = version.fullversion
        self.pipelineinfo.vcs_githash = version.gitrevision

        self.pipelineinfo.hostname = socket.gethostname()
        self.pipelineinfo.user = getpass.getuser()
    def Process(self, fr):
        if self.infoemitted:
            if fr.type == G3FrameType.PipelineInfo and \
              self.originalpi is not None and \
              self.originalpi == list(fr.keys()):
                # Deduplicate PipelineInfo frames identical to one that we
                # added to earlier, which avoids false semi-duplicates down
                # the line when processing multiple files.
                return False
            return fr

        # Allow limited reordering of metadata
        if fr.type in [G3FrameType.Observation, G3FrameType.Wiring,
          G3FrameType.Calibration]:
            self.buffer.append(fr)
            return []

        self.infoemitted = True
        if fr.type == G3FrameType.PipelineInfo:
            self.originalpi = list(fr.keys())
            fr[str(G3Time.Now())] = self.pipelineinfo
            self.buffer.append(fr)
        else:
            self.originalpi = None
            f = G3Frame(G3FrameType.PipelineInfo)
            f[str(G3Time.Now())] = self.pipelineinfo
            self.buffer += [f, fr]
        rv = self.buffer
        del self.buffer
        return rv

def PipelineAddCallable(self, callable, name=None, subprocess=False, **kwargs):
    '''
    Add a processing module to the pipeline. It can be any subclass of
    spt3g.core.G3Module or any Python callable, either an instance or
    a class. Positional and keyword arguments are passed through to the
    argument's constructor (if a class) or as additional arguments to a
    function. If subprocess is set to True, this module will be
    run in a separate process.
    '''

    addpipelineinfo = False
    if not hasattr(self, '_pipelineinfo'):
        self._pipelineinfo = _add_pipeline_info()
        addpipelineinfo = True

    if not hasattr(self, 'nameprefix'):
        self.nameprefix = ''
    if (hasattr(callable, '__name__')):
        callable_name = callable.__name__
    elif (hasattr(callable, '__class__')):
        callable_name = callable.__class__.__name__
    else:
        raise RuntimeError("Cannot establish name of pipeline module")
    if name is None:
        name = '%s.%s' % (callable.__module__, callable_name)
    name = self.nameprefix + name

    # Record module configuration for root objects
    if self.nameprefix == '': 
        modconfig = G3ModuleConfig()
        modconfig.instancename = name
        modconfig.modname = '%s.%s' % (callable.__module__, callable_name)
        for k,v in kwargs.items():
            tostore = v
            try:
                if v.npix_allocated > 0:
                    # Don't store full sky maps as configuration options. It
                    # just wastes a ton of disk space with simulations.
                    tostore = v.clone(False)
            except:
                # If that threw an exception, it either isn't a map or dropping
                # data didn't work, so just don't bother.
                pass
            modconfig[k] = tostore
        self._pipelineinfo.pipelineinfo.modules.append(modconfig)

    # Deal with the segment case
    if hasattr(callable, '__pipesegment__'):
        # Prefix module names added by segments with the segment name
        oldnameprefix = self.nameprefix
        self.nameprefix = name + '/'

        rv = callable(self, **kwargs)

        self.nameprefix = oldnameprefix

    else:
        # Otherwise it's a module
        pymod = build_pymodule(callable, **kwargs)
        if subprocess:
            from .multiprocess import Subproc
            pymod = build_pymodule(Subproc(pymod, name=name))

        rv = self._Add_(pymod, name=name)

    if addpipelineinfo:
        self._Add_(self._pipelineinfo, name='_pipelineinfo')

    return rv

# Add this as G3Pipeline's Add method so it takes any Python callable
G3Pipeline.Add = PipelineAddCallable
G3Pipeline.__repr__ = lambda self: repr(self._pipelineinfo.pipelineinfo) if hasattr(self, '_pipelineinfo') else 'pipe = {}.G3Pipeline()'.format(G3Pipeline.__module__)

