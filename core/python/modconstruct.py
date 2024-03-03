from spt3g.core import G3Module, G3Pipeline, G3PipelineInfo, G3Frame, G3FrameType, G3Time, G3ModuleConfig, log_fatal
try:
    from spt3g.core import multiprocess
    multiproc_avail = True
except ImportError:
    multiproc_avail = False
    pass
import types
import re

def pipesegment(func, autodoc=True):
    '''
    Use as a decorator for a pre-assembled set of pipeline modules. Makes a
    pseudo-module consisting of several inputs. If autodoc is True (the
    default), will attempt to introspect the segment to find out what it
    does. Set this to False if your module does anything complicated.

    For example:

    @core.pipesegment
    def standardfiltering(pipe, PolyOrder=4, MaskedHighPassEll=6000, Input='CalTimestreams', Output='FilteredTimestreams'):
        pipe.Add(analysis.PolyFilter, PolyOrder=PolyOrder, Input=Input,
            Output='__Temp' + Output)
        pipe.Add(analysis.MaskedHighPass, MaskedHighPassEll=MaskedHighPassEll, Input='__Temp' + Output, Output=Output)
        def cleanup(frame):
            del frame['__Temp' + Output]
        pipe.Add(cleanup)

    pipe.Add(standardfiltering, PolyOrder=3)
    '''

    func.__pipesegment__ = True

    if autodoc:
        introdoc = ''
        if hasattr(func, '__doc__') and func.__doc__ is not None:
            introdoc = func.__doc__ + '\n\n'
        from .docparser import format_doc
        rstintrodoc = format_doc(introdoc)
        introdoc += 'Equivalent to:\n'
        rstintrodoc += '\n*Equivalent to:*\n\n'
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
            func(fake)
        except Exception as e:
            doclines.append('Exception evaluating equivalence (%s)' % (str(e), ))
        introdoc += '\n'.join(doclines)
        rstintrodoc += '.. code-block:: python\n\n    '
        rstintrodoc += '\n    '.join(doclines)
        rstintrodoc += '\n'
        func.__doc__ = introdoc
        func.__rstdoc__ = rstintrodoc

    return func

def pipesegment_nodoc(func):
    return pipesegment(func, autodoc = False)

def indexmod(func):
    '''
    Mark argument as a processing module that can be found by automated documentation tools. For example:

    @core.indexmod
    def dostuff(frame):
        dosomestuff()
    '''

    func.__g3module__ = True

    return func

def build_pymodule(pycallable, *args, **kwargs):
    '''Convert a python callable and arguments into a core.G3Module by hook or by crook'''

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
        from spt3g import version
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

    # Otherwise it's a module
    elif subprocess:
        if not multiproc_avail:
            raise ImportError('Multiprocess not available')
        rv = self._Add_(build_pymodule(multiprocess.Subproc(build_pymodule(callable, **kwargs), name=name)), name=name)
    else:
        rv = self._Add_(build_pymodule(callable, **kwargs), name=name)

    if addpipelineinfo:
        self._Add_(self._pipelineinfo, name='_pipelineinfo')

    return rv

# Add this as G3Pipeline's Add method so it takes any Python callable
G3Pipeline.Add = PipelineAddCallable
G3Pipeline.__repr__ = lambda self: repr(self._pipelineinfo.pipelineinfo) if hasattr(self, '_pipelineinfo') else 'pipe = spt3g.core.G3Pipeline()'

