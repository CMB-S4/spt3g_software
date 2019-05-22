from spt3g.core import G3Module, G3Pipeline, log_fatal
try:
    from spt3g.core import multiprocess
    multiproc_avail = True
except ImportError:
    multiproc_avail = False
    pass
import types

def pipesegment(func, autodoc=True):
    '''
    Use as a decorator for a pre-assembled set of pipeline modules. Makes a pseudo-module consisting of several inputs. If autodoc is True (the default), will attempt to introspect the segment to find out what it does. Set this to False if your module does anything complicated.

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
        rstintrodoc += 'Equivalent to:\n\n'
        doclines = []
        class PotemkinPipe(object):
            def Add(self, thing, *args, **kwargs):
                doc = 'pipe.Add(%s.%s' % (thing.__module__, thing.__name__)
                for arg in args:
                    doc += ', %s' % repr(arg)
                for arg in kwargs:
                    doc += ', %s=%s' % (arg, repr(kwargs[arg]))
                doc += ')'
                doclines.append(doc)
        fake = PotemkinPipe()
        try:
            func(fake)
            introdoc += '\n'.join(doclines)
            rstintrodoc += '.. code-block:: python\n\n    '
            rstintrodoc += '\n    '.join(doclines)
            rstintrodoc += '\n'
        except Exception as e:
            introdoc += 'Exception evaluating equivalence (%s)' % (str(e), )
            rstintrodoc += 'Exception evaluating equivalence (%s)' % (str(e), )
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

def PipelineAddCallable(self, callable, name=None, subprocess=False, *args, **kwargs):
    '''
    Add a processing module to the pipeline. It can be any subclass of
    spt3g.core.G3Module or any Python callable, either an instance or
    a class. Positional and keyword arguments are passed through to the
    argument's constructor (if a class) or as additional arguments to a
    function. If subprocess is set to True, this module will be
    run in a separate process.
    '''

    if not hasattr(self, 'nameprefix'):
        self.nameprefix = ''
    if name is None:
        if (hasattr(callable, '__name__')):
            callable_name = callable.__name__
        elif (hasattr(callable, '__class__')):
            callable_name = callable.__class__.__name__
        else:
            raise RuntimeError("Cannot establish name of pipeline module")
        name = '%s.%s' % (callable.__module__, callable_name)
    name = self.nameprefix + name

    # Deal with the segment case
    if hasattr(callable, '__pipesegment__'):
        # Prefix module names added by segments with the segment name
        oldnameprefix = self.nameprefix
        self.nameprefix = name + '/'

        rv = callable(self, *args, **kwargs)

        self.nameprefix = oldnameprefix
        return rv

    # Otherwise it's a module
    if subprocess:
        if not multiproc_avail:
            raise ImportError('Multiprocess not available')
        return self._Add_(build_pymodule(multiprocess.Subproc(build_pymodule(callable, *args, **kwargs), name=name)), name=name)
    else:
        return self._Add_(build_pymodule(callable, *args, **kwargs), name=name)

# Add this as G3Pipeline's Add method so it takes any Python callable
G3Pipeline.Add = PipelineAddCallable

