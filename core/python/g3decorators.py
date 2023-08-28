from spt3g.core import G3FrameType
from copy import copy
import inspect
import textwrap
import re

def get_function_signature(f, replacement_kwargs=None):
    '''
    Returns the function signature of f.  If replacement_kwargs is supplied
    replaces the default values with the ones in replacement_kwargs.
    '''
    sig = inspect.signature(f)

    if replacement_kwargs is not None:
        params = []
        for name, par in sig.parameters.items():
            if name in replacement_kwargs:
                par = par.replace(default=replacement_kwargs[name])
            params.append(par)
        sig = sig.replace(parameters=params)

    return f.__name__ + re.sub('<(.*) at (.*)>', '<\\1>', str(sig))


class cache_frame_data(object):
    '''
    This is a decorator for use with G3Modules written as functions.  It enables
    a function to use cached values from other types of frames in the processing
    of a frame.
    
    To make that confusing sentence clearer with an example, in a lot of cases
    we want to have a module that works on Scan frames, but have access to the
    BolometerProperties.  This decorator allows you to specify the information
    to cache.  This case looks like:

    .. code-block:: python

        @core.cache_frame_data(type=core.G3FrameType.Scan, bolo_props='BolometerProperties')
        def FlagSomeStuff(frame, flag_key='Flags', bolo_props=None):
            pass

    You specify the type of frame the function is operating on with the type
    argument.  Any additional keyword arguments specifies information to cache
    and send to the function.

    For the keyword args passed to the decorator having the format: Key = Value.
    Key specifies the name of the argument that we pass the infromation to in
    the decorated function.  If the decorated function is called with Key as an
    argument it will overwrite the value specified in the decorator.
       
    Value specifies the default path to look for the cached data.  It will
    search all of the frames that do not have the frame type 'type' for a key
    with that value.  This can be overridden when calling the decorated
    function.
    '''

    def __init__(self_outer, type, **kwargs):
        self_outer.keyargs = kwargs
        self_outer.type = type
    def __call__(self_outer, f):
        func_sig = get_function_signature(f, replacement_kwargs = self_outer.keyargs)
        doc_prepend = '\nFunction Signature:  %s\n\n' % func_sig
        func_doc = textwrap.dedent(str(f.__doc__ or ''))
        doc_append = '''\n\n
This function automatic caches values from frames that contain the keys to be
cached.  It will only operate on %s frames.  If the cached keys are not found it
will pass None to the function.  You can change the default key by supplying the
following arguments:

Cached Values
-------------
''' % (str(self_outer.type or 'any type of') )
        for k,v in self_outer.keyargs.items():
            doc_append += '    %s = "%s"\n' % ( k, v )

        class WrappedFunc:
            def __init__(self, *args, **kwargs): 
                self.args = args
                self.kwargs = kwargs

                self.argument_map = copy(self_outer.keyargs)
                pop_ks = []
                for k in self.kwargs.keys():
                    if k in self.argument_map:
                        self.argument_map[k] = self.kwargs[k]
                        pop_ks.append(k)
                for k in pop_ks:
                    self.kwargs.pop(k)

            def __call__(self, frame):
                for vname, stored_key in self.argument_map.items():
                    if stored_key and stored_key in frame:
                        self.kwargs[vname] = frame[stored_key]
                if self_outer.type is None or frame.type == self_outer.type:
                    return f(frame, *(self.args), **(self.kwargs))
        WrappedFunc.__wrapped__ = f
        WrappedFunc.__name__ = f.__name__
        WrappedFunc.__doc__ = doc_prepend + func_doc + doc_append
        WrappedFunc.__g3module__ = True
        return WrappedFunc


class scan_func_cache_data(cache_frame_data):
    '''
    This is a simple wrapper around cache_frame_data where the type argument has been set to
    core.G3FrameType.Scan.
    '''
    def __init__(self, **kwargs):
        super(scan_func_cache_data, self).__init__(G3FrameType.Scan, **kwargs)
