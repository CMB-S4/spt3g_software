import sys, inspect, re, textwrap
from spt3g.core import G3Module, G3FrameObject

def format_doc(x, simple=False):
    """
    Apply some formatting cleanup to a docstring to make it play nice with the sphinx layout.
    """
    if isinstance(x, str):
        doc = x
    else:
        doc = inspect.getdoc(x)
    if not doc:
        return doc
    try:
        if doc == inspect.getdoc(type.__init__):
            return None
    except:
        pass
    doc = textwrap.dedent(doc)
    lines = doc.split('\n')
    doclines = []
    head = False
    for line in lines:
        if line.strip() and all([x == '-' for x in line.strip()]):
            if doclines:
                head = True
                hdr = doclines[-1].strip()
                if hdr[-1] != ':':
                    hdr += ':'
                doclines[-1] = '*{}*'.format(hdr)
        else:
            if simple and '->' in line and not line.startswith('    '):
                line = '``{}``'.format(line)
            doclines.append('\t' * head + line)
    return '\n'.join(doclines)

def get_doc_for_module(module_path, include_link_list = True):
    try:
        mod = __import__(module_path)
        if '.' in module_path and hasattr(mod, module_path.split('.')[-1]):
            mod = mod.__dict__[module_path.split('.')[-1]]
    
    except ImportError as e:
        print('Could not import module %s:' % module_path)
        print(e)
        sys.exit(1)
    

    def format_object(modname, x):
        name = '%s.%s' % (modname, x)
        return '\n.. _%s:\n\n%s\n%s\n' % (name, name, '-' * len(name))
    def format_name(modname, x):
        return '\n.. _%s.%s:\n\n**%s.%s**\n' % (modname, x, modname, x)
    def format_signature(obj):
        return re.sub('<(.*) at (.*)>', '<\\1>', str(inspect.signature(obj)))
    def format_definition(name, obj):
        argdef = '%s%s' % (name, format_signature(obj))
        return '\n\n*Definition:*\n        ``%s``\n' % argdef.strip()
    def add_str(s0, s1):
        return s0 + s1 + '\n'
    
    #G3Module documentation
    anti_recursion_protector = []
    def iterate_through_mod(mod, modname):
        mod_dict = {}
        for x, obj in mod.__dict__.items():
            if x.startswith('_') or x.endswith('_'):
                continue
            out_str = ''
            try:
                if hasattr(obj, '__module__'):
                    if hasattr(obj, '__wrapped__'):
                        modname = obj.__wrapped__.__module__
                    else:
                        modname = obj.__module__
                    itemname = '%s.%s' % (modname, obj.__name__)
                else:
                    itemname = obj.__name__ #helps with files imported in __init__
                if not itemname.startswith(module_path):
                    continue
                if itemname in anti_recursion_protector:
                    continue
                anti_recursion_protector.append(itemname)
            except:
                continue
            
            subclasstest = False
            try:
                subclasstest = issubclass(obj, G3Module)
            except TypeError:
                pass
            if hasattr(obj, '__g3module__') or subclasstest:
                out_str = add_str(out_str, format_name(modname, x))

                if format_doc(obj) is not None:
                    baredoc = format_doc(obj)
                    out_str = add_str(out_str, baredoc)
                try:
                    # in py3, all functions have an __init__, so check if we are dealing with a function first
                    # this logic should be backward-compatible with py2
                    if inspect.isfunction(obj):
                        out_str = add_str(out_str, format_definition(x, obj))
                    else:
                        out_str = add_str(out_str, '\n*Constructor:*\n\t``%s%s``\n' %
                                          ((x, format_signature(obj.__init__)) ) )
                        if format_doc(obj.__init__):
                            con_str = '\n*Constructor:*\n\t%s\n' % format_doc(obj.__init__).replace('\n', '\n\t')
                            con_str = con_str.replace(' -> None', ' -> None``')
                            con_str = con_str.replace('__init__', '``__init__')
                            out_str += con_str + '\n'
                except:
                    pass
                out_str = add_str(out_str, '')
                mod_dict[itemname] = out_str
            elif hasattr(obj, '__pipesegment__'):
                out_str = add_str(out_str, format_name(modname, x))
                if format_doc(obj) is not None:
                    if hasattr(obj, '__rstdoc__'):
                        out_str = add_str(out_str, obj.__rstdoc__)
                    else:
                        baredoc = format_doc(obj)
                        out_str = add_str(out_str, baredoc)
                else:
                    out_str = add_str(out_str, '\nNo documentation\n')
                out_str = add_str(out_str, format_definition(x, obj))
                mod_dict[itemname] = out_str
            elif hasattr(obj, '__dict__'):
                mod_dict.update(iterate_through_mod(obj, itemname))
        return mod_dict

    #Useful Function documentation    
    other_anti_recursion_protector = []
    def iterate_through_func(mod, modname):
        mod_dict = {}
        for x, obj in mod.__dict__.items():
            out_str = ''
            if x.startswith('_') or x.endswith('_'):
                continue
            try:
                if hasattr(obj, '__module__'):
                    itemname = '%s.%s' % (obj.__module__, obj.__name__)
                    modname = obj.__module__
                else:
                    itemname = obj.__name__ #helps with files imported in __init__
                if not itemname.startswith(module_path):
                    continue
                if itemname in other_anti_recursion_protector:
                    continue
                other_anti_recursion_protector.append(itemname)
            except:
                continue
            subclasstest = False
            try:
                subclasstest = 'Boost.Python.function' in str(type(obj))
            except TypeError:
                pass
            if hasattr(obj, '__g3usefulfunc__') or subclasstest:
                out_str = add_str(out_str, format_name(modname, x))
                if format_doc(obj) is not None:
                    if subclasstest:
                        tmp_str = format_doc(obj, simple=True)
                        out_str = out_str +tmp_str
                    else:
                        tmp_str = format_doc(obj)
                        out_str = out_str + tmp_str
                    if subclasstest:
                        out_str = out_str 
                try:
                    out_str = add_str(out_str, format_definition(x, obj))
                except:
                    pass
                out_str = add_str(out_str, '')
                mod_dict[itemname] = out_str
            elif hasattr(obj, '__dict__'):
                mod_dict.update(iterate_through_func(obj, itemname))
        return mod_dict

    #Frame Object documentation    
    other_other_anti_recursion_protector = []
    def iterate_through_frame_object(mod, modname):
        mod_dict = {}
        for x, obj in mod.__dict__.items():
            out_str = ''
            try:
                if hasattr(obj, '__module__'):
                    itemname = '%s.%s' % (obj.__module__, obj.__name__)
                    modname = obj.__module__
                else:
                    itemname = obj.__name__ #helps with files imported in __init__
                if not itemname.startswith(module_path):
                    continue
                if itemname in other_other_anti_recursion_protector:
                    continue
                other_other_anti_recursion_protector.append(itemname)
            except:
                continue
            subclasstest = False
            try:
                subclasstest = issubclass( obj, G3FrameObject)
            except TypeError:
                pass
            if hasattr(obj, '__g3frameobject__') or subclasstest:
                out_str = add_str(out_str, format_object(modname, x))
                if format_doc(obj) is not None:
                    tmp_str = format_doc(obj).strip()
                    out_str = out_str + tmp_str
                    if format_doc(obj.__init__):
                        con_str = '\n\n*Constructors:*\n\t%s\n' % format_doc(obj.__init__).replace('\n', '\n\t')
                        con_str = con_str.replace(' -> None', '``')
                        con_str = con_str.replace(' -> object', '``')
                        con_str = con_str.replace('__init__( (object)arg1,', '``{}('.format(obj.__name__))
                        con_str = con_str.replace('__init__( (object)arg1)', '``{}()'.format(obj.__name__))
                        con_str = con_str.replace('__init__( (object)arg1 [,', '``{}( ['.format(obj.__name__))
                        out_str += con_str + '\n'
                    #after we have gotten the documention, find the properties and load their documentation
                    prop_str = ''
                    for p_name, p_obj in obj.__dict__.items():
                        if isinstance( p_obj, property):
                            if p_name in ['value', '__g3frameobject__']:
                                continue
                            prop_str += '\n'
                            tmp_doc = format_doc(p_obj)
                            if tmp_doc is None:
                                tmp_doc = 'No Doc (Shame!)'
                            prop_str += '* **%s**: %s\n' % (p_name, tmp_doc.strip() )
                    if prop_str:
                        out_str += '\n\n*Members:*\n'
                        out_str += prop_str
                    meth_str = ''
                    for p_name, p_obj in obj.__dict__.items():
                        if p_name.startswith('_') or p_name.endswith('_'):
                            continue
                        if 'Boost.Python.function' in str(type(p_obj)):
                            meth_str = add_str(meth_str, format_name('%s.%s' % (modname, x), p_name))
                            tmp_doc = format_doc(p_obj, simple=True)
                            meth_str += tmp_doc.strip()
                            meth_str = add_str(meth_str, '')
                    if meth_str:
                        out_str += '\n\n*Methods:*\n'
                        out_str += meth_str
                out_str = add_str(out_str, '')
                mod_dict[itemname] = out_str
            elif hasattr(obj, '__dict__'):
                mod_dict.update(iterate_through_frame_object(obj, itemname))
        return mod_dict

    out_str = ''

    obj_dict = iterate_through_frame_object(mod, module_path)
    obj_lst = sorted(obj_dict)
    if len(obj_dict) > 0:
        title_str = 'Frame Objects in %s' % module_path
        out_str = add_str(out_str, title_str)
        out_str = add_str(out_str,'='*len(title_str))
        if include_link_list:
            for m in obj_lst:
                out_str = add_str(out_str,'* %s_\n' % m)
        doc_str = '\n'.join([obj_dict[k] for k in obj_lst])
        out_str = add_str(out_str,doc_str)

    mod_dict = iterate_through_mod(mod, module_path)
    mod_lst = sorted(mod_dict)
    if len(mod_dict) > 0:
        title_str = 'Modules in %s' % module_path
        out_str = add_str(out_str, title_str)
        out_str = add_str(out_str,'='*len(title_str))
        if include_link_list:
            for m in mod_lst:
                out_str = add_str(out_str,'* %s_\n' % m)
        doc_str = '\n'.join([mod_dict[k] for k in mod_lst])
        out_str = add_str(out_str,doc_str)

    #now fun!
    fun_dict = iterate_through_func(mod, module_path)
    fun_lst = sorted(fun_dict)
    if len(fun_dict) > 0:
        title_str = 'Functions in %s' % module_path
        out_str = add_str(out_str,title_str)
        out_str = add_str(out_str,'='*len(title_str))
        if include_link_list:
            for m in fun_lst:
                out_str = add_str(out_str,'* %s_\n' % m)
        doc_str = '\n'.join([fun_dict[k] for k in fun_lst])
        out_str = add_str(out_str,doc_str)
    return out_str, mod_lst, fun_lst, obj_lst
