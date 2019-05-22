import sys, inspect, re
from spt3g.core import G3Module, G3FrameObject

def format_doc(x):
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
    
    
    def format_name(modname, x):
        return '\n.. _%s.%s:\n\n**%s.%s**\n' % (modname, x, modname, x)
    def format_arguments(s0,s1):
        argdef = '%s%s'%(s0,s1)
        return '\n\n*Definition:*\n        ``%s``\n' % argdef.strip()
    def add_str(s0, s1):
        return s0 + s1 + '\n'
    
    #G3Module documentation
    anti_recursion_protector = []
    def iterate_through_mod(mod, modname):
        out_str = ''
        mod_lst = []
        for x in sorted(mod.__dict__.keys()):
            try:
                if hasattr(mod.__dict__[x], '__module__'):
                    itemname = '%s.%s' % (mod.__dict__[x].__module__, mod.__dict__[x].__name__)
                else:
                    itemname = mod.__dict__[x].__name__ #helps with files imported in __init__
                if not itemname.startswith(module_path):
                    continue
                if itemname in anti_recursion_protector:
                    continue
                anti_recursion_protector.append(itemname)
            except:
                continue
            
            subclasstest = False
            try:
                subclasstest = issubclass(mod.__dict__[x], G3Module)
            except TypeError:
                pass
            if hasattr(mod.__dict__[x], '__g3module__') or subclasstest:
                out_str = add_str(out_str, format_name(modname, x))
                mod_lst.append('%s.%s'%(modname, x))

    
                if format_doc(mod.__dict__[x]) is not None:
                    baredoc = format_doc(mod.__dict__[x])
                    out_str = add_str(out_str, ' '.join(re.split('\n[ \t]+', baredoc))) # Join new lines starting with whitespace with the previous line
                try:
                    # in py3, all functions have an __init__, so check if we are dealing with a function first
                    # this logic should be backward-compatible with py2
                    if inspect.isfunction(mod.__dict__[x]):
                        out_str = add_str(out_str, format_arguments(x, inspect.formatargspec(*inspect.getargspec(mod.__dict__[x]))))
                    else:
                        out_str = add_str(out_str, '\n*Constructor:*\n\t``%s%s``\n' %
                                          ((x, inspect.formatargspec(*inspect.getargspec(mod.__dict__[x].__init__))) ) )
                        if format_doc(mod.__dict__[x].__init__):
                            con_str = '\n*Constructor:*\n\t%s\n' % format_doc(mod.__dict__[x].__init__).replace('\n', '\n\t')
                            con_str = con_str.replace(' -> None', ' -> None``')
                            con_str = con_str.replace('__init__', '``__init__')
                            out_str += con_str + '\n'
                except:
                    pass
                out_str = add_str(out_str, '')
            elif hasattr(mod.__dict__[x], '__pipesegment__'):
                out_str = add_str(out_str, format_name(modname, x))
                mod_lst.append('%s.%s'%(modname, x))
                if format_doc(mod.__dict__[x]) is not None:
                    if hasattr(mod.__dict__[x], '__rstdoc__'):
                        out_str = add_str(out_str, mod.__dict__[x].__rstdoc__)
                    else:
                        baredoc = format_doc(mod.__dict__[x])
                        out_str = add_str(out_str, ' '.join(re.split('\n[ \t]+', baredoc))) # Join new lines starting with whitespace with the previous line
                else:
                    out_str = add_str(out_str, '\nNo documentation\n')
                out_str = add_str(out_str, format_arguments(x, inspect.formatargspec(*inspect.getargspec(mod.__dict__[x]))))
            elif hasattr(mod.__dict__[x], '__dict__'):
                out_str_tmp, mod_lst_tmp = iterate_through_mod(mod.__dict__[x], modname + '.' + x)
                out_str += out_str_tmp
                mod_lst += mod_lst_tmp
        return out_str, mod_lst

    #Useful Function documentation    
    other_anti_recursion_protector = []
    def iterate_through_func(mod, modname):
        out_str = ''
        mod_lst = []
        for x in sorted(mod.__dict__.keys()):
            if x.startswith('_'):
                continue
            try:
                if hasattr(mod.__dict__[x], '__module__'):
                    itemname = '%s.%s' % (mod.__dict__[x].__module__, mod.__dict__[x].__name__)
                else:
                    itemname = mod.__dict__[x].__name__ #helps with files imported in __init__
                if not itemname.startswith(module_path):
                    continue
                if itemname in other_anti_recursion_protector:
                    continue
                other_anti_recursion_protector.append(itemname)
            except:
                continue
            subclasstest = False
            try:
                subclasstest = str(type(mod.__dict__[x])) == "<type 'Boost.Python.function'>"
            except TypeError:
                pass
            if hasattr(mod.__dict__[x], '__g3usefulfunc__') or subclasstest:
                out_str = add_str(out_str, format_name(modname, x))
                mod_lst.append('%s.%s'%(modname, x))
                if format_doc(mod.__dict__[x]) is not None:
                    if subclasstest:
                        tmp_str = '``' + format_doc(mod.__dict__[x]).strip().replace(':\n',':``\n')
                        if (tmp_str.count('``') <2):
                            tmp_str += '``'
                        out_str = out_str +tmp_str
                    else:
                        tmp_str = format_doc(mod.__dict__[x])
                        out_str = out_str + tmp_str
                    if subclasstest:
                        out_str = out_str 
                try:
                    out_str = add_str(out_str, format_arguments(x, inspect.formatargspec(*inspect.getargspec(mod.__dict__[x]))))
                except:
                    pass
                out_str = add_str(out_str, '')
            elif hasattr(mod.__dict__[x], '__dict__'):
                out_str_tmp, mod_lst_tmp = iterate_through_func(mod.__dict__[x], modname + '.' + x)
                out_str += out_str_tmp
                mod_lst += mod_lst_tmp
        return out_str, mod_lst

    #Frame Object documentation    
    other_other_anti_recursion_protector = []
    def iterate_through_frame_object(mod, modname):
        out_str = ''
        mod_lst = []
        for x in sorted(mod.__dict__.keys()):
            try:
                if hasattr(mod.__dict__[x], '__module__'):
                    itemname = '%s.%s' % (mod.__dict__[x].__module__, mod.__dict__[x].__name__)
                else:
                    itemname = mod.__dict__[x].__name__ #helps with files imported in __init__
                if not itemname.startswith(module_path):
                    continue
                if itemname in other_other_anti_recursion_protector:
                    continue
                other_other_anti_recursion_protector.append(itemname)
            except:
                continue
            subclasstest = False
            try:
                subclasstest = issubclass( mod.__dict__[x], G3FrameObject)
            except TypeError:
                pass
            if subclasstest:
                out_str = add_str(out_str, format_name(modname, x))
                mod_lst.append('%s.%s'%(modname, x))
                if format_doc(mod.__dict__[x]) is not None:
                    tmp_str = format_doc(mod.__dict__[x]).strip()
                    out_str = out_str + tmp_str
                    #after we have gotten the documention, find the properties and load their documentation
                    out_str += '\n\n    Members:\n'
                    for p_name, p_obj in mod.__dict__[x].__dict__.items():
                        if isinstance( p_obj, property):
                            if p_name == 'value':
                                continue
                            out_str += '\n'
                            tmp_doc = format_doc(p_obj)
                            if tmp_doc is None:
                                tmp_doc = 'No Doc (Shame!)'
                            out_str += '    * **%s**: %s\n' % (p_name, tmp_doc.strip() )
                    
                out_str = add_str(out_str, '')
            elif hasattr(mod.__dict__[x], '__dict__'):
                out_str_tmp, mod_lst_tmp = iterate_through_frame_object(mod.__dict__[x], modname + '.' + x)
                out_str += out_str_tmp
                mod_lst += mod_lst_tmp
        return out_str, mod_lst

    out_str = ''


    doc_str, obj_lst = iterate_through_frame_object(mod, module_path)
    if len(doc_str.strip()) > 0:
        title_str = 'Frame Objects in %s' % module_path
        out_str = add_str(out_str, title_str)
        out_str = add_str(out_str,'='*len(title_str))
        if include_link_list:
            for m in obj_lst:
                out_str = add_str(out_str,'* %s_\n' % m)
        out_str = add_str(out_str,doc_str)



    doc_str, mod_lst = iterate_through_mod(mod, module_path)
    if len(doc_str.strip()) > 0:
        title_str = 'Modules in %s' % module_path
        out_str = add_str(out_str, title_str)
        out_str = add_str(out_str,'='*len(title_str))
        if include_link_list:
            for m in mod_lst:
                out_str = add_str(out_str,'* %s_\n' % m)
        out_str = add_str(out_str,doc_str)

    #now fun!
    doc_str, fun_lst = iterate_through_func(mod, module_path)

    if len(doc_str.strip()) > 0:
        title_str = 'Functions in %s' % module_path
        out_str = add_str(out_str,title_str)
        out_str = add_str(out_str,'='*len(title_str))
        if include_link_list:
            for m in fun_lst:
                out_str = add_str(out_str,'* %s_\n' % m)
        out_str = add_str(out_str,doc_str)
    return out_str, mod_lst, fun_lst, obj_lst
