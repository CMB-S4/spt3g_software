def usefulfunc(func):
    '''
    Mark argument as a useful function that can be found by automated documentation tools. 

    For example:

    @core.usefulfunc
    def do_some_science(data):
        science(data)
    '''

    func.__g3usefulfunc__ = True

    return func

