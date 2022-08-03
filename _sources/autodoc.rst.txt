Auto-documentation
==================

The SPT3G software can automatically generate documentation on the modules and
functions in all directories in the repository in a variety of formats. To generate documentation in the default format (HTML), run ``make docs`` in your build directory. Note that this must be run after the software has been built and with PYTHONPATH set appropriately (i.e. after running ``env-shell.sh``).

Getting it documented
---------------------

To ensure that your functions and modules are documented properly, you need to
tell the document generator that you want it to parse them. In Python,
``G3Module`` objects will automatically be parsed. For functions and non-``G3Module``
inherited classes, you'll need to decorate them with ``@core.indexmod``,
``@core.pipesegment`` or ``@core.usefulfunc``, depending on the type of object to
be documented.

In C++, functions exported as ``bp::def()`` will automatically be documented.
Additionally, any modules exported with the ``EXPORT_G3MODULE`` macro will
automatically be documented as well.

Documentation needs to be valid RST_.  Improperly formatted RST may result in 
really weird html.  To be 100% sure that you are not generating RST warnings, 
you can run ``make docs`` in your build directory and check the output.

.. _RST: http://docutils.sourceforge.net/rst.html

Viewing the docs
----------------

``spt3g-inspect`` generates the documentation in an rst format. ``make docs``
will authomatically generate an HTML browseable doc. If a project includes a README.rst file at the root of its directory tree, the contents of this file will be prepended to the manual page for the project.
