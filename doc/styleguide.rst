-----------
Style Guide
-----------

Version Control Hygiene
-----------------------

You can use two mechanisms to access the repository: git and SVN. The following is a brief overview of how to use these in a way that your collaborators will appreciate.

Git
===

To initially check out the repository:

.. code-block:: sh

	git clone https://user@github.com/CMB-S4/spt3g_software.git

To update your checkout (the --rebase is important, especially if you have local changes):

.. code-block:: sh

	git pull --rebase

To send your changes back:

.. code-block:: sh

	git diff files_to_commit <- Examine this
	git commit files_to_commit
	git push


SVN
===

To initially check out the repository:

.. code-block:: sh

	svn co https://user@github.com/CMB-S4/spt3g_software/trunk spt3g_software

To update your checkout:

.. code-block:: sh

	svn up

To send your changes back:

.. code-block:: sh

	svn diff files_to_commit <- Examine this
	svn ci files_to_commit

Coding Style
------------

We have style guidelines for both Python and C code. These are intended to keep code readable and avoid jarring transitions before code written by different people. These are summarized below. In general, the guiding principle is to make new code look like existing code. If you run into a file that obviously violates these style guidelines, please (a) hit the author with a stick and (b) do *not* change the style yourself. When adding to files that conform to some consistent, but different, style, prefer the style of the existing code to these guidelines unless you are replacing more than 50% of the file.

Python Style Guide
==================
For python we will be following PEP8_ (the style guide for the python standard library).
PEP8 is a fairly significant document, so I will summarize the important points.
Nonetheless, it is a generally good idea to read at least the section on readability_.
These are what we consider the most important points of PEP8:

* Function and variables names should be lower case, with words separated by underscores,
  when it improves readability.
* Class names should use CapitalizedWords.
* Constants should be in ALL_CAPS.
* Non-public class methods and instance variables should begin with a single underscore
  (this is baked into the way python implements classes).
* Lines should be short.  The standard convention is 79 characters.  We're not going to be anal
  about this, but long lines should not be significantly longer than that, as they are hard to read.
  For reference, most terminals open at 80 characters wide.
* Indentation should be done with 4 spaces.  At least one person is going to continue to violate
  this rule, so we should generally

For example:

.. code-block:: python

    CONSTANT_VALUE = 10

    def function_name(variable_name, other):
        do_some_things

    class MyClass:
        CONSTANCE_INSTANCE_VARIABLE = 100
        def __init__(self, input, other_input):
            self.instance_variable = input + other_input

        def class_method(self, input):
            do_things

For classes, modules, and functions individuals are HIGHLY encouraged to use the standard SPT3G docstring format. 
This includes individually defining each argument and output of a function along with their default values and types.
Below is an example of a well written docstring using the numpy_ docstring format.
Please adhere to this format for all spt3g_software code.

Docsting example:

.. code-block:: python 

    def get_fft_scale_fac(res=None, n1=None, n2=None,
                          apod_mask=None, parent=None):
        """
        Returns a scale factor to convert a 2d fft of a map into
        sqrt(C_l) normalized units by scaling it by the sky area.
        In particular it returns:
          np.mean(apod_mask**2)**0.5 * (n1 * n2)**0.5 / reso_rad

        For forward transforms: fft2(mp)  / scale_fac
        For inverse transforms: ifft2(mp) * scale_fac

        Arguments
        ---------
        res : float, default=None
            Resolution of the map in G3Units. 
            Must be specified if parent is None.
        n1 : int, default=None
            Number of map pixels in x direction. 
            Must be specified if parent is None.
        n2 : int, default=None
            Number of map pixels in y direction.
            Must be specified if parent is None.
        apod_mask : FlatSkyMap or ndarray, default=None
            An apodization mask.
        parent : FlatSkyMap, default=None
            Parent Map object containing resolution and shape info.
            Must be specified if res, n1, and n2 are None.
    
        Returns
        -------
        float
	    FFT normalization factor

        See Also
        -------
        ft_to_map
        """

.. _PEP8: https://www.python.org/dev/peps/pep-0008/
.. _readability: https://www.python.org/dev/peps/pep-0008/#a-foolish-consistency-is-the-hobgoblin-of-little-minds
.. _numpy: https://numpydoc.readthedocs.io/en/latest/format.html

C++ Style Guide
===============

C and C++ code in the tree is supposed to follow the KNF_ conventions. On Mac OS X systems, you have a local copy of this guide installed in a man page named ``style``. Like PEP8, the document is long, so a brief summary follows. The general rules are similar to PEP8 in most cases.

* Variable and global function names should be lower case, with words separated by underscores when it improves readability.
* Class names and class member function should use CapitalizedWords.
* Constants should be in ALL_CAPS, as should all macros created with ``#define``.
* Non-public class variables should end with a single underscore.
* Indentation for control flow marking (e.g. inside an ``if``) should be done with hard tabs. Indentation for alignment should be done with spaces following the number of hard tabs indicated for control flow. This allows control flow indentation to be uniform no matter what the setting for tab width is on people's editors.
* Lines should be short.  The standard convention is 79 characters.  We're not going to be anal
  about this, but long lines should not be significantly longer than that, as they are hard to read.
  For reference, most terminals open at 80 characters wide. If you need to continue a line, it should be indented four spaces deeper than the beginning of the line. For line widths, tabs are considered to represent 8 spaces.
* For functions, the opening brace (``{``) should be on the line *following* the function definition. For classes, conditionals, and loops, it should be on the *same* line.
* Where functions are not used outside of a single C++ file, they should be marked ``static``.

Conforming code looks like this:

.. code-block:: c++

	const int CONSTANT_VALUE = 10;

	void
	function_name(int variable_name, double other)
	{
		do_some_things(variable_name);
	}

	class MyClass : public Parent {
	public:
		MyClass(int input, int other_input);
		void ClassMethod(int input);
	private:
		int instance_variable_;
	};

	MyClass::MyClass(int input, int other_input) : Parent(input, 15, 12,
	    "long line")
	{
		instance_variable_ = input + other_input;
	}

	void
	MyClass::ClassMethod(int input)
	{
		do_things();
	}

.. _KNF: https://www.freebsd.org/cgi/man.cgi?query=style&sektion=9

