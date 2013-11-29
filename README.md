part_frac
=========

A partial fraction decomposition algorithm coded in the Sage Computer Algebra System.

This is module of functions written in SAGE that can be used to find the
partial fraction decomposition of a proper or improper rational function.

All of the methods in this module are based on exercises in Joel S. Cohen's
book Computer Algebra And Symbolic Computation: Mathematical Methods.

This module contains various methods with the main purpose to be used by the
part_frac method which returns the partial fraction expansion of a rational
function in one indeterminant.

To run the file use, the following command in the SAGE shell.

sage: runfile 'part_frac.sage'

Then you are free to use the part_frac function. See the docstring for the
part_frac function for examples.

To run the file as a module, you need to preparse the module as a python file;
i.e., use the command in a linux shell (NOT in the SAGE shell):

$ sage --preparse part_frac.sage

This will create a file part_frac.py, which can then be imported as a module
like in normal python.

Next, just call the part_frac method with the desired inputs inside the SAGE
shell. For examples, see the docstrings of the various methods.

Because the purpose of this module is for educational purposes, I used methods
for polynomial division, polynomial GCD, etc., that I wrote instead of Sage's
built-in methods. It should be noted though that I do use Sage's factor
method.

Name:       John Kluesner
Date:       28 Nov, 2013
Email:      stringman45@gmail.com
