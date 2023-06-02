# Compilation

The code can be compiled using the included makefile.

# Input

Input files must be in exp format, see the file 333-27-mod2.exp or the schemes directory for an
example.

# Usage

./flip filename.exp n m p pathlength restart

If n=m=p then m and p can be omitted.

The parameter pathlength sets a limit for the number of flips before aborting. The parameter restart
can be set to 0 or 1, Setting it to 0 will produce a scheme of rank one less than the input scheme,
setting it to 1 will restart the search procedure everytime a reduction is found, with the newly
found scheme as input until no reduction is found within the given limit on the pathlength.

# Example

Try:

./flip 333-27-mod2.exp 3 1000000 1
