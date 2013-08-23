vasp2cif
========

A program that can read VASP output files and generate crystal structure files in CIF file format. The CIF files can then be used for visualization in graphical applications. Note that vasp2cif will not preserve symmetries from VASP – the output CIF is always in P1 symmetry. If you have the [FINDSYM](http://stokes.byu.edu/iso/isotropy.php) program by Harold Stokes et al. installed, vasp2cif will attempt to use that to find any symmetries.

Installation
------------

1. Download or clone the repository.
2. Install by putting the vasp2cif.py in a folder in your $PATH (such as ~/bin).
3. Change permissions: ```mv vasp2cif.py vasp2cif; chmod u+x vasp2cif```.

This should allow you to write e.g. "vasp2cif CONTCAR" in any job folder and get the corresponding CONTCAR.cif file.

Basic usage
-----------

Convert CONTCAR to CIF in a job directory:

    vasp2cif CONTCAR

Convert a lone POSCAR without a corresponding POTCAR file to CIF:

    vasp2cif --elements="Zn,O" POSCAR

Generate CIF file with multiple data blocks from a relaxation run:

    vasp2cif OUTCAR

Authors
-------

* [Peter Larsson](http://www.nsc.liu.se/~pla/)
* [Torbjörn Björkman](http://physics.aalto.fi/personnel/?id=538)

License
---------------------

Copyright 2008-2013 Peter Larsson under the Apache License 2.0. See the LICENSE file for more detailed information.
