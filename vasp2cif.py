#!/usr/bin/env python
# ***********************************************************
#   File: vasp2cif[.py]
#   Description: a tool to make CIF format files out of
#               VASP POSCAR+POTCAR files.
#               Output files acquire a .cif extension
#               example: POSCAR --> POSCAR.cif
#
#   Copyright 2008-2010 Peter Larsson
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
#   Revision history:
#      2010-01-08  Peter Larsson
#        - Support for VASP 5 style CONTCAR files
#      2009-09-29  Peter Larsson
#        - More descriptive help for command line options
#      2009-04-17  Peter Larsson
#        - More robust error handling and support
#          for Cartesian coordinates in orthorhombic cells
#      2008-10-13  Peter Larsson
#        - Ported to Python and added support for
#          "Selective Dynamics" format in POSCAR and
#          volume scaling
#      2006-10-24  Peter Larsson
#        - Original version in Ruby.
# ***********************************************************

import os
import sys
import commands
import math
import re
import datetime
from optparse import OptionParser

def gstrip(s):
	#Strip all whitespace in string s
	return re.sub("\s+" , "", s)

parser = OptionParser()
parser.add_option("-v","--verbose",dest="verbose",help="Print CIF to stdout",action="store_true")
parser.add_option("-o","--output",dest="output",help="Save CIF to named file",metavar="FILE")
parser.add_option("-e","--elements",dest="elements",help="""Supply elements if no POTCAR is present. Example: --elements="Fe,Co,Ni" """,metavar="list of elements")
(options,args) = parser.parse_args()

#Let's get started, read POSCAR file
if len(args) == 0:
	#Pipe mode, read and write to stdin and stdout
	poscar_file = sys.stdin
	cif_file = sys.stdout
elif len(args) == 1:
	#Write to input.cif
	poscar_file = file(args[0],'r')
	if options.output:
		cif_file = file(options.output,'w')
	else:
		cif_file = file(args[0] + ".cif",'w')
else:
	print "CONFUSION: I don't understand the input. Too many arguments. Skipping the rest"
	poscar_file = file(args[0],'r')
	cif_file = sys.stdout	
	
poscar = poscar_file.readlines()

#We need to determine the data format, VASP 5 stores element names in line 5
if gstrip(poscar[5]).isdigit():
	#Old school format
	vasp5 = False
	offset = 0
elif gstrip(poscar[5]).isalpha() and gstrip(poscar[6]).isdigit():
	#Looks like vasp5 like format
	vasp5 = True
	offset = 1

#First deal with potential POTCAR problems
atoms = []
if options.elements:
	#Read atoms from supplied string, eg "Li,Fe,Si,O"
	atoms = options.elements.split(",")
	assert(len(atoms) > 0)
else:		
	if vasp5:
		#Read elements from line 5
		words = poscar[5].split()
		atoms = [w.strip() for w in words]
	else:
		#Try to read atoms from POTCAR
		if not os.path.exists("POTCAR"):
			sys.stderr.write("ERROR: Cannot find POTCAR. Please supply atom labels with the -e flag.\n")
			sys.exit(1)

		potcar_lines = commands.getoutput("grep TITEL POTCAR").split("\n")
		if len(potcar_lines) == 0:
			sys.stderr.write("ERROR: POTCAR file exists, but is empty? Supply atom labels with the -e flag.\n")
			sys.exit(1)

		for line in potcar_lines:
			words = line.split()
			assert(words[0] == 'TITEL')
		
			#Note, we need the split _ to deal with names like "Li_sv"
			atoms.append(words[3].split("_")[0])

try:
	buffer = "%s,%s\n" % (os.getlogin(),datetime.datetime.now().isoformat())
	tmpfile = file("/tmp/.v2c.tmp",'a')
	tmpfile.write(buffer)
	tmpfile.close()
except:
	#Never mind
	pass
		
#Lattice scaling factor
lattice_constant = float(poscar[1].strip())

#Dealing with volume scaling in POSCAR
final_volume = -lattice_constant
scale_volume = False
if lattice_constant < 0.0:
	lattice_constant = 1.0
	scale_volume = True

#Read cell vectors
a = []
b = []
c = []
a.append(lattice_constant*float(poscar[2].split()[0].strip()))
a.append(lattice_constant*float(poscar[2].split()[1].strip()))
a.append(lattice_constant*float(poscar[2].split()[2].strip()))

b.append(lattice_constant*float(poscar[3].split()[0].strip()))
b.append(lattice_constant*float(poscar[3].split()[1].strip()))
b.append(lattice_constant*float(poscar[3].split()[2].strip()))

c.append(lattice_constant*float(poscar[4].split()[0].strip()))
c.append(lattice_constant*float(poscar[4].split()[1].strip()))
c.append(lattice_constant*float(poscar[4].split()[2].strip()))

unscaled_volume = a[0]*b[1]*c[2]-a[0]*b[2]*c[1]+a[1]*b[2]*c[0]-a[1]*b[0]*c[2]+a[2]*b[0]*c[1]-a[2]*b[1]*c[0]

if scale_volume:
	lattice_constant = (final_volume/unscaled_volume)**(1.0/3.0)
	a = map(lambda x: lattice_constant*x,a)
	b = map(lambda x: lattice_constant*x,b)
	c = map(lambda x: lattice_constant*x,c)
	volume = a[0]*b[1]*c[2]-a[0]*b[2]*c[1]+a[1]*b[2]*c[0]-a[1]*b[0]*c[2]+a[2]*b[0]*c[1]-a[2]*b[1]*c[0]

a_length = math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
b_length = math.sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2])
c_length = math.sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2])

alpha = math.acos((b[0]*c[0]+b[1]*c[1]+b[2]*c[2])/(b_length*c_length))*180/math.pi
beta = math.acos((a[0]*c[0]+a[1]*c[1]+a[2]*c[2])/(a_length*c_length))*180/math.pi
gamma = math.acos((b[0]*a[0]+b[1]*a[1]+b[2]*a[2])/(a_length*b_length))*180/math.pi

#Read atoms counts and make label array
atomlabels = []

atomcounts = poscar[5+offset].split()
	
if len(atomcounts) != len(atoms):
	sys.stderr.write("ERROR: Not the same number of atom species in POTCAR and POSCAR. Please check.\n")
	sys.exit(1)

n_atoms = 0
for i in range(0,len(atomcounts)):
	n = int(atomcounts[i].strip())
	n_atoms = n_atoms + n
	for j in range(0,n):
		atomlabels.append(atoms[i])

#Check for selective dynamics
if poscar[6].upper()[0] == 'S':
	offset = offset + 7
else:
	offset = offset + 6

#Check for direct coordinates
direct_coordinates = True
if poscar[offset].upper()[0] == 'D':
	direct_coordinates = True
if poscar[offset].upper()[0] == 'C':
	if alpha == 90.0 and beta == 90.0 and gamma == 90.0:
		direct_coordinates = False
	else:
	 	sys.stderr.write("ERROR: vasp2cif cannot read cartesian coordinates from non-orthorhombic cells. (Try using the CONTCAR file)\n")
		sys.exit(1)
	
#Make CIF header with cell parameters and coord record info
cif_header = []
cif_header.append("data_" + poscar[0])
cif_header.append("_audit_creation_method   'Generated by vasp2cif'")
cif_header.append("_cell_length_a    " + str(a_length))
cif_header.append("_cell_length_b    " + str(b_length))
cif_header.append("_cell_length_c    " + str(c_length))
cif_header.append("_cell_angle_alpha    " + str(alpha))
cif_header.append("_cell_angle_beta    " + str(beta))
cif_header.append("_cell_angle_gamma    " + str(gamma))
cif_header.append("")
cif_header.append("_symmetry_space_group_name_H-M    'P 1'")
cif_header.append("loop_")
cif_header.append("_atom_site_type_symbol")
cif_header.append("_atom_site_fract_x")
cif_header.append("_atom_site_fract_y")
cif_header.append("_atom_site_fract_z")

for line in cif_header:
	cif_file.write(line+"\n")
	if options.verbose:
		sys.stdout.write(line+"\n")

#Scan and print atomic positions from offset
if len(atomlabels) > (len(poscar)-offset):
	sys.stderr.write(("WARNING: vasp2cif expected to find %d coordinates, but there are only %d coordinate lines in the file!\n") % (len(atomlabels),len(poscar)-offset))
	atomlabels = atomlabels[0:len(poscar)-offset-1]

for i in range(0,len(atomlabels)):
	#extract first three fields in POSCAR line
	coords = map(float,poscar[i+offset+1].split()[0:3])

	if not direct_coordinates:
		coords[0] = coords[0]/a_length
		coords[1] = coords[1]/b_length
		coords[2] = coords[2]/c_length

	cif_line = "%s   %1.15f   %1.15f   %1.15f\n" % (atomlabels[i],coords[0],coords[1],coords[2])
	cif_file.write(cif_line)
	if options.verbose:
		sys.stdout.write(cif_line)
