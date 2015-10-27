from dolfin import *
import numpy as np
import math as math
import sys as sys
import pandas as pandas
from dolfin.cpp._mesh import Cell_get_cell_data, Cell_get_vertex_coordinates,\
    Cell_normal, Cell_cell_normal, Cell_contains
from cmath import pi

mesh = Mesh("./final_xml/glob1.xml")

df = pandas.read_csv('glob1_properties.csv', sep=',')

BED = float(df['BED'][0])
SED = float(df['SED'][0])
slit_length = float(df['slit_length'][0])
opening = float(df['opening'][0])
log_len = float(df['log_len'][0])
split_strain = opening*BED/(1.74*slit_length**2)/1000000
El = 10000000000.0
split_stress = El*split_strain
df['rad_pos'] = (df['rad_pos']-1)*(pi/2)
df['cirad_strain'] = 11.6*df['cirad']/1000
df['cirad_stress'] = El*df['cirad_strain']

#function to convert cylindrical to cart and visaversa

# need to some how set the stress/strain at each point we have, 
# then interpolate that over the whole mesh? setting a netural 
# plane and then calculating the internal stresses based on sum to zero

#also probably try to make it for a generic number of holes, an have an import for a file with the mesh data in it

#Later will need to port it all to a new mesh that has the cut in it

plot(mesh, interactive = True)