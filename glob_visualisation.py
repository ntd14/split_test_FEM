from dolfin import *
import numpy as np
import math as math
import sys as sys
import pandas as pandas
from dolfin.cpp._mesh import Cell_get_cell_data, Cell_get_vertex_coordinates,\
    Cell_normal, Cell_cell_normal, Cell_contains

mesh = Mesh("./final_xml/glob1.xml")

df = pandas.read_csv('glob1_properties.csv', sep=',')

print(df)
print(df['log_ID'])

BED = df['BED'][0]
SED = df['SED'][0]
slit_length = df['slit_length'][0]
opening = df['opening'][0]
log_len = df['log_len'][0]

# need to some how set the stress/strain at each point we have, 
# then interpolate that over the whole mesh? setting a netural 
# plane and then calculating the internal stresses based on sum to zero

#also probably try to make it for a generic number of holes, an have an import for a file with the mesh data in it

#Later will need to port it all to a new mesh that has the cut in it

plot(mesh, interactive = True)