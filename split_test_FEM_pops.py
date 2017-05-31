from dolfin import *
import numpy as np
import math as math
import sys as sys
from dolfin.cpp._mesh import Cell_get_cell_data, Cell_get_vertex_coordinates, Cell_normal, Cell_cell_normal, Cell_contains

mesh_path = sys.argv[1]
theta = float(sys.argv[2])
cutp = float(sys.argv[3])
cutn = float(sys.argv[4])
seed_num = int(sys.argv[5])
stress_l = float(sys.argv[6])

r_devider = -5.0 #  jacobs 1945 values in psi come to about 10 - stresses and strains in tree tunks as they grow
t_devider = -5.0 # stress_l = 5.23* stress_transverse  -- Patterns of longitudinal and tangential maturation stresses in Eucalyptus nitens plantation trees -- find boyd 1950 it has a lot more species in it also need to find Emods

#stress_l = 19439510 #for argo samples, + is the residual strain which is left in the stem,    #2014 = 7650562 #2012 = 1985117 #2013 -> 15928235; jacobs stresses and strian in tree trunks as they grow in length and width found the same pith stress value (ie expantion in the pith of 2800psi)
stress_sd_l = 8353013 # for argo samples #2014 = 6771055 #2012 = 4573265 #2013 -> 7700536
#define these if l and r/t are realted by the means etc their random distrobutions rather than by the specific l value for the cell
#stress_t = stress_l/-3.0
#stress_sd_t = stress_sd_l/3.0
#stress_r = stress_l/-1000.0
#stress_sd_r = stress_sd_l/1000.0

#assumed parameters about samples, this is what the meshes are
big_end_height = 400 #assumed small end centred on origin
rad_slices = 8
vert_slices = 4
num_of_radial_devisions = 3
big_rad = 14.5 # diameter is 29
small_rad = 12.5
#slit length = 200

# small rad = 12.5, big rad = 14.5
# mind cut = 8.33, 10.25

num_of_sd = num_of_radial_devisions * rad_slices * vert_slices
mesh = Mesh(mesh_path)
np.random.seed(seed_num)


norm_l = np.zeros(num_of_sd)
l_vals = np.random.normal(stress_l, stress_sd_l, 1)
for ii in range(0, num_of_sd):
	if ii%num_of_radial_devisions == 0:
		norm_l[ii] = -1*np.random.normal(stress_l, stress_sd_l, 1)[0]/2.0 #as what we measure with the splitting test is the stress gradiant, the pith expands and the pirphery contracts in the test, below we generate a pith value, then reverse the sign for piriphery, hence the full splitting test value is represented from half of the input. 
	elif ii%num_of_radial_devisions == 2:
		norm_l[ii] = -1*norm_l[ii - 2]
	elif ii%num_of_radial_devisions != 1:
		print('check "num_of_radial_devisions" parameter is devisable by 3')

norm_t = norm_l/t_devider
norm_r = norm_l/r_devider



def irads(max_rad):
    iR = float(max_rad)/math.sqrt(num_of_radial_devisions)
    iA = math.pi*iR**2
    rad_array = np.zeros(num_of_radial_devisions-1)
    for ii in range(1,num_of_radial_devisions):
        rad_array[ii-1] = math.sqrt(max_rad**2-(iA*ii)/math.pi)
    return rad_array[::-1]


def calc_rad(h, small_end_rad, big_end_rad):
    grad = big_end_height/(small_end_rad-big_end_rad)
    intersept = grad*small_end_rad
    return (h + intersept)/grad

#create subdomains

z_bound_array = np.zeros(vert_slices+1)
for ii in range(1,vert_slices+1):
    z_bound_array[ii] = ii*float(big_end_height)/vert_slices
    z_bound_array[vert_slices] = z_bound_array[vert_slices] + 0.001

off_set_array = np.zeros(vert_slices+1)
for ii in range(0, vert_slices+1):
    off_set_array[ii] = ii*((2*math.pi/rad_slices)/(vert_slices+1))

def tba(cur_height):
    for ii in range(0, vert_slices):
        if cur_height >= z_bound_array[ii] and cur_height < z_bound_array[ii+1]:
            chind = ii
            break

    thi_bound_array = np.zeros(rad_slices+1)
    thi_bound_array[0] = 0
    for ii in range(1,rad_slices+1):
        thi_bound_array[ii] = thi_bound_array[ii-1]+2*math.pi/rad_slices 
    thi_bound_array = thi_bound_array + off_set_array[chind]
    return(thi_bound_array)

def rba(cur_height):
    crad = calc_rad(cur_height, small_rad, big_rad)
    rad_array = irads(crad)
    r_bound_array = np.zeros(num_of_radial_devisions+1)
    r_bound_array[0] = 0
    for ii in range(1,len(rad_array)+1):
        r_bound_array[ii] = rad_array[ii-1]
    r_bound_array[len(r_bound_array)-1] = 100000
    return r_bound_array

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}
   

class trans_angles(Expression):
    def eval(self, values, x):
        tol = 10**(-10)
        if x[0] > -tol and x[0] < tol:
            if x[1] > tol:
                w = math.pi
            if x[1] >= -tol and x[1] <= tol:
                w = 0
            if x[1] < tol:
                w = -math.pi
        else:
            w = math.atan(x[1]/x[0])
            if x[0]<0 and x[1]>=0:
                w = w+math.pi
            if x[0]<0 and x[1]<0:
                w = w+math.pi
            if x[0]>0 and x[1]<0:
                w = w+2*math.pi
        
        values[0] = math.cos(w)
        values[1] = math.cos(w - math.pi/2)    
        #constant as rotation about z axis so always perpendicular
        values[2] = 0# math.cos(math.pi/2)
        # roation angle from y to r print CMM
        values[3] = math.cos(w + math.pi/2)
        values[4] = math.cos(w)
        values[5] = 0#math.cos(math.pi/2)
        values[6] = 0#math.cos(math.pi/2)
        values[7] = 0#math.cos(math.pi/2)
        values[8] = 1.0#math.cos(0.0)

    def value_shape(self):
        return (9,)
trans = trans_angles()

def end_boundary(x, on_boundary):
    return abs(x[2]) < 0.001

#Elastic constants of wood determined by ultrasound using three geometries of specimens cover all 9 elastic constants for E. siligna
E3 = 11.33*1000000000
E2 = E3/6.5
E1 = E3/4.3
V21 = 0.36##Swapping these goes from +x, -y to +x +y.
V12 = V21*E1/E2##
V31 = 0.36#### no change
V13 = V31*E1/E3####
V32 = 0.56##Small change near 0,0
V23 = V32*E2/E3##
Ge23 = E2
Ge31 = Ge23/10.0
Ge12 = Ge23/4.0

cA11 = trans[0]
cA12 = trans[1]
cA13 = 0
cA21 = trans[3]
cA22 = trans[4]
cA23 = 0
cA31 = 0
cA32 = 0
cA33 = 1

#stiffness matrix entries
CM11 = E1 - E1*V12*V21/(E2*(V12*V21/E2 - 1/E2)) - (V13 - V12*(V13*V21/E2 + V23/E2)/(V12*V21/E2 - 1/E2))*(E1*V21*(V12*V31/E3 + V32/E3)/(E2*(V12*V21/E2 - 1/E2)) - E1*V31/E3)/((V13*V21/E2 + V23/E2)*(V12*V31/E3 + V32/E3)/(V12*V21/E2 - 1/E2) - V13*V31/E3 + 1/E3) 
CM12 = -E1*V21/(E2*(V12*V21/E2 - 1/E2)) + (V13*V21/E2 + V23/E2)*(E1*V21*(V12*V31/E3 + V32/E3)/(E2*(V12*V21/E2 - 1/E2)) - E1*V31/E3)/((V12*V21/E2 - 1/E2)*((V13*V21/E2 + V23/E2)*(V12*V31/E3 + V32/E3)/(V12*V21/E2 - 1/E2) - V13*V31/E3 + 1/E3)) 
CM13 = -(E1*V21*(V12*V31/E3 + V32/E3)/(E2*(V12*V21/E2 - 1/E2)) - E1*V31/E3)/((V13*V21/E2 + V23/E2)*(V12*V31/E3 + V32/E3)/(V12*V21/E2 - 1/E2) - V13*V31/E3 + 1/E3) 
CM14 = 0.0
CM15 = 0.0
CM16 = 0.0
CM21 = -V12/(V12*V21/E2 - 1/E2) - (V13 - V12*(V13*V21/E2 + V23/E2)/(V12*V21/E2 - 1/E2))*(V12*V31/E3 + V32/E3)/((V12*V21/E2 - 1/E2)*((V13*V21/E2 + V23/E2)*(V12*V31/E3 + V32/E3)/(V12*V21/E2 - 1/E2) - V13*V31/E3 + 1/E3)) 
CM22 = -1/(V12*V21/E2 - 1/E2) + (V13*V21/E2 + V23/E2)*(V12*V31/E3 + V32/E3)/((V12*V21/E2 - 1/E2)**2*((V13*V21/E2 + V23/E2)*(V12*V31/E3 + V32/E3)/(V12*V21/E2 - 1/E2) - V13*V31/E3 + 1/E3))
CM23 = -(V12*V31/E3 + V32/E3)/((V12*V21/E2 - 1/E2)*((V13*V21/E2 + V23/E2)*(V12*V31/E3 + V32/E3)/(V12*V21/E2 - 1/E2) - V13*V31/E3 + 1/E3)) 
CM24 = 0.0
CM25 = 0.0
CM26 = 0.0
CM31 = (V13 - V12*(V13*V21/E2 + V23/E2)/(V12*V21/E2 - 1/E2))/((V13*V21/E2 + V23/E2)*(V12*V31/E3 + V32/E3)/(V12*V21/E2 - 1/E2) - V13*V31/E3 + 1/E3)
CM32 = -(V13*V21/E2 + V23/E2)/((V12*V21/E2 - 1/E2)*((V13*V21/E2 + V23/E2)*(V12*V31/E3 + V32/E3)/(V12*V21/E2 - 1/E2) - V13*V31/E3 + 1/E3)) 
CM33 = 1/((V13*V21/E2 + V23/E2)*(V12*V31/E3 + V32/E3)/(V12*V21/E2 - 1/E2) - V13*V31/E3 + 1/E3) 
CM34 = 0.0
CM35 = 0.0
CM36 = 0.0
CM41 = 0.0
CM42 = 0.0
CM43 = 0.0
CM44 = (2*Ge23)
CM45 = 0.0
CM46 = 0.0
CM51 = 0.0
CM52 = 0.0
CM53 = 0.0
CM54 = 0.0
CM55 = (2*Ge31)
CM56 = 0.0
CM61 = 0.0
CM62 = 0.0
CM63 = 0.0
CM64 = 0.0
CM65 = 0.0
CM66 = (2*Ge12)

#assembling stiffness matrix
CMM = as_matrix([[CM11, CM12, CM13, CM14, CM15, CM16], [CM21, CM22, CM23, CM24, CM25, CM26], [CM31, CM32, CM33, CM34, CM35, CM36], \
                 [CM41, CM42, CM43, CM44, CM45, CM46], [CM51, CM52, CM53, CM54, CM55, CM56], [CM61, CM62, CM63, CM64, CM65, CM66]])

G11 = cA11**2
G12 = cA12**2
G13 = cA13**2
G14 = cA12*cA13
G15 = cA11*cA13
G16 = cA11*cA12
G21 = cA21**2
G22 = cA22**2
G23 = cA23**2
G24 = cA22*cA23
G25 = cA21*cA23
G26 = cA21*cA22
G31 = cA31**2
G32 = cA32**2
G33 = cA33**2
G34 = cA32*cA33
G35 = cA31*cA33
G36 = cA31*cA32
G41 = 2*cA21*cA31
G42 = 2*cA22*cA32
G43 = 2*cA23*cA33
G44 = cA22*cA33+cA23*cA32 
G45 = cA21*cA33+cA23*cA31 
G46 = cA21*cA32+cA22*cA31
G51 = 2*cA11*cA31
G52 = 2*cA12*cA32
G53 = 2*cA13*cA33
G54 = cA12*cA33+cA13*cA32 
G55 = cA11*cA33+cA13*cA31
G56 = cA11*cA32+cA12*cA31
G61 = 2*cA11*cA21
G62 = 2*cA12*cA22
G63 = 2*cA13*cA23
G64 = cA12*cA23+cA13*cA22
G65 = cA11*cA23+cA13*cA21
G66 = cA11*cA22+cA12*cA21

G = as_matrix(((G11, G12, G13, G14, G15, G16), (G21, G22, G23, G24, G25, G26), (G31, G32, G33, G34, G35, G36), (G41, G42, G43, G44, G45, G46), (G51, G52, G53, G54, G55, G56), (G61, G62, G63, G64, G65, G66)))

CMF = G.T*CMM*G

V = VectorFunctionSpace(mesh, "Lagrange", 1)
VV = VectorFunctionSpace(mesh, "Lagrange", 1)
VF = FunctionSpace(mesh, "Lagrange", 1)

c = Expression(("0.0", "0.0", "0.0"))
bcs = DirichletBC(VV, c, end_boundary)

domains = CellFunction("size_t", mesh)
domains.set_all(0)

def gsm(x):
        thi_bound_array = tba(x[2])
        r_bound_array = rba(x[2])
        for ii in range(0,vert_slices):
            if x[2] >= z_bound_array[ii] and x[2] < z_bound_array[ii+1]:
                vh = ii
                break
        for jj in range(0,rad_slices):
                tol = 10**(-10)
                if x[0] > -tol and x[0] < tol:
                    if x[1] > tol:
                        w = math.pi
                    if x[1] >= -tol and x[1] <= tol:
                        w = 0
                    if x[1] < tol:
                        w = -math.pi
                else:
                    w = math.atan(x[1]/x[0])
                    if x[0]<0 and x[1]>=0:
                        w = w+math.pi
                    if x[0]<0 and x[1]<0:
                        w = w+math.pi
                    if x[0]>0 and x[1]<0:
                        w = w+2*math.pi
                if w >= thi_bound_array[jj] and w < thi_bound_array[jj+1]:
                    th = jj
                    break
                elif w < thi_bound_array[0]:
                    th = 0
        for kk in range(0,num_of_radial_devisions):
                    r = math.sqrt(x[0]**2+x[1]**2)
                    if r >= r_bound_array[kk] and r < r_bound_array[kk+1]:
                        rh = kk
                        break
        dom = vh*num_of_radial_devisions*rad_slices + th*num_of_radial_devisions + rh
        return(dom)

#class get_GS(Expression):
#    def eval(self, values, x):
#            values[0] = 0#np.random.normal(stress_r, stress_sd_r, 1)
#            values[1] = 0#np.random.normal(stress_t, stress_sd_t, 1)
#            values[2] = 10#norm_l[cd]#np.random.normal(stress_l, stress_sd_l, 1)
#    def value_shape(self):
#        return (3,)
#get_gs = get_GS()     


for c in cells(mesh):
    cla = Cell_get_vertex_coordinates(c)
    cl = np.zeros(3)
    cl[0] = (cla[0] + cla[3] + cla[6] + cla[9])/4
    cl[1] = (cla[1] + cla[4] + cla[7] + cla[10])/4
    cl[2] = (cla[2] + cla[5] + cla[8] + cla[11])/4
    dom = gsm(cl)
    domains[c] = dom

dx = Measure("dx")[domains]

#for ii in range(0,num_of_sd):
#    sm = SubMesh(mesh, domains, ii)
#    plot(sm, interactive = True)

du = TrialFunction(VV)            # Incremental displacement
v  = TestFunction(VV)             # Test function
u  = Function(VV)                 # Displacement from previous iteration
E = 0.5*(grad(u)+(grad(u)).T) 

stress1s = CMF[0,0]*E[0,0] + CMF[0,1]*E[1,1] + CMF[0,2]*(E[2,2]) + CMF[0,3]*E[1,2] + CMF[0,4]*E[0,2] + CMF[0,5]*E[0,1]
stress2s = CMF[1,0]*E[0,0] + CMF[1,1]*E[1,1] + CMF[1,2]*(E[2,2]) + CMF[1,3]*E[1,2] + CMF[1,4]*E[0,2] + CMF[1,5]*E[0,1]
stress3s = CMF[2,0]*E[0,0] + CMF[2,1]*E[1,1] + CMF[2,2]*(E[2,2]) + CMF[2,3]*E[1,2] + CMF[2,4]*E[0,2] + CMF[2,5]*E[0,1]
stress4s = CMF[3,0]*E[0,0] + CMF[3,1]*E[1,1] + CMF[3,2]*(E[2,2]) + CMF[3,3]*E[1,2] + CMF[3,4]*E[0,2] + CMF[3,5]*E[0,1]
stress5s = CMF[4,0]*E[0,0] + CMF[4,1]*E[1,1] + CMF[4,2]*(E[2,2]) + CMF[4,3]*E[1,2] + CMF[4,4]*E[0,2] + CMF[4,5]*E[0,1]
stress6s = CMF[5,0]*E[0,0] + CMF[5,1]*E[1,1] + CMF[5,2]*(E[2,2]) + CMF[5,3]*E[1,2] + CMF[5,4]*E[0,2] + CMF[5,5]*E[0,1] 

#psi = 0.5*((stress1s+0)*E[0,0]+(stress2s+0)*E[1,1]+(stress3s+1000)*(E[2,2])+stress4s*E[1,2]+stress5s*E[0,2]+stress6s*E[0,1])
# Total potential energy
Pi = 0
for ii in range(0,num_of_sd):
    Pi = Pi + (0.5*((stress1s+norm_r[ii])*E[0,0]+(stress2s+norm_t[ii])*E[1,1]+(stress3s+norm_l[ii])*(E[2,2])+stress4s*E[1,2]+stress5s*E[0,2]+stress6s*E[0,1]))*dx(ii)

# Compute first variation of Pi (directional derivative about u in the direction of v)
F = derivative(Pi, u, v)

# Compute Jacobian of F
J = derivative(F, u, du)

# Solve variational problem
solve(F == 0, u, bcs, J=J,
      form_compiler_parameters=ffc_options)


#plot(u, mode = "displacement", interactive = True)

u_vec = project(u, V).vector().array()
coor_int = interpolate(Expression(("x[0]", "x[1]", "x[2]")), V).vector().array()
new_coor = coor_int+u_vec
num_of_verts = mesh.num_vertices()
num_of_cells = mesh.num_cells()
dofs_to_vert = np.zeros(num_of_verts, dtype=np.uintp)
vectordofs_to_vert = np.zeros((num_of_verts*3), dtype=np.uintp)
vectordofs_to_subvert = np.zeros((num_of_verts*3),dtype=np.uintp)
cellinds = np.zeros((num_of_cells,4), dtype=np.uintp)
cell_inds = np.zeros((num_of_cells,1), dtype=np.uintp)
dm = VF.dofmap()
    
dms = [V.sub(i).dofmap() for i in range(3)]
for cell in cells(mesh):
        cell_ind = cell.index()
        cell_inds[cell_ind] = cell.index()
        vert_inds = cell.entities(0)
        cellinds[cell_ind,:] = vert_inds #retuns the local indacies for each cell      
        dofs_to_vert[dm.cell_dofs(cell_ind)] = vert_inds
        for i, (dms_i, dmcs_i) in enumerate(zip(dms, dms)):
            vectordofs_to_vert[dms_i.cell_dofs(cell_ind)] = vert_inds #gives map for coords to cells, not sencitive to xyz position
            vectordofs_to_subvert[dms_i.cell_dofs(cell_ind)] = i #gives map for each x,y,z to to particular cells for map above
map_mat = np.zeros((len(new_coor),3), dtype=float)
coor_cur = np.zeros((len(new_coor)/3,3), dtype=float)
map_mat[:,0] = vectordofs_to_vert
map_mat[:,1] = vectordofs_to_subvert
map_mat[:,2] = new_coor

for ij in range(0,len(map_mat)):
        if map_mat[ij,1] == 0:
            coor_cur[map_mat[ij,0],0] = map_mat[ij,2]
        if map_mat[ij,1] == 1:
            coor_cur[map_mat[ij,0],1] = map_mat[ij,2]
        if map_mat[ij,1] == 2:
            coor_cur[map_mat[ij,0],2] = map_mat[ij,2] 


mesh_coor = mesh.coordinates()

mind = np.where(mesh_coor[:,2]>399)
mind = mind[0]

top_points = np.zeros((len(mind),2))
cur_points = np.zeros((len(mind),2))
for ii in range(0,len(top_points)):
    top_points[ii,0] = mesh_coor[mind[ii],0]*math.cos(theta) + mesh_coor[mind[ii],1]*(-math.sin(theta))
    top_points[ii,1] = mesh_coor[mind[ii],0]*math.sin(theta) + mesh_coor[mind[ii],1]*math.cos(theta)
    cur_points[ii,0] = coor_cur[mind[ii],0]*math.cos(theta) + coor_cur[mind[ii],1]*(-math.sin(theta))
    cur_points[ii,1] = coor_cur[mind[ii],0]*math.sin(theta) + coor_cur[mind[ii],1]*math.cos(theta)

mct = mesh_coor
mcc = coor_cur
for ii in range(0,len(top_points)):
    mct[mind[ii],0] = top_points[ii,0]
    mct[mind[ii],1] = top_points[ii,1]
    mcc[mind[ii],0] = cur_points[ii,0]
    mcc[mind[ii],1] = cur_points[ii,1]


mindp = np.where((mct[:,2]>399) & (mct[:,0]<cutp) & (mct[:,0]>cutp-1))
#mindp = np.where((mesh_coor[:,2]>399) & (mesh_coor[:,0]<cutp) & (mesh_coor[:,0]>cutp-1))

mindp = mindp[0]
mindn = np.where((mct[:,2]>399) & (mct[:,0]>cutn) & (mct[:,0]<cutn+1))

#mindn = np.where((mesh_coor[:,2]>399) & (mesh_coor[:,0]>cutn) & (mesh_coor[:,0]<cutn+1))
mindn = mindn[0]
if len(mindn) > len(mindp):
    print("diferent num of elements on each side of cut")
    ef = 0
    while len(mindn) > len(mindp):
        print(len(mindn), len(mindp))
        if ef==0:
            mindn = np.delete(mindn, 0)
            ef = 1
        elif ef==1:
            mindn = np.delete(mindn, len(mindn)-1)
            ef = 0
        else:
            print("ef val != 0, 1")
        
elif len(mindn) < len(mindp):
    print("diferent num of elements on each side of cut")
    ef = 0
    while  len(mindp) > len(mindn):
        print(len(mindp), len(mindn))
        if ef==0:
            mindp = np.delete(mindp, 0)
            ef = 1
        elif ef==1:
            mindp = np.delete(mindp, len(mindp)-1)
            ef = 0
        else:
            print("ef val != 0, 1")        


#if(abs(theta) <= math.pi/2):
disp_vec = mcc[mindp, 0] - mcc[mindn,0]
#elif(abs(theta) > math.pi/2):
#    disp_vec = coor_cur[mindn,0] - coor_cur[mindp, 0]
#else:
#    print("theta value poorly defined")
measure = np.mean(disp_vec)

print(disp_vec)
print(stress_l)
print(np.mean(disp_vec))

with open('tfile_var.csv', 'wb') as fh:
    fh.write(str(measure))
fh.close()

#note opening for average sample should be ~4.2

#subdomain tut http://fenicsproject.org/documentation/tutorial/materials.html
