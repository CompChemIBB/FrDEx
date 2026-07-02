#Python script to calculate coordinate of methyl group to add onto Guanine on G-Quadruplex

import math
import sys

#Cartesian coordinates of C8-N9-Me atoms read from Bash script passed as arguments, read from PDB

c8_x=float(sys.argv[1])
c8_y=float(sys.argv[2])
c8_z=float(sys.argv[3])
n9_x=float(sys.argv[4])
n9_y=float(sys.argv[5])
n9_z=float(sys.argv[6])
cme_x=float(sys.argv[7])
cme_y=float(sys.argv[8])
cme_z=float(sys.argv[9])

#Dihedral angle phi and N9-CMe-H bond angle theta passed from bash as argument, calculated on isolted guanine in PCM

phi=float(sys.argv[10])
theta=float(sys.argv[11])

#Definition of the bond vectors: a=C8N9, b=N9CMe, c=CMeH (will be calculated to get H coordinates), d=axb (vector normal to a and b), e=bxc (vector normal to b and c)

a_x=n9_x-c8_x
a_y=n9_y-c8_y
a_z=n9_z-c8_z
b_x=cme_x-n9_x
b_y=cme_y-n9_y
b_z=cme_z-n9_z

#Definition of the modulus of vectors a-e. Mod_c is passed from bash as an argument

mod_a=math.sqrt(a_x**2+a_y**2+a_z**2)
mod_b=math.sqrt(b_x**2+b_y**2+b_z**2)
mod_c=float(sys.argv[12])
mod_d=math.sqrt((a_y*b_z-a_z*b_y)**2+(-a_x*b_z+a_z*b_x)**2+(a_x*b_y-a_y*b_x)**2)
mod_e=math.sqrt(b_x**2+b_y**2+b_z**2)*mod_c*math.sin(theta)

#calculation of the c vector (equation computed in Maple, from the simulataneoius equations of the definition of dihedral angle and bond angle in vector picture). 

c_x=-(-abs(math.cos(theta))*a_x**2*b_x*b_y**2*mod_b*mod_c-abs(math.cos(theta))*a_x**2*b_x*b_z**2*mod_b*mod_c+2*abs(math.cos(theta))*a_x*a_y*b_x**2*b_y*mod_b*mod_c+2*abs(math.cos(theta))*a_x*a_z*b_x**2*b_z*mod_b*mod_c-abs(math.cos(theta))*a_y**2*b_x**3*mod_b*mod_c-abs(math.cos(theta))*a_y**2*b_x*b_z**2*mod_b*mod_c+2*abs(math.cos(theta))*a_y*a_z*b_x*b_y*b_z*mod_b*mod_c-abs(math.cos(theta))*a_z**2*b_x**3*mod_b*mod_c-abs(math.cos(theta))*a_z**2*b_x*b_y**2*mod_b*mod_c+math.cos(phi)*a_x*b_y**2*mod_d*mod_e+math.cos(phi)*a_x*b_z**2*mod_d*mod_e-math.cos(phi)*a_y*b_x*b_y*mod_d*mod_e-math.cos(phi)*a_z*b_x*b_z*mod_d*mod_e-math.sin(phi)*a_y*b_z*mod_b*mod_d*mod_e+math.sin(phi)*a_z*b_y*mod_b*mod_d*mod_e)/(a_x**2*b_x**2*b_y**2+a_x**2*b_x**2*b_z**2+a_x**2*b_y**4+2*a_x**2*b_y**2*b_z**2+a_x**2*b_z**4-2*a_x*a_y*b_x**3*b_y-2*a_x*a_y*b_x*b_y**3-2*a_x*a_y*b_x*b_y*b_z**2-2*a_x*a_z*b_x**3*b_z-2*a_x*a_z*b_x*b_y**2*b_z-2*a_x*a_z*b_x*b_z**3+a_y**2*b_x**4+a_y**2*b_x**2*b_y**2+2*a_y**2*b_x**2*b_z**2+a_y**2*b_y**2*b_z**2+a_y**2*b_z**4-2*a_y*a_z*b_x**2*b_y*b_z-2*a_y*a_z*b_y**3*b_z-2*a_y*a_z*b_y*b_z**3+a_z**2*b_x**4+2*a_z**2*b_x**2*b_y**2+a_z**2*b_x**2*b_z**2+a_z**2*b_y**4+a_z**2*b_y**2*b_z**2)

c_y=(abs(math.cos(theta))*a_x**2*b_y**3*mod_b*mod_c+abs(math.cos(theta))*a_x**2*b_y*b_z**2*mod_b*mod_c-2*abs(math.cos(theta))*a_x*a_y*b_x*b_y**2*mod_b*mod_c-2*abs(math.cos(theta))*a_x*a_z*b_x*b_y*b_z*mod_b*mod_c+abs(math.cos(theta))*a_y**2*b_x**2*b_y*mod_b*mod_c+abs(math.cos(theta))*a_y**2*b_y*b_z**2*mod_b*mod_c-2*abs(math.cos(theta))*a_y*a_z*b_y**2*b_z*mod_b*mod_c+abs(math.cos(theta))*a_z**2*b_x**2*b_y*mod_b*mod_c+abs(math.cos(theta))*a_z**2*b_y**3*mod_b*mod_c+math.cos(phi)*a_x*b_x*b_y*mod_d*mod_e-math.cos(phi)*a_y*b_x**2*mod_d*mod_e-math.cos(phi)*a_y*b_z**2*mod_d*mod_e+math.cos(phi)*a_z*b_y*b_z*mod_d*mod_e-math.sin(phi)*a_x*b_z*mod_b*mod_d*mod_e+math.sin(phi)*a_z*b_x*mod_b*mod_d*mod_e)/((b_x**2+b_y**2+b_z**2)*(a_x**2*b_y**2+a_x**2*b_z**2-2*a_x*a_y*b_x*b_y-2*a_x*a_z*b_x*b_z+a_y**2*b_x**2+a_y**2*b_z**2-2*a_y*a_z*b_y*b_z+a_z**2*b_x**2+a_z**2*b_y**2))

c_z=(abs(math.cos(theta))*a_x**2*b_y**2*b_z*mod_b*mod_c+abs(math.cos(theta))*a_x**2*b_z**3*mod_b*mod_c-2*abs(math.cos(theta))*a_x*a_y*b_x*b_y*b_z*mod_b*mod_c-2*abs(math.cos(theta))*a_x*a_z*b_x*b_z**2*mod_b*mod_c+abs(math.cos(theta))*a_y**2*b_x**2*b_z*mod_b*mod_c+abs(math.cos(theta))*a_y**2*b_z**3*mod_b*mod_c-2*abs(math.cos(theta))*a_y*a_z*b_y*b_z**2*mod_b*mod_c+abs(math.cos(theta))*a_z**2*b_x**2*b_z*mod_b*mod_c+abs(math.cos(theta))*a_z**2*b_y**2*b_z*mod_b*mod_c+math.cos(phi)*a_x*b_x*b_z*mod_d*mod_e+math.cos(phi)*a_y*b_y*b_z*mod_d*mod_e-math.cos(phi)*a_z*b_x**2*mod_d*mod_e-math.cos(phi)*a_z*b_y**2*mod_d*mod_e+math.sin(phi)*a_x*b_y*mod_b*mod_d*mod_e-math.sin(phi)*a_y*b_x*mod_b*mod_d*mod_e)/(a_x**2*b_x**2*b_y**2+a_x**2*b_x**2*b_z**2+a_x**2*b_y**4+2*a_x**2*b_y**2*b_z**2+a_x**2*b_z**4-2*a_x*a_y*b_x**3*b_y-2*a_x*a_y*b_x*b_y**3-2*a_x*a_y*b_x*b_y*b_z**2-2*a_x*a_z*b_x**3*b_z-2*a_x*a_z*b_x*b_y**2*b_z-2*a_x*a_z*b_x*b_z**3+a_y**2*b_x**4+a_y**2*b_x**2*b_y**2+2*a_y**2*b_x**2*b_z**2+a_y**2*b_y**2*b_z**2+a_y**2*b_z**4-2*a_y*a_z*b_x**2*b_y*b_z-2*a_y*a_z*b_y**3*b_z-2*a_y*a_z*b_y*b_z**3+a_z**2*b_x**4+2*a_z**2*b_x**2*b_y**2+a_z**2*b_x**2*b_z**2+a_z**2*b_y**4+a_z**2*b_y**2*b_z**2)

#Calculate H atom coordinates from c vector

hme_x = c_x+cme_x
hme_y = c_y+cme_y
hme_z = c_z+cme_z

#print values, to be used as output in bash script

print hme_x, hme_y, hme_z

