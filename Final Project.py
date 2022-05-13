#!/usr/bin/env python
# coding: utf-8

# In[53]:


import numpy as np
import os
import math
file_location = os.path.abspath('H2O.xyz') #This is the test file location
molecule_file = np.genfromtxt(fname=file_location, skip_header=2, dtype='unicode')
elements = molecule_file[:,0] #these are the elements included in the file
coordinates = molecule_file[:,1:] #these are the coordinates of each element
coordinates_float = coordinates.astype(float)
middle_atom=coordinates_float [0] #since this is a three atom molecule, this is the middle molecule
atom_2=coordinates_float [1] #This is one of the non-middle molecules
atom_3=coordinates_float [2] #This is the other
def distance_finder(one,two) : #This will find the distance between the middle molecule and the non-middle molecules
    [x1,y1,z1] = one
    [x2,y2,z2] = two
    return np.sqrt(((x2-x1)**2)+((y2-y1)**2)+((z2-z1)**2)) #this is the formula for finding distance
middle_atom_coords = np.array(middle_atom) #these next three lines form the distance array
atom_2_coords = np.array(atom_2)
atom_3_coords = np.array(atom_3)
distance_array = [middle_atom_coords,atom_2_coords,atom_3_coords]
final_distance=[] #this is for the final distance without repeats and just the molecules that actually connect
x=0
for value in distance_array:
    distance=distance_finder(value,atom_2_coords)
    if 1.5>distance>0:
        final_distance.append((elements[x], 'to H1:',  F'{distance:.3f}'))
        distance_string= str(distance)
        hypotenus= distance
    distance=distance_finder(value,atom_3_coords)
    if 1.5>distance>0:
        final_distance.append((elements[x], 'to H2:', F'{distance:.3f}'))
#above is calculating the distance from the center atom to the atoms that are connected to it
    distance=distance_finder(atom_2_coords,atom_3_coords)
    if distance>1.5:
        final_distance.append((elements[x], F'{distance:.3f}'))
        distance_string=str(distance)
        opposite=distance/2
        bond_angle= math.asin(opposite/hypotenus)
        bond_angle_final=(bond_angle*180/math.pi)*2
    #Above is using trigonometry to find the bond angle between the three atoms
print("The bond angle for the molecule is", F'{bond_angle_final:.3f}', "degrees")
#The bond angle is then printed
x=x+1


# In[ ]:
