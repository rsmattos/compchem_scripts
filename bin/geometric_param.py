#! /usr/bin/env python3

##############################################################################
#  Calculates the geometric parameters used to identify the points in the scan
#  can be the distance fo two atoms, the angle between three or the dihedral
#  angle btween four atoms.
#  Input: array ot two, three or four atoms
#  Output: scalar value with the distance, bond or dihedral angles
##############################################################################

import numpy as np

###################      CALCULATE PARAMETER VARIATION     ###################
def distance(p):
    return np.linalg.norm(p[0]-p[1])

def angle(p):
    x = p[0] - p[1]
    y = p[2] - p[1]

    # Making unitary vectors from arbitrary ones
    xu = x/np.linalg.norm(x)
    yu = y/np.linalg.norm(y)

    return np.degrees(np.arccos(np.dot(xu, yu)))

# functino to calculate the dihedral angle from cartesian coordinates, taken from
# https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
def dihedral(p):
    """Praxeolitic formula 1 sqrt, 1 cross product"""
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))