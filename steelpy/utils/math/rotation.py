'''
Quaternion  and tensor rotations.
'''

import numpy as np

def rodrot(theta, rotaxis=[0, 0, 1], style='row'):
    """
    Establishes 3D rotation matrix based on Euler-Rodrigues formula.
    See https://en.wikipedia.org/wiki/Euler-Rodrigues_formula.


    Arguments
    -------------
        theta : float
            the rotation angle (in radians)
        rotaxis : [0, 0, 1]
            vector defining rotation axis
        style : 'row'
            unit vectors as stacked as 'row' or 'column'

    Returns
    -----------
        T : float
            transformation matrix in NumPy format

    """
    axis = np.asarray(rotaxis)
    axis = rotaxis/np.sqrt(np.dot(rotaxis, rotaxis))    # Normalize
    a = np.cos(theta/2.0)
    b, c, d = axis*np.sin(theta/2.0)
    a2, b2, c2, d2 = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    T = np.array([[a2+b2-c2-d2, 2*(bc-ad), 2*(bd+ac)],
                  [2*(bc+ad), a2+c2-b2-d2, 2*(cd-ab)],
                  [2*(bd-ac), 2*(cd+ab), a2+d2-b2-c2]])
    
    if style=='row':
        T = T.T
    
    return T

#
def rot_from_R(R):
    '''
    Calculate rotation vector from (deformation part of) rotation tensor.

    Arguments
    -------------
    R : float
        numpy 3x3 matrix describing to rotation transformation
        matrix equivalent to the incremental rotation vector

    Returns
    -----------
    rot : float
        3-by-1 numpy array with three rotations of a node    

    '''

    rot = 0.5*np.array([R[2,1] - R[1,2],
                        R[0,2] - R[2,0],
                        R[1,0] - R[0,1]])
    
    theta = np.linalg.norm(rot)

    if theta<1e-12:
        scale = 1.0
    else:
        scale = theta/np.sin(theta)
    
    return scale * rot

#
def R_from_drot(drot):
    '''
    Calculate the rotation tensor from an incremental rotation vector.

    Arguments
    -----------
    drot : float
        3-by-1 numpy array with three rotation increments of a node

    Returns
    -------------
    R : float
        numpy 3x3 matrix describing to rotation transformation
        matrix equivalent to the incremental rotation vector

    '''

    theta = np.linalg.norm(drot)
    if theta<1e-12:
        return np.eye(3)
    else:
        n = drot/theta
        N = np.zeros([3,3])
        N[0, 1] =- n[2]
        N[0, 2] = n[1]
        N[1, 2] =- n[0]  
        I = np.eye(3)
        #
        N = N - N.T
        R = I + N * np.sin(theta) + (np.outer(n,n) - I) * (1 - np.cos(theta))
        #
        return R

#
def quat_from_R(R):
    '''
    Calculate quaternions from matrix [R]

    Arguments
    ----------
    R : float
        numpy 3x3 matrix describing to rotation transformation
        matrix equivalent to the input quaternion representation
        
    Returns 
    ---------
     r0 : float
         scalar defining the trace of the rotation tensor
     r : float
         numpy array with three quaternions r1,r2,r3

    '''
    e = lambda i,j,k:(i-j)*(j-k)*(k-i)/2

    r0 = (np.trace(R)+1)/4
    r = np.zeros(3)
    for l in range(3):
        e_mat = np.zeros([3,3])
        for i in range(3):
            for j in range(3):
                e_mat[i,j] = e(l, i, j)
        
        r[l] = -np.sum(e_mat*R)/4/r0
    
    return r0, r

#
def R_from_quat(r0, r, row_wise=True): 
    '''
    Calculate transformation matrix [R] according to
    Eq 3.50 in Krenk [1].

    Arguments
    ----------
    r0 : float
        scalar defining the trace of the rotation tensor
    r : float
        numpy array with three quaternions r1,r2,r3
    row_wise : {True, False}
        if converting the result such that basis vectors are
        stacked column-wise (instead of the column-wise used in [1])

    Returns 
    ---------
    R : float
        numpy 3x3 matrix describing to rotation transformation
        matrix equivalent to the input quaternion representation
        
    References
    ------------
    [[1]](../#1) Krenk, 2009.
              
    '''
    r_hat = np.array([[0, -r[2], r[1]], 
                     [r[2], 0, -r[0]], 
                     [-r[1], r[0], 0]]) # Equation 3.7 in Krenk [1]
    R = (r0**2 - np.dot(r,r)) * np.eye(3) + 2*r0*r_hat + 2 * np.outer(r,r)
    
    # Convert such that unit vectors are stacked row-wise
    if row_wise:    
        R = R.T

    return R

#
def increment_from_drot(drot):
    '''
    Establish linearized quaternion increments from 
    given rotation increments, as in Eq. 5.123 in Krenk [1].

    Arguments
    -----------
    drot : float
        3-by-1 or 1-by-3 numpy array with three rotation increments of a node

    Returns
    -----------
    dr0 : float
        scalar defining the trace of the rotation tensor
        of the incremental rotation
    dr : float
        numpy array with three quaternions r1,r2,r3, corresponding
        to the incremental rotation
        
    References
    ------------
    [[1]](../#1) Krenk, 2009.
          
    '''

    dr = 0.5 * drot
    dr0 = np.sqrt(1 - np.dot(dr,dr))

    return dr0, dr

#
def add_increment_from_quat(r0, r, dr0, dr):
    '''
    Add incremental rotation quaternions to initial rotation using Eq. 5.124
    in Krenk [1].

    Arguments
    ----------
    r0 : float
        scalar defining the trace of the initial rotation tensor
    r : float
        numpy array with three quaternions r1,r2,r3 representing
        the inital rotation
    dr0 : float
        scalar defining the trace of the rotation tensor
        of the incremental rotation
    dr : float
        numpy array with three quaternions r1,r2,r3, corresponding
        to the incremental rotation

    Returns 
    ---------
    r0 : float
        scalar defining the trace of the final rotation tensor
    r : float
        numpy array with three quaternions r1,r2,r3 representing
        the final rotation

        
    References
    ------------
    [[1]](../#1) Krenk, 2009.
    '''

    r0 = dr0*r0 - np.dot(dr, r)
    r = dr0*r + r0*dr + np.cross(dr, r)

    return r0, r

#
def add_increment_from_rot(r0, r, drot):
    '''
    Add incremental rotation to initial rotation using Eq. 5.124
    in Krenk [1]. Combining functions inc_from_rot and increment_from_quat.
    
    Arguments
    -----------
    r0 : float
        scalar defining the trace of the initial rotation tensor
    r : float
        numpy array with three quaternions r1,r2,r3 representing
        the inital rotation
    drot : float
        3-by-1 or 1-by-3 numpy array with three rotation increments of a node

    Returns 
    ---------
    r0 : float
        scalar defining the trace of the final rotation tensor
    r : float
        numpy array with three quaternions r1,r2,r3 representing
        the final rotation
        
        
    References
    ------------
    [[1]](../#1) Krenk, 2009.
    '''

    dr0, dr = increment_from_drot(drot)     
    r0, r = add_increment_from_quat(r0, r, dr0, dr)

    return r0,r

#
def quat_mean_and_diff(ra0, ra, rb0, rb):
    '''
    Calculate the mean and difference quaternions from two sets of quaternions
    {ra0, ra} and {rb0, rb}.

    Arguments
    ----------
    ra0 : float
        scalar defining the trace of the rotation tensor
        of the first point
    ra : float
        numpy array with three quaternions r1,r2,r3 representing
        the rotation of the first point
    rb0 : float
        scalar defining the trace of the rotation tensor
        of the second point
    rb : float
        numpy array with three quaternions r1,r2,r3 representing
        the rotation of the second point

    Returns 
    ---------
    r0 : float
        scalar defining the trace of the final rotation tensor
    r : float
        numpy array with three quaternions r1,r2,r3 representing
        the final rotation
    s0 : float
        scalar part of difference quaternion
    s : float
        numpy array with three quaternions s1,s2,s3 representing
        the difference in rotation
        
    References
    ------------
    [[1]](../#1) Krenk, 2009.
    '''    

    s0 = 0.5 * np.sqrt((ra0+rb0)**2 + np.linalg.norm(ra+rb)**2)
    r0 = 0.5 * (ra0+rb0)/s0
    r = 0.5 * (ra+rb)/s0
    s = 0.5 * (ra0*rb - rb0*ra + np.cross(ra,rb))/s0

    return r0, r, s0, s
#
#
#
def transform_unit(e1:list, e2:list|None=None, e3:list|None=None,
                   warnings:bool=False):
    '''
    Establish transformation matrix from e1 and temporary e2 or e3 vectors.

    Arguments
    -----------
    e1 : unit vector describing element longitudinal direction
    e2 : temporary unit vector describing a chosen vector that's perpendicular to the longitudinal direction (approximate y-direction)
    e3 : temporary unit vector describing a chosen vector that's perpendicular to the longitudinal direction (approximate z-direction)
        if both e2 and e3 are different from None, e2 is used (e3 disregarded)

    Returns:
        T : transformation matrix
    '''

    e1 = np.array(e1).flatten()

    if (e2 is not None) and (e3 is not None):
        e3 = None

    if e2 is not None:
        e2 = np.array(e2).flatten()
        e3_tmp = np.cross(e1, e2)  # Direction of the third unit vector
        if np.all(e3 == 0):
            if warnings:
                print('Warning: e1 and e2 identical. Check orientations carefully!')

            e2 = e2 + 0.1
            e3 = np.cross(e1, e2)

        e2 = np.cross(e3_tmp, e1)  # Direction of the second unit vector

        e1 = e1 / np.linalg.norm(e1)  # Normalize the direction vectors to become unit vectors
        e2 = e2 / np.linalg.norm(e2)
        e3 = np.cross(e1, e2)

    elif e3 is not None:
        e3 = np.array(e3).flatten()
        e2_tmp = np.cross(e3, e1)  # Direction of the third unit vector
        e3 = np.cross(e1, e2_tmp)
        e1 = e1 / np.linalg.norm(e1)  # Normalize the direction vectors to become unit vectors
        e3 = e3 / np.linalg.norm(e3)
        e2 = np.cross(e3, e1)

    if e2 is None and e3 is None:
        raise ValueError('Specify either e2 or e3')

    T = np.vstack([e1, e2, e3])

    return T
#
#