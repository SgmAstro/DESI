"""
Rotates the (ra,dec) coordinates of the BGS 1% mocks
"""
import numpy as np
from astropy.modeling.rotations import EulerAngleRotation


def radec2xyz(ra, dec):
    """
    Convert ra and dec (in degrees) to 2d array of x, y, z unit vectors

    Args:
        ra: array of ra (in degrees)
        dec: array of dec (in degrees)
    Returns:
        2D array of unit vectors with shape (3,N)
    """
    theta = (90-dec) * np.pi / 180.
    phi = ra * np.pi / 180.

    x = np.sin(theta)*np.cos(phi)
    y = np.sin(theta)*np.sin(phi) 
    z = np.cos(theta)

    return np.vstack((x, y, z))


def xyz2radec(vector):
    """
    Convert 2d array of x, y, z unit vectors to ra and dec (in degrees)

    Args:
        vector: 2D array of unit vectors with shape (3,N)
    Returns:
        array of ra (in degrees)
        array of dec (in degrees)
    """
    x, y, z = np.vsplit(vector, 3)

    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z/r)
    phi = np.arctan2(y,x)

    ra = phi * 180 / np.pi
    dec = 90 - (theta* 180 / np.pi)

    ra[ra<0] += 360
    ra[ra>360] -= 360

    return ra, dec


def rotate(ra, dec, rotation_matrix):
    """
    Rotates ra and dec coordinates on the sky, using a rotation matrix
    
    Args:
        ra:  array of ra (in degrees)
        dec: array of dec (in degrees)
        rotation_matrix: 3x3 array containing the rotation matrix
    Returns:
        array of ra after rotation (in degrees)
        array of dec after rotation (in degrees)
    """

    vector = radec2xyz(ra, dec)
    vector_rot = np.dot(rotation_matrix, vector)
    ra_new, dec_new = xyz2radec(vector_rot)
    
    return ra_new[0], dec_new[0]


def get_rotation_matrix(theta, ra, dec):
    """
    Returns the rotation matrix needed to rotate by angle theta about a point
    on the sky
    
    Args:
        theta: rotation angle (in degrees)  
        ra:    ra coordinate that the rotation is performed around (in degrees)
        dec:   dec coordinate that the rotation is performed around (in degrees)
    Returns:
        3x3 array containing the rotation matrix
    """
    # convert theta to radians
    theta *= np.pi/180.
    
    # convert ra, dec to a 3D Cartesian vector
    x,y,z = radec2xyz(ra, dec)
    rotation_axis=[x[0],y[0],z[0]]
    
    # get rotation matrix
    R = np.zeros((3,3))
    for j in range(3):
        for k in range(3):
            if j==k:
                R[j][k] = np.cos(theta/2.)**2 + \
                    np.sin(theta/2.)**2 * (2*rotation_axis[j]**2 - 1)
            else:
                l = 3-(j+k)
                R[j][k] = 2*rotation_axis[j]*rotation_axis[k]*\
                    np.sin(theta/2.)**2 + \
                    (-1)**((l-k+3)%3)*rotation_axis[l]*np.sin(theta)
                    
    return R


def multiply_rotation_matrices(*Rs):
    """
    Performs matrix multiplication to combine a series of rotation matrices into
    a single matrix
    
    Args:
        *Rs: rotation matrices
    Returns:
        3x3 array containing the combined rotation matrix
    """
    R = Rs[0]
    for i in range(1,len(Rs)):
        R = np.matmul(R, Rs[i])
    return R


def __rotate_orientation_v1_1(ra, dec, N, flip=False, inverse=False):
    """
    Takes the ra, dec coordinates of one of the 1% mocks, read directly from the file,
    and rotates the coordinates so that they cover the 1% footprint.
    This function is for the first orientation of version 1 of the footprint, 
    which is near dec=0.
    
    Args:
        ra:  array of ra (in degrees)
        dec: array of dec (in degrees)
        N:   integer from 0 to 11, which specifies which replication this is, along ra
        flip:    If True, also flips the coordinates, to cover dec<0. Default is False
        inverse: If True, applies the inverse rotation. Default is False
    Returns:
        array of ra after rotation (in degrees)
        array of dec after rotation (in degrees)
    """
    if flip and not inverse:
        ra, dec = 360-ra, -dec
    
    if inverse:
        R1 = get_rotation_matrix(theta=-14,  ra=20, dec=15)
        R2 = get_rotation_matrix(theta=N*30, ra=0,  dec=90)
        R = multiply_rotation_matrices(R2,R1)
    else:
        R1 = get_rotation_matrix(theta=-N*30, ra=0,  dec=90)
        R2 = get_rotation_matrix(theta=14,    ra=20, dec=15)
        R = multiply_rotation_matrices(R2,R1)
        
    ra_new, dec_new = rotate(ra, dec, R)
        
    if flip and inverse:
        return 360-ra_new, -dec_new
    else:
        return ra_new, dec_new

    
def __rotate_orientation_v1_2(ra, dec, N, flip=False, inverse=False):
    """
    Takes the ra, dec coordinates of one of the 1% mocks, read directly from the file,
    and rotates the coordinates so that they cover the 1% footprint.
    This function is for the second orientation of version 1 of the footprint, 
    which is similar to the first, but with the two regions swapped
    
    Args:
        ra:  array of ra (in degrees)
        dec: array of dec (in degrees)
        N:   integer from 0 to 11, which specifies which replication this is, along ra
        flip:    If True, also flips the coordinates, to cover dec<0. Default is False
        inverse: If True, applies the inverse rotation. Default is False
    Returns:
        array of ra after rotation (in degrees)
        array of dec after rotation (in degrees)
    """
    if flip and not inverse:
        ra, dec = 360-ra, -dec
    
    if inverse:
        R1 = get_rotation_matrix(theta=-14,         ra=20, dec=15)
        R2 = get_rotation_matrix(theta=181,         ra=57, dec=20)
        R3 = get_rotation_matrix(theta=N*30+14.435, ra=0,  dec=90)
        R = multiply_rotation_matrices(R3,R2,R1)
    else:
        R1 = get_rotation_matrix(theta=-N*30-14.435, ra=0,  dec=90)
        R2 = get_rotation_matrix(theta=-181,         ra=57, dec=20)
        R3 = get_rotation_matrix(theta=14,           ra=20, dec=15)
        R = multiply_rotation_matrices(R3,R2,R1)
    
    ra_new, dec_new = rotate(ra, dec, R)
    
    if flip and inverse:
        return 360-ra_new, -dec_new
    else:
        return ra_new, dec_new

    
def __rotate_orientation_v1_3(ra, dec, N, flip=False, inverse=False):
    """
    Takes the ra, dec coordinates of one of the 1% mocks, read directly from the file,
    and rotates the coordinates so that they cover the 1% footprint.
    This function is for the third orientation of version 1 of the footprint, 
    where one of the regions is close to dec=90
    
    Args:
        ra:  array of ra (in degrees)
        dec: array of dec (in degrees)
        N:   integer from 0 to 11, which specifies which replication this is, along ra
        flip:    If True, also flips the coordinates, to cover dec<0. Default is False
        inverse: If True, applies the inverse rotation. Default is False
    Returns:
        array of ra after rotation (in degrees)
        array of dec after rotation (in degrees)
    """
    if flip and not inverse:
        ra, dec = 360-ra, -dec
    
    if inverse:
        R1 = get_rotation_matrix(theta=-14,     ra=20,  dec=15)
        R2 = get_rotation_matrix(theta=-55,     ra=40,  dec=0)
        R3 = get_rotation_matrix(theta=-45,     ra=110, dec=-30)
        R4 = get_rotation_matrix(theta=30*N-60, ra=0,   dec=90)
        R = multiply_rotation_matrices(R4,R3,R2,R1)
    else:
        R1 = get_rotation_matrix(theta=-30*N+60, ra=0,   dec=90)
        R2 = get_rotation_matrix(theta=45,       ra=110, dec=-30)
        R3 = get_rotation_matrix(theta=55,       ra=40,  dec=0)
        R4 = get_rotation_matrix(theta=14,       ra=20,  dec=15)
        R = multiply_rotation_matrices(R4,R3,R2,R1)
        
    ra_new, dec_new = rotate(ra, dec, R) 
        
    if flip and inverse:
        return 360-ra_new, -dec_new
    else:
        return ra_new, dec_new
    
    
def __rotate_orientation_v1_4(ra, dec, N, flip=False, inverse=False):
    """
    Takes the ra, dec coordinates of one of the 1% mocks, read directly from the file,
    and rotates the coordinates so that they cover the 1% footprint.
    This function is for the fourth orientation of version 1 of the footprint, 
    where the two regions are close to dec=+-50
    
    Args:
        ra:  array of ra (in degrees)
        dec: array of dec (in degrees)
        N:   integer from 0 to 11, which specifies which replication this is, along ra
        flip:    If True, also flips the coordinates, to cover dec<0. Default is False
        inverse: If True, applies the inverse rotation. Default is False
    Returns:
        array of ra after rotation (in degrees)
        array of dec after rotation (in degrees)
    """
    if flip and not inverse:
        ra, dec = 360-ra, -dec
    
    if inverse:
        R1 = get_rotation_matrix(theta=-14,     ra=20, dec=15)
        R2 = get_rotation_matrix(theta=-75,     ra=40, dec=0)
        R3 = get_rotation_matrix(theta=30*N-10, ra=0,  dec=90)
        R = multiply_rotation_matrices(R3,R2,R1)
    else:
        R1 = get_rotation_matrix(theta=-30*N+10, ra=0,  dec=90)
        R2 = get_rotation_matrix(theta=75,       ra=40, dec=0)
        R3 = get_rotation_matrix(theta=14,       ra=20, dec=15)
        R = multiply_rotation_matrices(R3,R2,R1)
        
    ra_new, dec_new = rotate(ra, dec, R)

    if flip and inverse:
        return 360-ra_new, -dec_new
    else:
        return ra_new, dec_new

    
def __rotate_orientation_v2(ra, dec, N, flip=False, inverse=False):
    """
    Takes the ra, dec coordinates of one of the 1% mocks, read directly from the file,
    and rotates the coordinates so that they cover the 1% footprint.
    This function is for version 2 of the footprint.
    
    Args:
        ra:  array of ra (in degrees)
        dec: array of dec (in degrees)
        N:   integer from 0 to 11, which specifies which replication this is, along ra
        flip:    If True, also flips the coordinates, to cover dec<0. Default is False
        inverse: If True, applies the inverse rotation. Default is False
    Returns:
        array of ra after rotation (in degrees)
        array of dec after rotation (in degrees)
    """
    
    angles = np.array([0,2,4,11,13,15,22,24,26,33,35,37,44,46,48,55,57,59]) * 5
    
    if flip and not inverse:
        ra, dec = 360-ra, -dec
    
    if inverse:
        R1 = get_rotation_matrix(theta=-30,        ra=217, dec=35)
        R2 = get_rotation_matrix(theta=-angles[N], ra=0,   dec=90)
        R = multiply_rotation_matrices(R2,R1)
    else:
        R1 = get_rotation_matrix(theta=angles[N], ra=0,   dec=90)
        R2 = get_rotation_matrix(theta=30,        ra=217, dec=35)
        R = multiply_rotation_matrices(R2,R1)
        
    ra_new, dec_new = rotate(ra, dec, R)

    if flip and inverse:
        return 360-ra_new, -dec_new
    else:
        return ra_new, dec_new


def rotate_mock(ra, dec, nmock, inverse=False, version=2):
    """
    Rotates the ra, dec coordinates of the 1% mocks on the sky, mapping them
    back to the original 1% survey footprint
    
    Args:
        ra:    array of ra (in degrees)
        dec:   array of dec (in degrees)
        nmock:   integer from 0 to 95, specifying the mock number
        inverse: if True, applies the inverse rotation. Default is False
        version: version of the footprint. Default is version 2
    Returns:
        array of ra after rotation (in degrees)
        array of dec after rotation (in degrees)
    """ 
    
    if version==2:
        
        if nmock < 0 or nmock >= 36:
            raise ValueError("nmock must be between 0 and 35")
    
        if nmock < 18:
            ra_new, dec_new = __rotate_orientation_v2(ra, dec, nmock, inverse=inverse)
        else:
            ra_new, dec_new = __rotate_orientation_v2(ra, dec, nmock-18, flip=True, inverse=inverse)
        
    
    elif version==1:
        
        if nmock < 0 or nmock >= 96:
            raise ValueError("nmock must be between 0 and 95")
            
        if nmock < 24:
            if nmock%2==0:
                ra_new, dec_new = __rotate_orientation_v1_1(ra, dec, nmock//2, inverse=inverse)
            else:
                ra_new, dec_new = __rotate_orientation_v1_2(ra, dec, nmock//2, inverse=inverse)

        elif nmock < 48:
            if nmock%2==0:
                ra_new, dec_new = __rotate_orientation_v1_1(ra, dec, (nmock-24)//2, flip=True, inverse=inverse)
            else:
                ra_new, dec_new = __rotate_orientation_v1_2(ra, dec, (nmock-24)//2, flip=True, inverse=inverse)

        elif nmock < 60:
            ra_new, dec_new = __rotate_orientation_v1_3(ra, dec, (nmock-48), inverse=inverse)
        elif nmock < 72:
            ra_new, dec_new = __rotate_orientation_v1_3(ra, dec, (nmock-60), flip=True, inverse=inverse)

        elif nmock < 84:
            ra_new, dec_new = __rotate_orientation_v1_4(ra, dec, (nmock-72), inverse=inverse)
        elif nmock < 96:
            ra_new, dec_new = __rotate_orientation_v1_4(ra, dec, (nmock-84), flip=True, inverse=inverse)
    
    else: 
        raise ValueError("Invalid version number")
    
    return ra_new, dec_new
    
    
    
def test(version=2):
    """
    Reads in BGS mock mock and makes some test plots
    """
    
    if version==1:
        Nmock = 96
    elif version==2:
        Nmock = 36
    
    import h5py
    import matplotlib.pyplot as plt
        
    # read ra, dec and nmock from BGS mock
    print("Reading mock")
    f = h5py.File("/project/projectdirs/desi/mocks/bgs/MXXL/full_sky/v0.0.4/BGS_r20.6.hdf5","r")
    ra = f["Data/ra"][...][::200]
    dec = f["Data/dec"][...][::200]
    f.close()
    
    f = h5py.File("/project/projectdirs/desi/mocks/bgs/MXXL/one_percent/one_percent_v%i.hdf5"%version,"r")
    nmock = f["nmock"][...][::200]
    f.close()
    
    # plot the ra, dec coordinates of the 96 1% mocks
    print("Plotting 1% mocks")
    plt.figure(dpi=200)
    for i in range(Nmock):
        idx = nmock==i
        ra_i, dec_i = ra[idx], dec[idx]
        plt.scatter(ra_i, dec_i, s=1, edgecolor="none")
    plt.xlim(0,360)
    plt.ylim(-90,90)
    plt.xlabel("ra (deg)")
    plt.ylabel("dec (deg)")
    plt.show()
    
    # plot the ra, dec coordinates of the 96 1% mocks after rotation
    print("Rotating mocks")
    plt.figure(dpi=200)
    for i in range(Nmock):
        idx = nmock==i
        ra_i, dec_i = ra[idx], dec[idx]
        ra_i, dec_i = rotate_mock(ra_i, dec_i, nmock=i, version=version)
        plt.scatter(ra_i, dec_i, s=1, edgecolor="none")
    plt.xlim(0,360)
    plt.ylim(-90,90)
    plt.xlabel("ra (deg)")
    plt.ylabel("dec (deg)")
    plt.show()
    
    # plot the ra, dec coordinates of the 96 1% mocks after inverse rotation
    print("Inverse rotation")
    plt.figure(dpi=200)
    for i in range(Nmock):
        idx = nmock==i
        ra_i, dec_i = ra[idx], dec[idx]
        ra_i, dec_i = rotate_mock(ra_i, dec_i, nmock=i, version=version)
        ra_i, dec_i = rotate_mock(ra_i, dec_i, nmock=i, inverse=True, version=version)
        plt.scatter(ra_i, dec_i, s=1, edgecolor="none")
    plt.xlim(0,360)
    plt.ylim(-90,90)
    plt.xlabel("ra (deg)")
    plt.ylabel("dec (deg)")
    plt.show()
    
    
    
if __name__ == "__main__":
    
    test()