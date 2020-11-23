import numpy as np

from enum import Enum
from geopy.distance import geodesic
from scipy.spatial.distance import squareform

_AVG_EARTH_RADIUS_KM = 6371.0088

class Unit(Enum):
    """
    Enumeration of supported units.
    The full list can be checked by iterating over the class; e.g.
    the expression `tuple(Unit)`.
    """

    KILOMETERS = 'km'
    METERS = 'm'
    MILES = 'mi'
    NAUTICAL_MILES = 'nmi'
    FEET = 'ft'
    INCHES = 'in'


# Unit values taken from http://www.unitconversion.org/unit_converter/length.html
_CONVERSIONS = {Unit.KILOMETERS:       1.0,
                Unit.METERS:           1000.0,
                Unit.MILES:            0.621371192,
                Unit.NAUTICAL_MILES:   0.539956803,
                Unit.FEET:             3280.839895013,
                Unit.INCHES:           39370.078740158}


def get_avg_earth_radius(unit):
    unit = Unit(unit)
    return _AVG_EARTH_RADIUS_KM * _CONVERSIONS[unit]


def geodesicDistance(X):
    ''' 
    Calculate the pair-wise differences between lat/lon points on an ellipsoid 
    
    Parameters
    ----------
    X - a N x (2,3) matrix that represents either (lat, lon) coordinates or 
         (lat, lon, heights) coordinates on the WGS84 ellipsoid

    Returns
    -------
    Y - a condensed distance matrix Y. For each i and j (where i < j < m), where
         m is the number of original obesrvations. the metric dist(X[i], X[j] is 
         computed and stored in entry ij
    '''
    [N, D] = X.shape
    Y = np.zeros((N, N))

    for i1 in range(N):
        for i2 in range(N):
            if (i1 < i2) and (i2 < N):
                Y[i1,i2] = geodesic(X[i1,:2], X[i2,:2]).km

    return Y[~(Y==0)]
    

def geodesic_distance(x,y):
    return geodesic(x, y).km

def haversine_vector(array1, array2, unit=Unit.KILOMETERS):
    '''
    The exact same function as "haversine", except that this
    version replaces math functions with numpy functions.
    This may make it slightly slower for computing the haversine
    distance between two points, but is much faster for computing
    the distance between two vectors of points due to vectorization.
    '''
    # unpack latitude/longitude
    lat1, lng1 = array1[:, 0], array1[:, 1]
    lat2, lng2 = array2[:, 0], array2[:, 1]

    # convert all latitudes/longitudes from decimal degrees to radians
    lat1 = np.radians(lat1)
    lng1 = np.radians(lng1)
    lat2 = np.radians(lat2)
    lng2 = np.radians(lng2)

    # calculate haversine
    lat = lat2 - lat1
    lng = lng2 - lng1
    d = np.square(np.sin(lat/2)) + np.cos(lat1) * np.cos(lat2) * np.square(np.sin(lng/2))

    return 2 * get_avg_earth_radius(unit) * np.arcsin(np.sqrt(d))


def explode(X):
    m, n = X.shape
    if n != 2:
        raise NotImplementedError
    Y1 = np.empty((int(m*(m-1)/2),n))
    Y2 = np.empty((int(m*(m-1)/2),n))
    ind = 0
    for k in range(m):
        for j in range(k+1,m):
            Y1[ind,:] = X[k,:]
            Y2[ind,:] = X[j,:]
            ind = ind+1
    return Y1, Y2


def distanceJ(X):
    P1, P2 = explode(X)
    dists = haversine_vector(P1, P2)
    return dists


