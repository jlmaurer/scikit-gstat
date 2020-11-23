
def deramp_bilinear(x, y, data):
    ''' Deramp data '''
    A = np.array([x, y, np.ones(len(x))]).T
    ramp = np.linalg.lstsq(A, data.T, rcond=None)[0]
    data = data - (np.matmul(A, ramp))
    return data


