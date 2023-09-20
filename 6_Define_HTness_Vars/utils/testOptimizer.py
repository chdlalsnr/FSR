import numpy as np
from scipy.optimize import minimize



ip = np.array([3, 1])

opti = minimize(function, ip, method='CG')

print(opti)


def function(x):

        return x[0]**2 + x[1]**2 + 1
