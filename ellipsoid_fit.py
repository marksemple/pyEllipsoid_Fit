#! python

import numpy as np


def ellipsoid_fit(point_data, mode=''):
    """ Fit an ellipsoid to a cloud of points using linear least squares
        Based on Yury Petrov MATLAB algorithm: "ellipsoid_fit.m"
    """

    X = point_data[:, 0]
    Y = point_data[:, 1]
    Z = point_data[:, 2]

    # AlGEBRAIC EQUATION FOR ELLIPSOID, from CARTESIAN DATA
    if mode == '':  # 9-DOF MODE
        D = np.array([X * X + Y * Y - 2 * Z * Z,
                      X * X + Z * Z - 2 * Y * Y,
                      2 * X * Y, 2 * X * Z, 2 * Y * Z,
                      2 * X, 2 * Y, 2 * Z,
                      1 + 0 * X]).T

    elif mode == 0:  # 6-DOF MODE (no rotation)
        D = np.array([X * X + Y * Y - 2 * Z * Z,
                      X * X + Z * Z - 2 * Y * Y,
                      2 * X, 2 * Y, 2 * Z,
                      1 + 0 * X]).T

    # THE RIGHT-HAND-SIDE OF THE LLSQ PROBLEM
    d2 = np.array([X * X + Y * Y + Z * Z]).T

    # SOLUTION TO NORMAL SYSTEM OF EQUATIONS
    u = np.linalg.solve(D.T.dot(D), D.T.dot(d2))
    # chi2 = (1 - (D.dot(u)) / d2) ^ 2

    # CONVERT BACK TO ALGEBRAIC FORM
    if mode == '':  # 9-DOF-MODE
        a = np.array([u[0] + 1 * u[1] - 1])
        b = np.array([u[0] - 2 * u[1] - 1])
        c = np.array([u[1] - 2 * u[0] - 1])
        v = np.concatenate([a, b, c, u[2:, :]], axis=0).flatten()

    elif mode == 0:  # 6-DOF-MODE
        a = u[0] + 1 * u[1] - 1
        b = u[0] - 2 * u[1] - 1
        c = u[1] - 2 * u[0] - 1
        zs = np.array([0, 0, 0])
        v = np.hstack((a, b, c, zs, u[2:, :].flatten()))

    else:
        pass

    # PUT IN ALGEBRAIC FORM FOR ELLIPSOID
    A = np.array([[v[0], v[3], v[4], v[6]],
                  [v[3], v[1], v[5], v[7]],
                  [v[4], v[5], v[2], v[8]],
                  [v[6], v[7], v[8], v[9]]])

    # FIND CENTRE OF ELLIPSOID
    centre = np.linalg.solve(-A[0:3, 0:3], v[6:9])

    # FORM THE CORRESPONDING TRANSLATION MATRIX
    T = np.eye(4)
    T[3, 0:3] = centre

    # TRANSLATE TO THE CENTRE, ROTATE
    R = T.dot(A).dot(T.T)

    # SOLVE THE EIGENPROBLEM
    evals, evecs = np.linalg.eig(R[0:3, 0:3] / -R[3, 3])

    # SORT EIGENVECTORS
    # i = np.argsort(evals)
    # evals = evals[i]
    # evecs = evecs[:, i]
    # evals = evals[::-1]
    # evecs = evecs[::-1]

    # CALCULATE SCALE FACTORS AND SIGNS
    radii = np.sqrt(1 / abs(evals))
    sgns = np.sign(evals)
    radii *= sgns

    return (centre, evecs, radii)


if __name__ == "__main__":
    pass

