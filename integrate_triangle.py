import numpy as np
from scipy.special import erf
from scipy.integrate import quad
import sys

def compute_integral(p=(0,0), std_dev=1, t1=(0,0), t2=(1,0), t3=(0,1), num_sub_ints=None):
    p = np.array(p)
    t1 = np.array(t1)
    t2 = np.array(t2)
    t3 = np.array(t3)

    # F maps the traingle T to the triangle Q = { (0,0), (1,0), (0,1) }

    # F(x) = Mx + b
    #
    # where M_inv = (t2-t1|t3-t2)

    M_inv = np.column_stack((t2-t1, t3-t2))

    try:
        M = np.linalg.inv(M_inv)
    except Exception:
        print('M_inv doesnt have an inverse', file=sys.stderr)
        # we shouldn't get to this point but if we do it is because it is a
        # degenerate triangle (Area = 0)
        # so the integral will be 0
        return 0

    b = - np.dot(M, t1)

    # print(M @ t1 + b)
    # print(M @ t2 + b)
    # print(M @ t3 + b)

    CoefOfTransformation = np.abs(np.linalg.det(M_inv))

    # t1 + Minv * <u,v> - p = <Au + Bv + C, Du + Ev + F>

    A=M_inv[0][0]
    B=M_inv[0][1]
    C=(t1-p)[0]
    D=M_inv[1][0]
    E=M_inv[1][1]
    F=(t1-p)[1]

    # print(A, B, C, D, E, F)

    # so we get
    #
    # (Au+Bv+C)^2 + (Du+Ev+F)^2 = T(v + Gu + H)^2 + Ju^2 + Ku + L

    T = B**2 + E**2

    if T == 0:
        G=1
        H=1
    else:
        G = (A*B + D*E) / T
        H = (B*C + E*F) / T

    L=C**2 + F**2 - T * (H ** 2)

    J = A**2 + D**2 - T * (G**2)

    K = 2*(A*C + D*F - T*G*H)

    # print(A, B, C, D, E, F)
    # print(T, G, H, J, K, L)

    def integrand(u):
        exponent=-(J*(u**2)+K*u+L)/(2 * (std_dev**2))
        var=np.sqrt(T)/(np.sqrt(2)*std_dev)
        erf1_arg=(var)*(u + G*u + H)
        erf2_arg=(var)*(G*u + H)

        return np.exp(exponent) * (erf(erf1_arg) - erf(erf2_arg))

    return (CoefOfTransformation / (2 * np.sqrt(2 * np.pi * T) * std_dev)) * quad(integrand, 0, 1)[0]
