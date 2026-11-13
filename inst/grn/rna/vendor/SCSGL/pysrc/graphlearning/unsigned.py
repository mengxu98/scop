
















import numpy as np
import scipy as sp

from ..utils import rowsum_matrix

def _project_hyperplane(v, n):
    """ Project v onto the hyperplane defined by np.sum(v) = -n
    """
    return v - (n + np.sum(v))/(len(v))

def _qp_admm(b, inv_mat_part1, inv_mat_part2, rho=1, max_iter=1000):
    m = len(b)
    n = (1 + np.sqrt(8*m + 1))//2


    w = np.zeros((m, 1))
    lambda_ = np.zeros((m, 1))

    for iter in range(max_iter):

        l_temp = b + rho*w + lambda_
        l = inv_mat_part1@l_temp + inv_mat_part2(l_temp)
        l = _project_hyperplane(l, n)


        w = l - lambda_/rho
        w[w>0] = 0


        lambda_ += rho*(w - l)

        residual = np.linalg.norm(w-l)
        if residual < 1e-4:
            break

    w[w>-1e-4] = 0


    return np.abs(w)

def learn(k, d, alpha):


    n = len(d)
    m = len(k)

    S = rowsum_matrix.rowsum_matrix(n)


    rho = 1


    a = 4*alpha + rho
    b = 2*alpha
    c1 = 1/a
    c2 = b/(a*(a+n*b-2*b))
    c3 = 4*b**2/(a*(a+2*n*b-2*b)*(a+n*b-2*b))
    inv_mat_part1 = c1*sp.sparse.eye(m) - c2*S.T@S
    inv_mat_part2 = lambda x : c3*np.sum(x)*np.ones((m, 1))

    b = S.T@d - 2*k
    if np.ndim(b) == 1:
        b = b[..., None]

    return _qp_admm(b, inv_mat_part1, inv_mat_part2, rho)

def learn_ladmm(k, d, alpha1, alpha2, degree_reg="l2"):



    n = len(d)
    m = len(k)


    if np.ndim(k) == 1:
        k = k[..., None]

    if np.ndim(d) == 1:
        d = d[..., None]


    if degree_reg not in ["lb", "l2"]:
        raise ValueError(
            "The input argument 'degree_reg' must be either 'l2' or 'lb'."
        )

    S = rowsum_matrix.rowsum_matrix(n)

    rho = .1
    mu = 1/(0.9/(rho*(2*(n-1))))

    rng = np.random.default_rng()
    w = _project_hyperplane(rng.uniform(low=0, high=1, size=(n, 1)), -n)
    y = np.zeros((n, 1))
    l = rng.uniform(low=-1, high=0, size=(m, 1))

    primal_res = np.zeros(1000)
    dual_res = np.zeros(1000)
    for iter in range(1000):

        l = np.asarray(S.T@(d - rho*w - y - rho*S@l) - 2*k + mu*l)/(2*alpha2+mu)
        l[l>0] = 0


        w_old = w.copy()
        if degree_reg == "l2":
            w = - np.asarray(rho*S@l + y)/(2*alpha1 + rho)
            w = _project_hyperplane(w, -n)
        elif degree_reg == "lb":
            b = np.asarray(y + rho*S@l)
            w = (- b + np.sqrt(b**2 + 4*rho*alpha1))/(2*rho)


        y += rho*np.asarray(S@l + w)


        primal_res[iter] = np.linalg.norm(rho*S.T@(w-w_old))
        dual_res[iter] = np.linalg.norm(S@l + w)

        if iter > 10:
            if primal_res[iter] < 1e-4 and dual_res[iter] < 1e-4:
                break

    l[l>-1e-4] = 0
    return np.abs(l)
