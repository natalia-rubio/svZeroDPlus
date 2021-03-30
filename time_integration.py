# coding=utf-8
import numpy as np
import scipy
from scipy.sparse import csr_matrix
import pdb
import copy
import matplotlib.pyplot as plt


def min_ydot_least_sq_init(neq, eps_min, yinit, block_list, args, dt, rho, eps_factor=5.0):
    # System : min (over y) ||Fy+C||^2 + eps||Ay-yinit||^2
    # Inversion equation:
    # y <-- inv(F'F+eps*D) (-F'C+eps*D*yinit)
    # yinit : Desired set of initial conditions
    # D : diag([(i in yinit) for i in range(neq) ] )
    # If yinit == zeros, set D = I

    # ydot is solved in an approximate sense: ydot <-- np.linalg.lstsq(E,-Fy-C)

    # Solve as a sequence of problems : eps ---> eps_min
    eps = 10.0
    iit = 0
    args['Time'] = 0
    # y0 = np.zeros(neq)

    y0 = yinit

    E, F, C, dE, dF, dC = initialize_solution_matrices(neq)

    if np.linalg.norm(yinit) == 0.:
        D = np.eye(neq)
    else:
        D = np.diag([_ != 0 for _ in yinit])

    print("Approximate consistent initialization : \n\n")

    while eps > eps_min:
        iit += 1
        args['Solution'] = y0
        assemble_structures(E, F, C, dE, dF, dC, args, block_list)

        M = np.dot(F.transpose(), F) + eps * D
        v = -np.dot(F.transpose(), C) + eps * np.dot(D, yinit)
        y0, _, _, _ = np.linalg.lstsq(M, v)

        ydot0, _, _, _ = np.linalg.lstsq(E, -np.dot(F, y0) - C)

        print("Iteration ", iit, ", Initializing residual: ", np.linalg.norm(form_rhs_NR(E, F, C, y0, ydot0)))
        eps = eps / eps_factor

    return y0, ydot0


def min_ydot_cons_least_sq_init(neq, eps_min, yinit, block_list, args, dt, rho, eps_factor=5.0):
    # System : min (over y) ||Fy+C||^2 + sum_j(lambda_j(y_j-yinit_j))
    # Inversion equation:
    # [2(F'F)  A' ][  y0  ] = [ -2F'C ]
    # [  A     0  ][lambda] = [ yinit ]

    # ydot is solved in an approximate sense: ydot <-- np.linalg.lstsq(E,-Fy-C)

    # Solve as a sequence of problems : eps ---> eps_min
    eps = 1.0
    iit = 0
    args['Time'] = 0

    y0 = yinit
    ydot0 = np.zeros(neq)

    indx = np.nonzero(yinit)

    E, F, C, dE, dF, dC = initialize_solution_matrices(neq)

    print("Approximate consistent initialization : \n\n")

    isize = indx[0].size

    if isize == 0:
        while eps > eps_min:
            iit += 1
            args['Solution'] = y0
            assemble_structures(E, F, C, dE, dF, dC, args, block_list)
            M = np.dot(F.transpose(), F)
            v = -np.dot(F.transpose(), C + np.dot(E, ydot0))
            y0, _, _, _ = np.linalg.lstsq(M, v)
            ydot0, _, _, _ = np.linalg.lstsq(E, -np.dot(F, y0) - C)
            print("Iteration ", iit, ", Initializing residual: ", np.linalg.norm(form_rhs_NR(E, F, C, y0, ydot0)))
            print("Iteration ", iit, ", Initializing time derivative size: ", np.linalg.norm(ydot0))
            eps = eps / eps_factor

    else:
        A = np.zeros((isize, neq))
        i = 0
        print(indx[0])
        for j in indx[0]:
            A[i, j] = 1
            i += 1

        while eps > eps_min:
            iit += 1
            args['Solution'] = y0
            assemble_structures(E, F, C, dE, dF, dC, args, block_list)

            T = 2 * np.dot(F.transpose(), F)

            M = np.block([
                [T, A.transpose()],
                [A, np.zeros((isize, isize))]
            ])

            v = np.block([-2 * np.dot(F.transpose(), C + np.dot(E, ydot0)), yinit[indx]]).transpose()

            yM, _, _, _ = np.linalg.lstsq(M, v)

            y0 = yM[:neq]
            ydot0, _, _, _ = np.linalg.lstsq(E, -np.dot(F, y0) - C)

            print("Iteration ", iit, ", Initializing residual: ", np.linalg.norm(form_rhs_NR(E, F, C, y0, ydot0)))
            print("Iteration ", iit, ", Initializing time derivative size: ", np.linalg.norm(ydot0))
            eps = eps / eps_factor

    return y0, ydot0


class GenAlpha:
    """
    Solves system E*ydot + F*y + C = 0 with generalized alpha and Newton-Raphson for non-linear residual
    """
    def __init__(self, rho, y):
        # Constants for generalized alpha
        self.alpha_m = 0.5 * (3.0 - rho) / (1.0 + rho)
        self.alpha_f = 1.0 / (1.0 + rho)
        self.gamma = 0.5 + self.alpha_m - self.alpha_f

        # problem dimension
        self.n = y.shape[0]

        # stores matrices E, F, vector C, and tangent matrices dE, dF, dC
        self.mat = {}

        # jacobian matrix
        self.M = []

        # residual vector
        self.res = []

        # initialize matrices in self.mat
        self.initialize_solution_matrices()

    def initialize_solution_matrices(self):
        """
        Create empty dense matrices and vectors
        """
        mats = ['E', 'F', 'dE', 'dF', 'dC']
        vecs = ['C']

        for m in mats:
            self.mat[m] = np.zeros((self.n, self.n))
        for v in vecs:
            self.mat[v] = np.zeros(self.n)

    def assemble_structures(self, block_list):
        """
        Assemble block matrices into global matrices
        """
        for bl in block_list:
            for n in self.mat.keys():
                # vectors
                if len(self.mat[n].shape) == 1:
                    for i in range(len(bl.mat[n])):
                        self.mat[n][bl.global_row_id[i]] = bl.mat[n][i]
                # matrices
                else:
                    for i in range(len(bl.mat[n])):
                        for j in range(len(bl.mat[n][i])):
                            self.mat[n][bl.global_row_id[i], bl.global_col_id[j]] = bl.mat[n][i][j]

    def form_matrix_NR(self, dt):
        """
        Create Jacobian matrix
        """
        self.M = (self.mat['F'] + (self.mat['dE'] + self.mat['dF'] + self.mat['dC'] + self.mat['E'] * self.alpha_m / (
                    self.alpha_f * self.gamma * dt)))

    def form_rhs_NR(self, y, ydot):
        """
        Create residual vector
        """
        self.res = - np.dot(self.mat['E'], ydot) - np.dot(self.mat['F'], y) - self.mat['C']
        # return - csr_matrix(E).dot(ydot) - csr_matrix(F).dot(y) - C

    def form_matrix_NR_numerical(self, res_i, ydotam, args, block_list, epsilon):
        """
        Numerically compute the Jacobian by computing the partial derivatives of the residual using forward finite differences
        """
        # save original values for restoration later
        yaf_original = copy.deepcopy(args['Solution']) # yaf_i

        # compute numerical Jacobian
        J_numerical = np.zeros((self.n, self.n))
        for jj in range(self.n):

            yaf_step_size = np.zeros(self.n)
            yaf_step_size[jj] = np.abs(yaf_original[jj])  * epsilon

            # get solution at the i+1 step
            args['Solution'] = yaf_original  + yaf_step_size # yaf_ip1

            for b in block_list:
                b.update_solution(args)
            self.initialize_solution_matrices()
            self.assemble_structures(block_list)
            self.form_rhs_NR(args['Solution'], ydotam)

            # use forward finite differences (multiply by -1 b/c form_rhs_NR creates the negative residual)
            J_numerical[:, jj] = (self.res - res_i) / yaf_step_size[jj] * -1

        # restore original quantities
        args['Solution'] = yaf_original

        for b in block_list:
            b.update_solution(args)
        self.initialize_solution_matrices()
        self.assemble_structures(block_list)
        self.form_rhs_NR(args['Solution'], ydotam)

        return J_numerical

    def check_jacobian(self, res_i, ydotam, args, block_list):
        """
        Check if the analytical Jacobian (computed from form_matrix_NR) matches the numerical Jacobian
        """

        epsilon_list = np.power(10, np.linspace(-6, 4, 25))

        fig, axs = plt.subplots(self.n, self.n, figsize = (20, 20))
        for epsilon in epsilon_list:
            J_numerical = self.form_matrix_NR_numerical(res_i, ydotam, args, block_list, epsilon)
            error = np.abs(self.M - J_numerical)
            for ii in range(self.n):
                for jj in range(self.n):
                    axs[ii, jj].loglog(epsilon, error[ii, jj], 'k*-')

        for ax in axs.flat:
            ax.set(xlabel='epsilon', ylabel='error')

        fig.suptitle('absolute error vs epsilon')
        plt.show()

    def step(self, y, ydot, t, block_list, args, dt, nit=30):
        """
        Perform one time step
        """
        # initial guess for time step
        curr_y = y.copy() + 0.5 * dt * ydot
        curr_ydot = ydot.copy() * ((self.gamma - 0.5) / self.gamma)

        # Substep level quantities
        yaf = y + self.alpha_f * (curr_y - y)
        ydotam = ydot + self.alpha_m * (curr_ydot - ydot)

        # initialize solution
        args['Time'] = t + self.alpha_f * dt
        args['Solution'] = yaf

        # initialize blocks
        for b in block_list:
            b.update_constant()
            b.update_time(args)

        self.res = [1e16]
        iit = 0
        while np.max(np.abs(self.res)) > 5e-4 and iit < nit:
            # update solution-dependent blocks
            for b in block_list:
                b.update_solution(args)

            # update residual and jacobian
            self.assemble_structures(block_list)
            self.form_rhs_NR(yaf, ydotam)
            self.form_matrix_NR(dt)

            # perform finite-difference check of jacobian if requested
            if args['check_jacobian']:
                if args['Time'] > dt:
                    self.check_jacobian(copy.deepcopy(self.res), ydotam, args, block_list)

            # solve for Newton increment
            dy = scipy.sparse.linalg.spsolve(csr_matrix(self.M), self.res)

            # update solution
            yaf += dy
            ydotam += self.alpha_m * dy / (self.alpha_f * self.gamma * dt)

            if np.any(np.isnan(self.res)):
                raise RuntimeError('Solution nan')

            args['Solution'] = yaf
            iit += 1

        if iit >= nit:
            print("Max NR iterations reached at time: ", t, " , max error: ", max(abs(self.res)))

        # update time step
        curr_y = y + (yaf - y) / self.alpha_f
        curr_ydot = ydot + (ydotam - ydot) / self.alpha_m

        args['Time'] = t + dt

        return curr_y, curr_ydot
