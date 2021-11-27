from __future__ import absolute_import, print_function, division
import matplotlib.pyplot as plt
import math


def bisection_method(f_bisection, a, b, tolerance):
    """
    :param f_bisection: function to determine the find the root of it
    :param a: left limit of the domain
    :param b: right limit of the domain
    :param tolerance: the tolerance which should suffice |z-a'|<tolerance where a' is the root and z is the returned value
    :return: z s.t |z-a'|<tolerance, # of iterations to find the root
    """
    i = 0
    fb = f_bisection(b)
    fa = f_bisection(a)
    if fa * fb > 0:
        raise ValueError("f(a)*f(b) should be negative but is {0}".format(fa * fb))
    elif fa == 0:
        return a, i
    elif fb == 0:
        return b, i

    limit = 2 * tolerance
    while math.fabs(a - b) >= limit:
        i += 1
        z = (a + b) / 2
        fz = f_bisection(z)
        if fz == 0:
            return z, i
        fa = f_bisection(a)
        if fa * fz < 0:
            b = z
        else:
            a = z

    return (a + b) / 2, i


def regula_falsi(f_regula, x0, x1, tolerance, goal):
    """
    :param f_regula: function to determine the find the root of it
    :param x0: initial guess s.t f(x0)*f(x1) < 0
    :param x1: initital guess s.t f(x0)*f(x1) < 0
    :param tolerance: the tolerance to be clos to the root
    :param goal: the goal to return
    :return: x s.t |f(x)| < tolerance
    """
    fx0 = f_regula(x0)
    fx1 = f_regula(x1)
    if fx0 == 0:
        return x0, 0
    elif fx1 == 0:
        return x1, 0
    elif fx0 * fx1 > 0:
        raise ValueError("f(x0)*f(x1) should be negative but is {0}".format(fx0 * fx1))
    i = 0
    xi_minus_1 = x1
    xi_minus_2 = x0
    xi = float('inf')
    while math.fabs(xi - goal) >= tolerance:
        i += 1
        fxi_minus_1 = f_regula(xi_minus_1)
        delta_f = fxi_minus_1 - f_regula(xi_minus_2)
        delta_xi = xi_minus_1 - xi_minus_2
        xi = xi_minus_1 - fxi_minus_1 * delta_xi / delta_f
        fxi = f_regula(xi)
        if fxi == 0:
            return xi, i
        if fxi * fxi_minus_1 < 0:
            xi_minus_2 = xi
        else:
            xi_minus_1 = xi

    return xi, i


def f(x):
    return x ** 2 - 0.2 * x - 3


def main():
    # start a
    x0 = -1
    x1 = 4
    tolerance = 10 ** -8
    bisection_result, bisection_iter_num = bisection_method(f, x0, x1, tolerance)
    print("got result {0} in bisection method with {1} iterations.".format(bisection_result, bisection_iter_num))
    # finish a
    # start b
    tolerance_axis = []
    bisection_iterations_axis = []
    falsi_iterations_axis = []
    for d in range(1, 6):
        d = -1 * d
        tolerance_axis.append(d)
        tolerance = 10 ** d
        bisection, bisection_iter = bisection_method(f, x0, x1, tolerance)
        bisection_iterations_axis.append(bisection_iter)
        falsi_result, falsi_iter = regula_falsi(f, x0, x1, tolerance, bisection_result)
        falsi_iterations_axis.append(falsi_iter)
    plt.plot(tolerance_axis, bisection_iterations_axis, color='orange', marker="o", label='Bisection Method')
    plt.plot(tolerance_axis, falsi_iterations_axis, color='blue', marker="o", label='Regula Falsi')
    plt.xlabel("d")
    plt.ylabel("num of iterations")
    plt.title("Regula Falsi and Bisection methods iterations per 10**d comparison")
    plt.legend()
    plt.show()
    # finish b


if __name__ == '__main__':
    main()
