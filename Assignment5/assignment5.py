import math


def euler_method(f, h, base_case, num_of_iterations):
    '''

    :param num_of_iterations:
    :param f:
    :param h:
    :param base_case: tuple such that y(base_case[0]) = base_case[1]
    :return:
    '''
    x = base_case[0]  # 0
    y = base_case[1]  # 1
    evaluation = 0
    for _ in range(num_of_iterations):
        evaluation += 1
        y = y + h * f(x, y)
        x += h
    return x, y, evaluation


def runge_kutta_2(f, lamda, h, base_case, num_of_iterations):
    '''

    :param f:
    :param lamda:
    :param h:
    :param base_case:
    :param num_of_iterations:
    :return:
    '''
    x = base_case[0]
    y = base_case[1]
    two_times_lamda = 2 * lamda
    half_time_lamda = (1 / two_times_lamda)
    alpha2 = half_time_lamda
    alpha1 = 1 - alpha2
    h_times_lamda = h * lamda
    evaluation = 0
    for _ in range(num_of_iterations):
        evaluation += 2
        k1 = f(x, y)
        k2 = f(x + h_times_lamda, y + h_times_lamda * k1)
        y = y + h * (alpha1 * k1 + alpha2 * k2)
        x += h

    return x, y, evaluation


def exact_function(f, x):
    return f(x)


def find_minimal_for_epsilon(f, epsilon):
    i = 0
    curr = float('inf')
    while curr > epsilon:
        i += 1
        curr = f(i)
    return i


if __name__ == "__main__":
    def f(x, y):
        return x + y


    def exact(x):
        e = math.e
        return x, 2 * (e ** x) - x - 1


    euler_results = []
    rk_results = []
    exact_results = []
    euler_better_results = []
    h = 0.1
    base_case = (0, 1)
    lamda = 2 / 3
    print("Eulers Method h=0.1:")
    for i in range(11):
        x, y, evalutaion = euler_method(f, h, base_case, i)
        euler_results.append(y)
        print("x={:.1f},\t\ty={:.8f},\t\tnum of evaluations={}".format(x, y, evalutaion))

    print("Eulers Method h=0.05:")
    for i in range(21):
        x, y, evalutaion = euler_method(f, 0.05, base_case, i)
        if i % 2 == 0:
            euler_better_results.append(y)
            print("x={:.1f},\t\ty={:.8f},\t\tnum of evaluations={}".format(x, y, evalutaion))

    print("RK2 method:")
    for i in range(11):
        x, y, evalutaion = runge_kutta_2(f, lamda, h, base_case, i)
        rk_results.append(y)
        print("x={:.1f},\t\ty={:.8f},\t\tnum of evaluations={}".format(x, y, evalutaion))

    print("Exact function 2*e**x - x - 1")
    for i in range(11):
        x, y = exact_function(exact, h * i)
        exact_results.append(y)
        print("x={:.1f},\t\ty={:.8f}".format(x, y))

    print("\t\teuler error\t\t\trk2 error")
    for i in range(len(exact_results)):
        euler = euler_results[i]
        rk = rk_results[i]
        exact_result = exact_results[i]
        print("x={:.1f}\t{:.8f}\t\t\t{:.8f}".format(0.1 * i, exact_result - euler, exact_result - rk))

    print("\t\tRK2 - euler")
    for i in range(len(exact_results)):
        euler = euler_results[i]
        rk = rk_results[i]
        print("x={:.1f}\t{:.8f}".format(0.1 * i, rk - euler))

    print("\t\teuler h=0.05 error\trk2 error")
    for i in range(len(exact_results)):
        euler = euler_better_results[i]
        rk = rk_results[i]
        exact_result = exact_results[i]
        print("x={:.1f}\t{:.8f}\t\t\t{:.8f}".format(0.1 * i, exact_result - euler, exact_result - rk))
