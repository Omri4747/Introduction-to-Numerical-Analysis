import numpy as np


def power_method(matrix, num_of_iterations=None, epsilon=None):
    if (not num_of_iterations) and (not epsilon):
        raise ValueError("Need to get num of iterations or epsilon difference from last eigenvalue")
    x = np.array([1, 0, 0])
    lamda = float('inf')
    if num_of_iterations:
        for _ in range(num_of_iterations):
            y = matrix @ x
            lamda = np.amax(y)
            y_factored = np.divide(y, lamda)
            x = y_factored
    else:
        while True:
            y = matrix @ x
            new_lamda = np.amax(y)
            if abs(new_lamda-lamda) < epsilon:
                break
            lamda = new_lamda
            y_factored = np.divide(y, lamda)
            x = y_factored
    return lamda, x


if __name__ == "__main__":
    matrix_to_check = np.asarray([[5, -1, 7], [-1, -1, 1], [7, 1, 5]])
    eigen_value, eigen_vector = power_method(matrix=matrix_to_check, num_of_iterations=9)
    print("eigenvalue is " + str(eigen_value))
    print("eigenvector is " + str(eigen_vector))
