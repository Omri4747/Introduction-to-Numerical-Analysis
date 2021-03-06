import math

# [i] python has a simple typing mechanism.
# ': <type>' signifies the type of the parameter (integer in the following)
def most_significant(num: int, digits_to_keep: int):
    """ return the digits_to_keep most significant digits of num
    :return (the kept digits, the number of non-significant (zeroed- out) digits)
    """
    # [i] 'assert' is a special keyword in python.
    # it verifies if the condition is True
    assert digits_to_keep > 0, 'digits_to_keep should be positive'
    num_digits = math.floor(math.log10(num)) + 1
    non_significant = max(0, num_digits - digits_to_keep)
    # [i] '//' is division without reminder in python 3+
    return (num // 10 ** non_significant), non_significant


def adder(man_a: int, exp_a: int, man_b: int, exp_b: int):
    if exp_a > exp_b:
        return adder(man_b, exp_b, man_a, exp_a)
    exp_diff = exp_b - exp_a
    man_b *= (10 ** exp_diff)
    man_output, exp_output = most_significant(man_a + man_b, 3)
    exp_output += exp_a
    return man_output, exp_output


def accumulate(n: int):
    man_acc = 100
    exp_acc = -5
    man_c = 400
    exp_c = -5
    for i in range(n):
        man_acc, exp_acc = adder(man_a=man_acc, exp_a=exp_acc, man_b=man_c, exp_b=exp_c)

    real_out = 0.001 + 0.004 * n
    approximate_out = man_acc * (10 ** exp_acc)
    absolute_err = math.fabs(real_out - approximate_out)
    return absolute_err


def main():
    err_70 = accumulate(70)
    err_72 = accumulate(72)
    err_7000 = accumulate(7000)
    err_8000 = accumulate(8000)
    err_8002 = accumulate(8002)
    print(f"err 70 = {err_70}")
    print(f"err 7000 = {err_7000}")
    print(f"err72 - err70 = {err_72 - err_70}")
    print(f"err8002 - err8000 = {err_8002 - err_8000}")


if __name__ == "__main__":
    main()