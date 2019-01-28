import math
import random as rd
K = 10

def miller_rabine_test(n):

    def extended_euclid(a, b):
        if (b == 0):
            return a, 1, 0
        d, x, y = extended_euclid(b, a % b)
        return d, y, x - a // b * y


    def pre_division_test(n):
        if n % 2 == 0 or n % 3 == 0 or n % 5 == 0 or n % 7 == 0 or n % 11 == 0 or n % 13 == 0:
            return False
        else:
            return True

    d = n - 1
    s = 0
    k = 0
    while d % 2 == 0:
        d = d // 2
        s += 1
    while k <= K:
        x = rd.randint(2, n - 1)
        if extended_euclid(x, n)[0] > 1:
            return False
        x_d = pow(x, d, n)
        if x_d == 1 or x_d == n - 1:
            k += 1
        else:
            fl = False
            for r in range(1, s):
                x_r = pow(x, d * 2**r, n)
                if x_r == n - 1:
                    fl = True
                    k += 1
                    break
                elif x_r == 1:
                    return False
                else:
                    continue
            if not fl:
                return False
    return True

def p35_search(bitlen):
    l3 = int(bitlen*math.log(2, 3))
    l5 = int(bitlen*math.log(2, 5))
    print("l3, l5: ", l3, l5)
    for k in range(l3-20, l3+20):
        for m in range(l5-20, l5+20):
            for p2 in [4]:
                p = p2 * (3**k) * (5**m) - 1
                if miller_rabine_test(p):
                    print("p found: ", p2, k, m, p)

def main():
    p35_search(384)


if __name__ == "__main__":
    main()