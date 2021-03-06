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

def p35_search(sec_lvl):
    bitlen = sec_lvl*3
    print("p35 primes near %d bit quantum security level in the problem of Finding Path in Isogeny Graph" % (sec_lvl))
    l3 = int(bitlen*math.log(2, 3))
    l5 = int(bitlen*math.log(2, 5))
    for k in range(l3-10, l3+10):
        for m in range(l5-10, l5+10):
            for p2 in [4]:
                p = p2 * (3**k) * (5**m) - 1
                if miller_rabine_test(p):
                    print("p = %d * 3^%d * 5^%d - 1 = 0x%x; ||p|| = %d bit; q_comp(E[3^%d]) = %d qbit; q_comp(E[5^%d]) = %d qbit." % ( p2, k, m, p, math.ceil(math.log2(p)), k, math.ceil(k*math.log2(3)/3), m, math.ceil(m*math.log2(5)/3) ) )

def main():
    p35_search(10)

if __name__ == "__main__":
    main()