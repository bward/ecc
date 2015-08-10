import math
from collections import Counter
from random import randint


def euclid(a, b):
    if abs(b) > abs(a):
        a, b = b, a
    p, q = divmod(a, b)
    if q == 0:
        return b
    else:
        return euclid(b, q)


def extended_euclid(a, b, coefficients=None, switched=False):
    if coefficients is None:
        coefficients = []
    if abs(b) > abs(a):
        a, b, switched = b, a, True
    p, q = divmod(a, b)
    coefficients.append(p)

    if type(a) is int:
        one = 1
    else:
        one = a.mult_identity()

    if q == 0:
        if len(coefficients) > 2:
            m, n = one, -coefficients[-2]
            for c in coefficients[-3::-1]:
                m, n = n, m - n * c
        elif len(coefficients) == 2:
            m, n = one, -coefficients[0]
        else:
            m, n = one, one - coefficients[0]
        return (n, m, b) if switched else (m, n, b)
    else:
        return extended_euclid(b, q, coefficients, switched)


def solve_congruence(a, b, m):
    a = a % m
    b = b % m
    p, q, hcf = extended_euclid(a, m)
    if b % hcf:
        return False
    elif hcf > 1:
        return solve_congruence(a // hcf, b // hcf, m // hcf)
    else:
        return (p * b) % m, m


def mod_mult_inv(x, m):
    x = x % m
    if x == 0:
        return False

    if solve_congruence(x, 1, m):
        return solve_congruence(x, 1, m)[0]
    else:
        return False


def prime_factorise(n, start=2, counted=False):
    if counted:
        return Counter(prime_factorise(n, start=start))
    for i in range(start, int(math.sqrt(n))+1):
        if n % i == 0:
            return [i] + prime_factorise(n//i, start=i)
    return [n]


def chinese_remainder_theorem(congruences):
    # congruences a list of the form [(c,n)] (x=c(n))
    cs, ns = zip(*congruences)
    N = 1
    for n in ns:
        N *= n
    x = 0
    for (c_i, n_i) in congruences:
        s_i = extended_euclid(n_i, N//n_i)[1]
        e_i = s_i * (N//n_i)
        x += c_i * e_i

    return x % N, N


def is_b_smooth(x, b):
    factors = prime_factorise(x)
    if max(factors) <= b:
        return True
    else:
        return False


def eratosthenes(limit):
    is_prime = [False] * 2 + [True] * (limit - 1)
    for n in range(int(limit**0.5 + 1.5)):
        if is_prime[n]:
            for i in range(n*n, limit+1, n):
                is_prime[i] = False
    return [i for i, prime in enumerate(is_prime) if prime]


def primitive_root(p):
    factors = set(prime_factorise(p-1))
    for i in range(1, p):
        for f in factors:
            if pow(i, (p-1)//f, p) == 1:
                break
        else:
            return i


def mod(M, n):
    return list(map(lambda l: list(map(lambda x: x % n, l)), M))


def reduced_row_echelon(M, p):
    M = mod(M, p)
    l = len(M[0]) - 2

    # first go to echelon form
    for i in range(l+1):
        if not mod_mult_inv(int(M[i][i]), p):
            # we have some linear dependence. search for a new row and swap in
            for k in range(i, len(M)):
                if mod_mult_inv(int(M[k][i]), p):
                    # swap the rows
                    row = M[i]
                    M[i] = M[k]
                    M[k] = row
                    break
            else:
                raise ValueError('Not enough LI relations. Try again?')

        inv = mod_mult_inv(int(M[i][i]), p)
        M[i] = [(M[i][j]*inv) % p for j in range(len(M[0]))]

        for j in range(i+1, len(M)):
            M[j] = [(M[j][k] - M[j][i]*M[i][k]) % p for k in range(len(M[0]))]

    # check if last row is 0
    if not mod_mult_inv(int(M[l][l]), p):
        for k in range(l+1, len(M)):
            if mod_mult_inv(int(M[k][l]), p):
                # swap the rows
                M[[k, l]] = M[[l, k]]
                break
        else:
            raise ValueError('Not enough LI relations. Try again?')

    # M is now in row echelon form: next, reduce by back-substitution
    for i in range(l, 0, -1):
        for j in range(i-1, -1, -1):
            M[j] = [(M[j][k] - M[j][i]*M[i][k]) % p for k in range(len(M[0]))]

    return M


def legendre_symbol(a, p):
    try:
        a = a % p
    except ZeroDivisionError:
        return 0
    factorisation = prime_factorise(a, counted=True)
    out = 1
    for factor in factorisation.keys():
        if factor == p-1:
            if p % 4 == 3:
                out *= (-1)**factorisation[factor]
        elif factor == 2:
            if p % 8 == 3 or p % 8 == 5:
                out *= (-1)**factorisation[factor]
        else:
            if factor % 4 == 3 and p % 4 == 3:
                out *= (-legendre_symbol(p, factor))**factorisation[factor]
            else:
                out *= (legendre_symbol(p, factor))**factorisation[factor]

    return out


def cmod(z, p, i):
    return (round(z.real) % p) + (round(z.imag/math.sqrt(i)) % p)*1j


def square(c, d, n):
    if d is 0:
        return 1
    elif d % 2:
        return (c * square((c**2) % n, (d - 1) // 2, n)) % n
    else:
        return square((c**2) % n, d // 2, n)


def square_p2(x, p, i):
    return (x[0]**2+i*x[1]**2) % p, (2*x[0]*x[1]) % p


def power_p2(x, n, p, i):
    if n is 0:
        return (1, 0)
    elif n % 2:
        new = power_p2(square_p2(x, p, i), (n - 1) // 2, p, i)
        out = ((x[0]*new[0]+i*x[1]*new[1]) % p, (x[1]*new[0] + x[0]*new[1]) % p)
        return out
    else:
        return power_p2(square_p2(x, p, i), n // 2, p, i)


def modular_sqrt(n, p):
    a = 1
    while legendre_symbol(a**2-n, p) != -1:
        a = randint(1, p)
    i = (a**2-n) % p
    x = (a, 1)
    return power_p2(x, (p+1)//2, p, i)[0]
