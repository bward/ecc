from math import floor, ceil, sqrt, exp, log
from random import randint
import groups
import number_theory as nt


def babystep_giantstep(g, h, order=0):  # discrete log problem g^x = h (p)
    if order == 0:
        order = g.p
    n = floor(sqrt(order)) + 1
    powers = [1]

    for i in range(n-1):
        powers.append(g * powers[i])

    y = h
    inv = g**-n

    for j in range(n):
        if y in powers:
            return j*n+powers.index(y)
        y = y * inv


def discrete_log_pp(g, h, q, e):  # solves discrete logs for pp order elts
    coefficients = []

    for i in range(1, e):
        if pow(g, q**i) == 1:
            e = i
            break

    for i in range(e):
        exp = 0

        for j in range(i):
            exp += coefficients[j]*q**j

        g_i = pow(g, q**(e-1))
        h_i = pow(h*pow(g, exp).inverse(), q**(e-i-1))
        coefficients.append(babystep_giantstep(g_i, h_i, order=q))

    x = 0

    for i in range(e):
        x += coefficients[i]*q**i

    return x


def pohlig_hellman(g, h):
    factorisation = nt.prime_factorise(g.order, counted=True)
    ys = []

    for q in factorisation.keys():

        order = q**factorisation[q]
        g_i = pow(g, g.order//order)
        h_i = pow(h, g.order//order)

        ys.append((discrete_log_pp(g_i, h_i, q, factorisation[q]), order))

    return nt.chinese_remainder_theorem(ys)[0]


def index_calculus(g, h, p, b=0):
    if not b:
        b = ceil(exp(sqrt(log(p)*(log(log(p)))))**(1/sqrt(2)))
    primes = nt.eratosthenes(b)
    crt_factors = nt.prime_factorise(p-1, counted=True)
    relations = []

    # we need to find a number of linearly dependent relations
    # first we'll just find 4 times as many as we need
    while len(relations) < 4*len(primes):
        x = randint(1, p-1)
        y = pow(g, x, p)

        if nt.is_b_smooth(y, b):
            factorisation = nt.prime_factorise(y, counted=True)
            relation = [factorisation.get(q, 0) for q in primes] + [x]
            relations.append(relation)

    # perform Gaussian elimination on relations modulo each factor
    sols = []

    for base in sorted(crt_factors.keys()):
        m = base**crt_factors[base]
        sol = nt.reduced_row_echelon(relations, m)
        sols.append([sol[i][-1] for i in range(len(sol[0])-1)])

    small_prime_logs = []

    # use the CRT to stitch together the solutions
    for i in range(len(primes)):
        congruences = []
        for j, base in enumerate(sorted(crt_factors.keys())):
            m = base**crt_factors[base]
            congruences.append((sols[j][i], m))

        small_prime_logs.append(nt.chinese_remainder_theorem(congruences)[0])

    # have the logs of the small primes
    g_inv = nt.mod_mult_inv(g, p)
    g_k = g_inv

    for k in range(1, p):
        v = (h * g_k) % p
        if nt.is_b_smooth(v, b):
            out = k
            factors = nt.prime_factorise(v, counted=True)
            for i, q in enumerate(primes):
                out += small_prime_logs[i] * factors[q]

            return out % (p-1)

        g_k *= g_inv


def _ints_mod_p_shuffle(x, a, b, g, h):
        p = x.p
        if x < p/3:
            return g*x, (a+1) % (p-1), b
        elif x >= p/3 and x < 2*p/3:
            return x**2, (2*a) % (p-1), (2*b) % (p-1)
        else:
            return h*x, a, (b+1) % (p-1)


def pollard_rho(g, h, f):
    x, y = g.identity(), h.identity()
    a, b, c, d = 0, 0, 0, 0
    for i in range(g.order):
        x, a, b = f(x, a, b, g, h)
        y, c, d = f(y, c, d, g, h)
        y, c, d = f(y, c, d, g, h)

        if x == y:
            # found a collision
            u = (a-c) % g.order
            v = (d-b) % g.order
            # know that g^u=h^v
            inv, gcd = nt.extended_euclid(v, g.order)[::2]

            log_mod_gcd = ((u*inv) % g.order)//gcd
            reduced = g.p//gcd

            for j in range(gcd):
                if pow(g, log_mod_gcd+j*reduced) == h:
                    return log_mod_gcd+j*reduced

if __name__ == '__main__':
    g = groups.IntModP(13, 11251)
    h = groups.IntModP(6909, 11251)
    print(pollard_rho(g, h, lambda x, a, b, q, h: _ints_mod_p_shuffle(x, a, b, q, h)))
    print(pohlig_hellman(g, h))
