import number_theory as nt
from random import randint


def elgamal_sign(m, key, x):
    p, g = key[0], key[1]
    while True:
        k = randint(2, p-2)
        if nt.euclid(k, p-1) == 1:
            break
    s_1 = pow(g, k, p)
    s_2 = ((m - x*s_1) * nt.mod_mult_inv(k, p-1)) % (p-1)
    if s_2 == 0:
        print('Try again')
        return
    return s_1, s_2


def elgamal_verify(m, key, sig):
    s_1, s_2 = sig[0], sig[1]
    p, g, y = key[0], key[1], key[2]

    if pow(g, m, p) == (pow(y, s_1, p) * pow(s_1, s_2, p)) % p:
        return True
    return False


def dsa_sign(m, key, x):
    q, g = key[0], key[1]
    k = randint(2, q-1)
    s_1 = g**k % q
    s_2 = ((m+x*s_1) * nt.mod_mult_inv(k, q)) % q
    return s_1, s_2


def dsa_verify(m, key, sig):
    s_1, s_2 = sig[0], sig[1]
    q, g, y = key[0], key[1], key[2]
    v_1 = m*nt.mod_mult_inv(s_2, q) % q
    v_2 = s_1 * nt.mod_mult_inv(s_2, q) % q
    if (g**v_1*y**v_2) % q == s_1:
        return True
    return False