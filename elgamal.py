def encrypt(plain, g, k, A):
    return g**k, A**k * plain


def decrypt(c, a):
    return c[0]**-a * c[1]
