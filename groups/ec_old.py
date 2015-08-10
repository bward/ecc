from cmath import sqrt
from random import randint
import number_theory as nt
from groups.ecpoint import ECPoint


def weierstrass(curve):
    a_1, a_2, a_3, a_4, a_6 = curve
    b_2 = a_1**2 + 4*a_2
    b_4 = a_1*a_3 + 2*a_4
    b_6 = a_3**2 + 4*a_6
    c_4 = b_2**2 - 24*b_4
    c_6 = -b_2**3 + 36*b_2*b_4 - 216*b_6
    return (-27*c_4, -54*c_6)


class EllipticCurve():
    def __init__(self, curve, p, e=1):
        self.p = p
        self.e = e
        if len(curve) > 2:
            self.curve = weierstrass(curve)
        else:
            self.curve = curve

    def point(self, value):
        return ECPoint(value, self.curve, self.p, self.e)

    def order(self, base_field=False):
        a, b = self.curve[0] % self.p, self.curve[1] % self.p
        if base_field:
            count = 1
            for x in range(self.p):
                y_squared = (x**3+a*x+b) % self.p
                if y_squared == 0:
                    count += 1
                elif nt.legendre_symbol(y_squared, self.p) == 1:
                    count += 2
            return count
        else:
            base_order = self.order(True)
            t = self.p + 1 - base_order
            a = t + sqrt(t**2 - 4*self.p)
            b = t - sqrt(t**2 - 4*self.p)
            return round((1 + self.p**self.e - a**self.e - b**self.e).real)

    def bad_reduction_primes(self):
        d = 4*self.curve[0]**3+27*self.curve[1]**2
        return set(nt.prime_factorise(abs(d)))

    def reduction_type(self):
        if self.p in self.bad_reduction_primes():
            if self.curve[0] % self.p == 0:
                return 'bad reduction, additive'
            else:
                return 'bad reduction, multiplicative'
        else:
            return 'good reduction'

    def weil_pairing(self, P, Q, m):
        def g(R, S):
            x_1, y_1 = R.value
            x_2, y_2 = S.value
            if R == -S:
                def h(T):
                    return (T.value[0] - x_1) % self.p
            else:
                if R == S:
                    m = (3*x_1**2 + self.curve[0]) * \
                        nt.mod_mult_inv(2*y_1, self.p)
                else:
                    m = (y_2 - y_1)*nt.mod_mult_inv(x_2 - x_1, self.p)

                def h(T):
                    num = T.value[1] - y_1 - m*(T.value[0] - x_1)
                    den = T.value[0] + x_1 + x_2 - m**2
                    return (num * nt.mod_mult_inv(den, self.p)) % self.p
            return h

        binary_m = bin(m)
        n = len(binary_m) - 2
        ms = [int(binary_m[i]) for i in range(n+1, 1, -1)]

        def f(R, x):
            T = R
            ops = []
            for i in range(n-2, -1, -1):
                ops.append((True, 2))
                ops.append((False, T, T))
                T = T+T
                if ms[i] == 1:
                    ops.append((False, T, R))
                    T = T+R
            out = 1
            for op in ops:
                if op[0]:
                    out = out**2 % self.p
                else:
                    out = (out * g(op[1], op[2])(x)) % self.p
            return out

        a, b = self.curve
        x = randint(1, self.p)

        while True:
            x = randint(1, self.p)
            if x == P.value[0] or x == -Q.value[0]:
                continue
            elif x == (P-Q).value[0]:
                continue
            elif nt.legendre_symbol(x**3+a*x+b, self.p) != 1:
                continue
            break

        y = nt.modular_sqrt(x**3+a*x+b, self.p)
        S = self.point((x, y))
        out = f(P, Q+S)*nt.mod_mult_inv(f(P, S), self.p) * \
            nt.mod_mult_inv(f(Q, P-S)*nt.mod_mult_inv(f(Q, -S), self.p), self.p)
        return out % self.p






