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
    def __init__(self, curve, field):
        self.field = field
        if len(curve) > 2:
            self.curve = weierstrass(curve)
            self.curve = list(field.elt([i]) for i in curve)
        else:
            self.curve = list(field.elt([i]) for i in curve)
        if self.field.p in self.bad_reduction_primes():
            print('Warning: over the field ' + str(field) +
                  ', this field has ' + self.reduction_type())

    def point(self, x, y):
        try:
            return ECPoint((self.field.elt(x), self.field.elt(y)), self)
        except:
            x = x.value
            y = y.value
            return ECPoint((self.field.elt(x), self.field.elt(y)), self)

    def order(self, base_field=False):
        a, b = self.curve[0].value[0].value, self.curve[1].value[0].value
        p = self.field.p
        e = self.field.e
        if base_field:
            count = 1
            for x in range(p):
                y_squared = (x**3+a*x+b) % p
                if y_squared == 0:
                    count += 1
                elif nt.legendre_symbol(y_squared, p) == 1:
                    count += 2
            return count
        else:
            base_order = self.order(True)
            t = p + 1 - base_order
            a = (t + sqrt(t**2 - 4*p))/2
            b = (t - sqrt(t**2 - 4*p))/2
            return round((1 + p**e - a**e - b**e).real)

    def bad_reduction_primes(self):
        a, b = self.curve[0].value[0].value, self.curve[1].value[0].value
        d = 4*a**3+27*b**2
        if d == 0:
            return [self.field.p]
        return set(nt.prime_factorise(abs(d)))

    def reduction_type(self):
        if self.field.p in self.bad_reduction_primes():
            if self.curve[0] == 0:
                return 'bad (additive) reduction'
            else:
                return 'bad (multiplication) reduction'
        else:
            return 'good reduction'

    def weil_pairing(self, P, Q):
        m = P.order()
        assert(m == Q.order())

        def g(R, S):
            x_1, y_1 = R.value
            x_2, y_2 = S.value
            if R == -S:
                def h(T):
                    return T.value[0] - x_1
            else:
                if R == S:
                    m = (3*x_1**2 + self.curve[0]) * (2*y_1).inverse()
                else:
                    m = (y_2 - y_1)*(x_2 - x_1).inverse()

                def h(T):
                    num = T.value[1] - y_1 - m*(T.value[0] - x_1)
                    den = T.value[0] + x_1 + x_2 - m**2

                    if not den.inverse() and num == 0:
                        return self.field.elt(1)
                    elif not den.inverse():
                        return 'inf'
                    return num * den.inverse()
            return h

        binary_m = bin(m)
        n = len(binary_m) - 2
        ms = [int(binary_m[i]) for i in range(n+1, 1, -1)]

        def f(R, x):
            if R == x:
                return self.field.elt(0)
            T = R
            ops = []
            for i in range(n-2, -1, -1):
                ops.append((True, 2))
                ops.append((False, T, T))
                T = T+T
                if ms[i] == 1:
                    ops.append((False, T, R))
                    T = T+R
            out = self.field.elt(1)

            for op in ops:
                if op[0]:
                    out = out**2
                else:
                    # zeroes and infinity cancel
                    if type(g(op[1], op[2])(x)) is str or g(op[1], op[2])(x) == 0:
                        continue
                    out = (out * g(op[1], op[2])(x))
            return out

        a, b = self.curve[0].value[0].value, self.curve[1].value[0].value
        x = randint(1, self.field.p)

        while True:
            x = randint(1, self.field.p)
            f_x = self.field.elt(x)
            if f_x == P.value[0] or f_x == -Q.value[0]:
                continue
            elif P-Q != 0 and f_x == (P-Q).value[0]:
                continue
            elif nt.legendre_symbol(x**3+a*x+b, self.field.p) != 1:
                continue
            y = nt.modular_sqrt(x**3+a*x+b, self.field.p)
            S = self.point(x, y)

            break

        out = f(P, Q+S)*f(P, S).inverse() * (f(Q, P-S)*f(Q, -S).inverse()).inverse()
        return out

    def __repr__(self):
        out = 'y^2=x^3'
        if self.curve[0]:
            out += '+' + str(self.curve[0]) + 'x'
        if self.curve[1]:
            out += '+' + str(self.curve[1])

        return 'Elliptic curve ' + out + ' over ' + str(self.field)







