from copy import deepcopy
from groups.groupelement import GroupElement


class ECPoint(GroupElement):
    id_value = 0
    test = False

    def __init__(self, value, curve, field):
        a, b = curve
        if value != 0:
            x, y = value
            if y**2 != x**3 + a*x + b:
                raise AttributeError('Not a point on curve', value)
            self.value = (x, y)
        else:
            self.value = 0
        self.curve = (a, b)
        self.field = field

    def __eq__(self, other):
        if other == 1 or other == 0:
            return self.id_value == self.value
        elif self.field != other.field:
            raise ArithmeticError('Curves over different fields')
        elif self.curve != other.curve:
            raise ArithmeticError('Points on different curves')
        else:
            return self.value == other.value

    def operate(self, other):
        if self.field != other.field:
            raise ArithmeticError('Curves over different fields')
        elif self.curve != other.curve:
            raise ArithmeticError('Points on different curves')

        if self.value == 0:
            return other
        elif other.value == 0:
            return self
        else:
            x_1, x_2 = self.value[0], other.value[0]
            y_1, y_2 = self.value[1], other.value[1]
            a = self.curve[0]

            if x_1-x_2 == 0 and y_1+y_2 == 0:
                return ECPoint(0, self.curve, self.field)
            elif self.value == other.value:
                m = (3*x_1**2 + a) * (2*y_1).inverse()
            else:
                m = (y_2-y_1) * (x_2-x_1).inverse()

            x_3 = m**2-x_1-x_2
            y_3 = m*(x_1-x_3) - y_1

            return ECPoint((x_3, y_3), self.curve, self.field)

    def inverse(self):
        if self.value == 0:
            return self
        else:
            new_value = (self.value[0], -self.value[1])
            return ECPoint(new_value, self.curve, self.field)

    def frobenius(self, n):
        new = deepcopy(self)
        for i in range(n):
            new.value = (new.value[0]**self.field.p, new.value[1]**self.field.p)
        return new

    def __mod__(self, q):
        return self.value[0] % q

    # def __pow__(self, n):
    #     if not self.test:
    #         return super(ECPoint, self).__pow__(n)
    #     p = self.p
    #     a = p + 1 - elliptic.order(self.curve, p)
    #     vs = [n % p, (n // p) * a, -(n // p)]

    #     while True:
    #         for i in range(len(vs)):
    #             if abs(vs[i]) >= p:
    #                 v = vs[i]
    #                 vs[i] = v % p
    #                 if i+1 >= len(vs):
    #                     vs.append((v // p) * a)
    #                     vs.append(-(v // p))
    #                 elif i+2 >= len(vs):
    #                     vs[i+1] += (v // p) * a
    #                     vs.append(-(v // p))
    #                 else:
    #                     vs[i+1] += (v // p) * a
    #                     vs[i+2] -= (v // p)

    #         if max(abs(v) for v in vs) < p:
    #             break

    #     out = self.identity()
    #     for (n, v) in enumerate(vs):
    #         if v != 0:
    #             out += super(ECPoint, self.frobenius(n)).__pow__(v)
    #     return out

