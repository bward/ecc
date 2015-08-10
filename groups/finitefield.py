from groups.fieldelement import FieldElement
from groups.intmodp import IntModP
from copy import copy, deepcopy
from random import randint
from number_theory import euclid, extended_euclid


class FiniteField():
    q = None

    def __init__(self, p, e, q=None):
        self.p = p
        self.e = e
        if q:
            self.q = self.elt(q)
        else:
            self.q = self.irreducible_polynomial()
        self.quotient = True

    def elt(self, value):
        return FiniteFieldElement(self, value)

    def irreducible_polynomial(self):
        while True:
            value = [IntModP(randint(1, self.p), self.p) for _ in range(self.e)]
            random_poly = self.elt(value+[IntModP(1, self.p)])
            if random_poly.is_irreducible():
                return random_poly

    def __eq__(self, other):
        return (self.p == other.p and self.q == other.q and self.e == other.e)

    def __repr__(self):
        return 'GF' + str(self.p) + '^' + str(self.e)


class FiniteFieldElement(FieldElement):

    def __init__(self, field, value):
        self.field = field
        self.value = []
        self.add_id_value = [IntModP(0, field.p)]
        self.mult_id_value = [IntModP(1, field.p)]
        if type(value) is tuple or type(value) is list:
            for v in value:
                if type(v) is int:
                    self.value.append(IntModP(v, self.field.p))
                else:
                    self.value.append(v)
            while value[1:] and not self.value[1:][-1]:
                value.pop()
            if self.field.q and self.field.quotient:
                self.value = divmod(self, self.field.q)[1].value
        else:
            self.value = [IntModP(value, self.field.p)]


    def add(self, other):
        l1 = len(self.value)
        l2 = len(other.value)
        if l1 >= l2:
            new_value = copy(self.value)
            for i, v in enumerate(other.value):
                new_value[i] += v
        else:
            return other + self

        while new_value[1:] and not new_value[1:][-1]:
            new_value.pop()
        return self.field.elt(new_value)

    def mult(self, other):
        if self.is_zero() or other.is_zero():
            return self.add_identity()
        degree = len(self.value) + len(other.value) - 1
        new_value = [IntModP(0, self.field.p) for _ in range(degree)]
        for i, a in enumerate(self):
            for j, b in enumerate(other):
                new_value[i+j] += a*b
        while new_value[1:] and not new_value[1:][-1]:
            new_value.pop()
        return self.field.elt(new_value)

    def inverse(self):
        if self.degree() == 0:
            new_value = self.value[0].inverse()
            new = deepcopy(self)
            new.value = [new_value]
            return new

        a, b, c = extended_euclid(self, self.field.q)

        return c.inverse() * a

    def is_irreducible(self):
        x = self.field.elt((0, 1))
        power_term = self.field.elt((0, 1))

        for _ in range(int(self.degree()/2)):
            power_term = power_term ** self.field.p
            gcd = euclid(self, power_term - x)
            if not gcd.is_unit():
                return False
        return True

    def is_unit(self):
        return self.degree() == 0

    def __divmod__(self, divisor):
        self.field.quotient = False
        quotient, remainder = self.add_identity(), self
        if divisor.degree() == 0:
            self.field.quotient = True
            return remainder, quotient

        while remainder.degree() >= divisor.degree():
            exp = remainder.degree() - divisor.degree()
            new_values = [IntModP(0, self.field.p) for _ in range(exp)]
            monomial = self.field.elt(new_values + [remainder[-1] / divisor[-1]])
            quotient += monomial
            remainder -= monomial * divisor
        self.field.quotient = True
        return quotient, remainder

    def add_inverse(self):
        new = deepcopy(self)
        new.value = [-x for x in new.value]
        return new

    def degree(self):
        return len(self.value) - 1

    def __repr__(self):
        out = ''
        first = True
        if self.value[0]:
            out += str(self.value[0])
            first = False
        if self.degree() > 0:
            if not first:
                out += '+'
            if self.value[1] > 1:
                out += str(self.value[1])+'x'
            elif self.value[1] > 0:
                out += 'x'
            for i, v in enumerate(self.value[2:]):
                if not first:
                    out += '+'
                if v > 1:
                        out += str(v)+'x^'+str(i+2)
                elif v > 0:
                    out += 'x^'+str(i+2)
        return out

    def __abs__(self):
        return self.degree()


