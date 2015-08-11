from groups.fieldelement import FieldElement
import number_theory as nt


class IntModP(FieldElement):
    id_value = 1
    add_id_value = 0

    def __init__(self, value, p):
        self.value = value % p
        self.p = p

    def mult(self, other):
        if self.p != other.p:
            raise ArithmeticError("Elements in different groups")

        new_value = (self.value * other.value) % self.p
        return IntModP(new_value, self.p)

    def add(self, other):
        if self.p != other.p:
            raise ArithmeticError("Elements in different groups")

        new_value = (self.value + other.value) % self.p
        return IntModP(new_value, self.p)

    def mult_inverse(self):
        return IntModP(nt.mod_mult_inv(self.value, self.p), self.p)

    def add_inverse(self):
        return IntModP((-self.value) % self.p, self.p)

    def order(self):
        return self.p - 1

    def __eq__(self, other):
        if other == 0:
            return self.add_id_value == self.value
        elif other == 1:
            return self.id_value == self.value
        else:
            return self.value == other.value and self.p == other.p

    def __mod__(self, q):
        return self.value % q
