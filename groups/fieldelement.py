from groups.groupelement import GroupElement
from copy import deepcopy


class FieldElement(GroupElement):
    def __neg__(self):
        return self.add_inverse()

    def __invert__(self):
        return self.mult_inverse()

    def __sub__(self, other):
        return self.add(other.add_inverse())

    def __truediv__(self, other):
        return self.mult(other.mult_inverse())

    def inverse(self):
        return self.mult_inverse()

    def __mul__(self, other):
        if other is 1:
            return self
        return self.mult(other)

    def __rmul__(self, other):
        out = self.add_identity()
        for i in range(other):
            out += self
        return out

    def __add__(self, other):
        return self.add(other)

    def add_identity(self):
        new = deepcopy(self)
        new.value = new.add_id_value
        return new

    def mult_identity(self):
        new = deepcopy(self)
        new.value = new.mult_id_value
        return new

    def identity(self):
        return self.mult_identity()

    def operate(self, other):
        return self.mult(other)

    def is_zero(self):
        return self.value == self.add_id_value

    def is_one(self):
        return self.value == self.mult_id_value

    def __bool__(self):
        if self.value == self.add_id_value:
            return False
        else:
            return True

    def __eq__(self, other):
        if other == 0:
            return self.add_id_value == self.value
        elif other == 1:
            return self.mult_id_value == self.value
        else:
            return self.value == other.value

    def __mod__(self, other):
        return divmod(self, other)[1]