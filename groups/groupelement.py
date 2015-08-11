from copy import deepcopy


class GroupElement():

    value = None
    id_value = None
    order = None

    def __init__(self, value):
        self.value = value

    def identity(self):
        new = deepcopy(self)
        new.value = new.id_value
        return new

    def __mul__(self, other):
        if other == 1:
            return self
        return self.operate(other)

    def __add__(self, other):
        return self.operate(other)

    def __sub__(self, other):
        return self.operate(other.inverse())

    def __neg__(self):
        return self.inverse()

    def __invert__(self):
        return self.inverse()

    def __pow__(self, n):
        if n < 0:
            return self.inverse()**-n
        elif n == 0:
            return self.identity()
        elif n == 1:
            return self
        elif n % 2:
            return self * (self*self)**((n-1)//2)
        else:
            return (self*self)**(n//2)

    def __rmul__(self, other):
        out = self.identity()
        for i in range(other):
            out += self
        return out

    def __repr__(self):
        return str(self.value)

    def __eq__(self, other):
        if other == 1 or other == 0:
            return self.id_value == self.value
        else:
            return self.value == other.value

    def __lt__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return self.value < other
        else:
            return self.value < other.value

    def __gt__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return self.value > other
        else:
            return self.value > other.value

    def __le__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return self.value <= other
        else:
            return self.value <= other.value

    def __ge__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return self.value >= other
        else:
            return self.value >= other.value

    def __getitem__(self, key):
        return self.value[key]