# ecc

Example of use:

```python
>>> from groups.finitefield import FiniteField
>>> from groups.ellipticcurve import EllipticCurve
>>> f = FiniteField(5,2) # Finite field of order 5^2
>>> c = EllipticCurve((0, 1), f) # Elliptic curve y^2 = x^3 + 1 over f
>>> p = c.point(2, 2) # Point (2, 2)
>>> q = c.point([0,1], [2,2]) # Point (x, 2+2x)
>>> c.weil_pairing(p, q, 6)
2+4x
>>> c.weil_pairing(p, q, 6)**6
1
```