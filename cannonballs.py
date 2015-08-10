import number_theory as nt


def add(p, q):
    if p == q:
        x, y = (p[0]+97) % 193, p[1]*14 % 193
        s = (3*pow(x, 2, 193)-145) * nt.mod_mult_inv(2*y, 193) % 193
        x_new = (pow(s, 2, 193) - 2*x) % 193
        y_new = (y + s*(x_new - x)) % 193
        return (x_new-97) % 193, -(y_new*69) % 193