""" 
    GFp2El - class representing element in the GF(p^2) field, p = 3 (mod 4)
""" 
class GFp2El:
    @staticmethod
    def unity(p):
        return GFp2El(1, 0 , p)
    
    @staticmethod
    def i(p):
        return GFp2El(0, 1 , p)

    def __init__(self, x, y, p):
        self.x = x % p
        self.y = y % p
        self.p = p % p

    def __add__(self, other):
        if type(other) == int:
            other = GFp2El(other, 0, self.p)
        return GFp2El( self.x + other.x % self.p, self.y + other.y % self.p, self.p )

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)

    def __neg__(self):
        return GFp2El( (-self.x) % self.p, (-self.y)% self.p , self.p)

    def __mul__(self, other):
        if type(other) == int:
            return GFp2El( (self.x * other) % self.p, (self.y * other) % self.p, self.p ) 
        return GFp2El( (self.x * other.x - self.y*other.y) % self.p, (self.x * other.y + self.y * other.x) % self.p, self.p )

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        return self * other.inv()

    def __str__(self):
        return hex(self.x) + " + " + hex(self.y) + " i"

    def __eq__(self, other):
        if type(other) == int:
            other = GFp2El(other, 0, self.p)
        return self.x == other.x and self.y == other.y and self.p == other.p

    def pow(self, n):
        a = 1
        t = self.clone()
        for i in range(n.bit_length()):
            if n & (1 << i):
                a = a * t
            t = t*t
        return a

    def clone(self):
        return GFp2El(self.x, self.y, self.p)

    def conj(self):
        return GFp2El(self.x, (-self.y) % self.p, self.p)

    def inv(self):
        s = pow((self.x ** 2 + self.y ** 2) % self.p, self.p - 2, self.p)
        return self.conj() * s
    
    def sqrt(self):
        if not self.is_quad_residue():
            raise RuntimeError("not residue")
        s = self.pow((self.p - 1) // 2)
        if s.x == (-1) % self.p and s.y == 0:
            return GFp2El.i(self.p) * self.pow( (self.p + 1)/4 )
        else:
            return (1 + s).pow( (self.p - 1)//2 ) * self.pow( (self.p + 1)/4 )

    def is_quad_residue(self) -> bool:
        return self.pow( (self.p * self.p - 1) // 2) == 1
    

"""
    SupersingularEllipticCurve - supersingular ellipric curve in Edwards form over GF(p^2)
"""
class SupersingularEllipticCurve:
    # x^2 + y^2 = 1 + d x^2 y^2
    def __init__(self, p, d):
        self.p = p
        self.d = d
        self.n = (p+1)**2
    
    def check_on_curve(self, P):
        Q = P.convert_to_affine()
        if (Q.x * Q.x + Q.y * Q.y) == (1 + self.d * Q.x * Q.x * Q.y * Q.y):
            return True
        return False

    def add(self, P, Q):
        if P.is_affine:
            return EcPoint( 
            (P.x * Q.x - P.y * Q.y) / ( 1 - self.d * P.x * Q.x * P.y * Q.y), 
            (P.x * Q.y + P.y * Q.x) / ( 1 + self.d * P.x * Q.x * P.y * Q.y) )

    def mul(self, P, k):
        pass


class EcPoint:

    def __init__(self, *args):
        if len(args) == 2:
            self.x = args[0]
            self.y = args[1]
            self.is_affine = True
        elif len(args) == 3:
            self.x = args[0]
            self.y = args[1]
            self.z = args[2]
            self.is_affine = False

    def convert_to_affine(self):
        if self.is_affine:
            return self
        else:
            return EcPoint(self.x * self.z.inv(), self.y * self.z.inv())

    def convert_to_proj(self):
        if self.is_affine:
            return EcPoint(self.x, self.y, GFp2El.unity(self.x.p))
        else:
            return self
    

if __name__ == "__main__":
    aa = GFp2El(9, 0, 59) 
    a = GFp2El(9, 10, 59)
    b = GFp2El(2, 39, 59)

    d = a * a.inv()
    print(a)
    print(a.inv())
    print(a * a.inv())
    print(b.is_quad_residue())
