from random import randint
from queue import Queue
from pprint import pprint
from graphviz import Digraph
K = 10

def miller_rabine_test(n):

    def extended_euclid(a, b):
        if (b == 0):
            return a, 1, 0
        d, x, y = extended_euclid(b, a % b)
        return d, y, x - a // b * y


    def pre_division_test(n):
        if n % 2 == 0 or n % 3 == 0 or n % 5 == 0 or n % 7 == 0 or n % 11 == 0 or n % 13 == 0:
            return False
        else:
            return True

    d = n - 1
    s = 0
    k = 0
    while d % 2 == 0:
        d = d // 2
        s += 1
    while k <= K:
        x = randint(2, n - 1)
        if extended_euclid(x, n)[0] > 1:
            return False
        x_d = pow(x, d, n)
        if x_d == 1 or x_d == n - 1:
            k += 1
        else:
            fl = False
            for r in range(1, s):
                x_r = pow(x, d * 2**r, n)
                if x_r == n - 1:
                    fl = True
                    k += 1
                    break
                elif x_r == 1:
                    return False
                else:
                    continue
            if not fl:
                return False
    return True

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
        self.p = p

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
        return str(self.x) + " + " + str(self.y) + " i"

    def __repr__(self):
        return self.__str__()

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
    
    def __pow__(self, n):
        return self.pow(n)

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
            return GFp2El.i(self.p) * self.pow( (self.p + 1)//4 )
        else:
            return (1 + s).pow( (self.p - 1)//2 ) * self.pow( (self.p + 1)//4 )

    def is_quad_residue(self) -> bool:
        return self.pow( (self.p * self.p - 1) // 2) == 1
    

"""
    SupersingularEllipticCurve - supersingular ellipric curve in Edwards form over GF(p^2), p = 4 3^a 5^b - 1
"""
class SupersingularEllipticCurve:
    # x^2 + y^2 = 1 + d x^2 y^2

    def __init__(self, a, b, d):
        self.a = a
        self.b = b
        
        self.p = 4 * (3**a) * (5**b) - 1
        if type(d) == int:
            self.d = GFp2El(d, 0, self.p)
        else:
            self.d = d

        self.n = (self.p+1)**2

        if self.mul(self.rand_point(), self.n) != EcPoint(1, 0):

            raise RuntimeError("Invalid Curve d={}".format(self.d))
    
    def j(self):
        return (16 * (14*self.d + self.d*self.d + 1).pow(3)) / (self.d * (1-self.d).pow(4))

    def check_on_curve(self, P):
        Q = P.convert_to_affine()
        if (Q.x * Q.x + Q.y * Q.y) == (1 + self.d * Q.x * Q.x * Q.y * Q.y):
            return True
        return False

    def add(self, P, Q):
        if P.is_affine:
            return EcPoint( 
                (P.x * Q.x - P.y * Q.y) / ( 1 - self.d * P.x * Q.x * P.y * Q.y), 
                (P.x * Q.y + P.y * Q.x) / ( 1 + self.d * P.x * Q.x * P.y * Q.y) 
            )
        else:
            return EcPoint(
                P.z * Q.z * ( P.z * P.z * Q.z * Q.z + self.d * P.x * P.y * Q.x * Q.y) * (P.x * Q.x - P.y * Q.y),
                P.z * Q.z * ( P.z * P.z * Q.z * Q.z - self.d * P.x * P.y * Q.x * Q.y) * (P.x * Q.y + P.y * Q.x),
                ( P.z * P.z * Q.z * Q.z + self.d * P.x * P.y * Q.x * Q.y ) * ( P.z * P.z * Q.z * Q.z - self.d * P.x * P.y * Q.x * Q.y)
            )

    def mul(self, P, k):
        Q = EcPoint(1, 0, 1)
        T = P.convert_to_proj()
        for i in range(k.bit_length()):
            if k & (1 << i):
                Q = self.add(Q, T)
                if Q.z == 0:
                    raise RuntimeError("Exceptional point: " + str(Q))
            T = self.add(T, T)
            if T.z == 0:
                raise RuntimeError("Exceptional point: " + str(T))
        return Q

    def rand_point(self):
        x = GFp2El(randint(1, self.p - 1), 0, self.p)
        while True:
            a = (1 - x*x) / (1 - self.d * x*x)
            if a.is_quad_residue():
                return EcPoint(x, a.sqrt())
            x += 1

    def find_order_dummy(self, P):
        for i in range(1, self.n+1):
            if self.mul(P, i).convert_to_affine() == EcPoint(1, 0):
                return i

    def generate_3_a_torsion_point(self, a = None):
        if a == None:
            a = self.a
        while True:
            P = self.rand_point()
            
            P = self.mul(P, 4)
            #print(P)
            #print(self.mul(P, self.n))
            #print(self.find_order_dummy(P))
            P = self.mul(P, 5**self.b)
            P = self.mul(P, 3**(self.a - a))
            
            if all([self.mul(P, 3**k) != EcPoint(1, 0) for k in range(a)]) and self.mul(P, 3**a) == EcPoint(1, 0):
                return P

    def generate_5_b_torsion_point(self, b = None):
        if b == None:
            b = self.b
        while True:
            P = self.rand_point()
            P = self.mul(P, 4)
            P = self.mul(P, 3**self.a)
            P = self.mul(P, 5**(self.b - b))
            if all([self.mul(P, 5**k) != EcPoint(1, 0) for k in range(b)]) and self.mul(P, 5**b) == EcPoint(1, 0):
                return P

    def distortion_map(self, P):
        Q = P.convert_to_affine()
        return EcPoint(Q.x * GFp2El(0, 1, self.p), Q.y.inv())

    def compute_3_isogeny(self, P3):
        P3 = P3.convert_to_affine()
        return (P3.x**8) * (self.d**3)

    def compute_3_isogeny_curve(self, P3):
        d = self.compute_3_isogeny(P3)
        return SupersingularEllipticCurve(self.a, self.b, d)
    
    def compute_5_isogeny(self, P5):
        P5 = P5.convert_to_affine()
        return (P5.x**8) * (self.d**5)

    def compute_5_isogeny_curve(self, P5):
        d = self.compute_5_isogeny(P5)
        return SupersingularEllipticCurve(self.a, self.b, d)

    def evaluate_3_isogeny(self, P3, P):
        P3 = P3.convert_to_affine()
        P = P.convert_to_y_less()
        return EcPoint(
            P.x * ( (P.x**2) - (P3.y**2) * (P.z**2) ), 
            None, 
            (P3.x**2) * P.z * (P.z**2 - self.d * (P3.y**2) * (P.x**2) ) 
        )
    
    def evaluate_5_isogeny(self, P5, P):
        P5_2 = self.mul(P5, 2).convert_to_affine()
        P5 = P5.convert_to_affine()
        P = P.convert_to_y_less()
        return EcPoint(
            P.x * ( P.x**2 - (P5.y**2) * (P.z**2) ) * ( P.x**2 - (P5_2.y**2) * (P.z**2) ),
            None,
            ((P5.x * P5_2.x) ** 2) * P.z * (P.z**2 - self.d * (P5.y**2) * (P.x**2)) * (P.z**2 - self.d * (P5_2.y**2) * (P.x**2))
        )

    def draw_3_isogeny_graph(self):
        edges = []
        vertexes = []
        bfsq = Queue(0)
        bfsq.put(self)
        vertexes.append(self.j())
        
        while not bfsq.empty(): # bfs...
            ec = bfsq.get()
            j0 = ec.j()
            P1 = ec.generate_3_a_torsion_point(1) # <a,0>
            P2 = P1
            while P2 == P1 or P2 == P1.neg():
                P2 = ec.generate_3_a_torsion_point(1)
            
            P3 = ec.add(P1, P2) # <a,b>
            P4 = ec.add(P3, P1) # <2a,b>
            
            #print("j = ", j0)
            #print(P1.convert_to_affine())
            #print(P2.convert_to_affine())
            #print(P3.convert_to_affine())
            #print(P4.convert_to_affine())
            
            kernels = [P1, P2, P3, P4]
            for kernel in kernels:
                isog = ec.compute_3_isogeny_curve(kernel)
                j1 = isog.j()
                print("j1 = ", j1)
                if j1 in vertexes:
                    edges.append((j0, j1))
                elif not j1 in vertexes:
                    vertexes.append(j1)
                    edges.append((j0, j1))
                    bfsq.put(isog)

        return (vertexes, edges)
                    

class EcPoint:

    def __init__(self, *args):
        if len(args) == 2:
            self.x = args[0]
            self.y = args[1]
            self.is_affine = True
            self.y_less = self.y is None

        elif len(args) == 3:
            self.x = args[0]
            self.y = args[1]
            self.z = args[2]
            self.is_affine = False
            self.y_less = self.y is None

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
    
    def convert_to_y_less(self):
        if self.y_less:
            return self
        self.y_less = True
        Q = self.convert_to_proj()
        return EcPoint(Q.x, None, Q.z)

    def clone(self):
        if self.is_affine:
            return EcPoint(self.x, self.y)
        else:
            return EcPoint(self.x, self.y, self.z)
    
    def neg(self):
        if self.is_affine:
            return EcPoint(self.x, -self.y)
        else:
            return EcPoint(self.x, -self.y, self.z) 

    def __str__(self):
        if self.is_affine:
            return "x = {}\ny = {}".format(self.x, self.y)
        else:
            return "X = {}\nY = {}\nZ = {}".format(self.x, self.y, self.z)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        a = self.convert_to_affine()
        b = other.convert_to_affine()
        return a.x == b.x and a.y == b.y

def render_graph(v, e):
    dot = Digraph()
    for vv in v:
        dot.node(vv.__str__())

    for ee in e:
        dot.edge(ee[0].__str__(), ee[1].__str__())

    dot.render("g1", view=True)

def test():
    ec = SupersingularEllipticCurve(2, 1, -1)
    v, e = ec.draw_3_isogeny_graph()
    render_graph(v, e)
    

def old_test():
    aa = GFp2El(9, 0, 59) 
    a = GFp2El(9, 10, 59)
    b = GFp2El(2, 39, 59)

    d = a * a.inv()
    print(a)
    print(a.sqrt())
    print(a.is_quad_residue())

    ec = SupersingularEllipticCurve(8, 3, -1)
        
    P3 = ec.generate_3_a_torsion_point(1)
    P5 = ec.generate_5_b_torsion_point(1)

    print(ec.j())
    print(P3.convert_to_affine())
    print(P5.convert_to_affine())

    P3_d = ec.mul(ec.distortion_map(P3), 4)
    P5_d = ec.mul(ec.distortion_map(P5), 4)

    print(P3_d.convert_to_affine())
    print(P5_d.convert_to_affine())

    print(ec.find_order_dummy(P3))
    print(ec.find_order_dummy(P3_d))
    print(ec.find_order_dummy(P5))
    print(ec.find_order_dummy(P5_d))
    
    print(ec.compute_3_isogeny(P3))
    print(ec.compute_5_isogeny(P3_d))

if __name__ == "__main__":
    test()
    