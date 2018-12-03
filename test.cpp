#include "ec.h"
#include <stdio.h>
#include <iostream>
int main() {
    PseudoRandomGenerator prng;
    EllipticCurve EC(&prng);
    BigInt p;
    BigInt n;
    GFElement d, J, t1, t2;
    GaloisField::InitFromString(p, "6FE5D541F71C0E12909F97BADC668562B5045CB25748084E9867D6EBE876DA959B1A13F7CC76E3EC968549F878A8EEAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
    GaloisField::InitFromString(d, "6FE5D541F71C0E12909F97BADC668562B5045CB25748084E9867D6EBE876DA959B1A13F7CC76E3EC968549F878A8EEAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE, 00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    GaloisField::InitFromString(n, "30E91D466DF5429960D2536B6AE0D99AA4835FED951F1D323FB4C115170A25E037E40347E3AB06E7A12F5FF7B20AD617C8DF437CFA421554FE2E49CA85BAB790796CF84D4D74319A6C9BCA37551AE9F5F5F6AFF3B7F731B89B2DA43F258BB9000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    GaloisField GF(p, 2, 752);

    GF.Sqr(d, t1);
    //std::cout << GF.Dump(d) << "\n";
    //std::cout << GF.Dump(t1) << "\n";

    GF.Copy(t1, GF.I);
    t1[0][0] = 1;
    GF.Copy(t2, t1);
    GF.BaseNeg(t2[1], t2[1]);
    std::cout << GF.Dump(t1) << " " << GF.IsQuadraticResidue(t1) << "\n";
    std::cout << GF.Dump(t2) << " " << GF.IsQuadraticResidue(t2) << "\n";
    GF.Mul(t1, t2, t2);
    std::cout << GF.Dump(t2) << "\n";

    GF.Pow(t1, *(GF.GetSize()), t2);
    std::cout << GF.Dump(t2) << "\n";
    EcPoint D;
    
    try {
        EC.InitAsEdwards(&GF, n, d, NULL);
        EC.PointCopy(&D, &(EC.BasePoint));
        
        std::cout << "Base Point: " << EC.PointDump(&D) << "\n";
        std::cout << "Base Point Valid: " << EC.IsPointOnCurve(&D) << " " ;
        std::cout << EC.CheckPointTorsion(&D, p) << "\n";
        EC.GetJInvariant(J);
        std::cout << "J invariant: " << GF.Dump(J) << "\n";
    } catch(GaloisFieldException& e) {
        printf(e.What().c_str());
    } catch(EllipticCurveException& e) {
        printf(e.What().c_str());
    }
    

   
    return 0;
}