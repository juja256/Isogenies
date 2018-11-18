#ifndef EC_H
#define EC_H

#define NORMAL_POINT 1
#define INFINITY_POINT 0

#include <string>
#include "gf.h"
/* 
    In SIDH we use prime p = 2^372 * 3^239 - 1 and Field F_p^2 for 124bit quantum security level
    p = 3 (mod 4)
    Resulting supersingular elliptic curve is E/F_{p^2} : y^2 = x^3 + x; #E = (2^237 * 3^239)^2
    
    For all x in F_p^2 x = a + i*b, where i is square root of QNR in F_p 
*/

typedef struct {
    GFElement X;
    GFElement Y;
} EcPoint;

typedef struct {
    GFElement X;
    GFElement Y;
    GFElement Z;
} EcPointProj;

#define PRNG_STATE_LEN 256 

/* Dummy realization of PRNG, new stronger implementations must be derived from this */
class PseudoRandomGenerator {
protected:
    unsigned char state[PRNG_STATE_LEN];
    virtual void Run(); 
    virtual void DeriveRandomFromState(unsigned char* ptr, int byteCnt);
public:
    PseudoRandomGenerator(unsigned char* seed, int len);
    virtual void GenerateSequence(int bit_len, unsigned char* dest);
    virtual void GenerateBaseEl(int bit_len, BaseEl a);
};

#define WEIERSTRASS 0
#define MONTHOMERRY 1
#define EDWARDS 2

#define NOT_SUPPORTED -1

class EllipticCurveException {
    int code;
public:
    EllipticCurveException(int c);
    std::string What();
};

/* 
Elliptic curve in Edwards:
    x^2 + y^2 = 1 + d*x^2*y^2
   or in Weierstrass form:
    y^2 = x^3 + a*x + b 
    
Scalar Multiplications(all in the projective coordinates):
    - AddAndDouble naive:
    Result: Q = kP
    Q <- 0P
    H <- 1P
    for i from 0 to m do
        if k[i] = 1
            Q <- Q + H
        H <- 2H

    - Montgomery constant-time:
    Result: Q = kP
    Q <- 0P
    H <- 1P
    for i from 0 to m do
        if k[i] = 0
            H <- H + Q
            Q <- 2Q
        else 
            Q <- H + Q
            H <- 2H 
    - Windowed constant time for Edwards:
    Precomputation: 0P, 1P, ... (2^w - 1)P
    Result: Q = kP
    Q <- 0P
    for i from m/w downto 0 do
        Q <- (2^w)Q
        Q <- Q + k[i]P

    - wNAF(not tested) 
*/

typedef void TScalarMul(const EcPointProj*, const BigInt, EcPointProj*, int);

class EllipticCurve {
    GaloisField* GF;
    u8 form;
    bool isSupersingular;

    GFElement d, a, b; // d for Edwards form; a,b for Weierstrass form
    BigInt n; // cardinality
    EcPoint BasePoint;
    bool isBasePointPresent;

    PseudoRandomGenerator* prng;
    EcPointProj* T; // for precomputations

    TScalarMul* scalarMulEngine;

    void AcquireEdwardsForm();
    void ScalarMulNaive(const EcPointProj*, const BigInt, EcPointProj*, int bitLen=0);
    void ScalarMulMontgomery(const EcPointProj*, const BigInt, EcPointProj*, int bitLen=0);
public:
    EcPoint UnityPoint;
    EcPointProj UnityPointProj; 
    EcPoint BasePoint;

    EllipticCurve();
    EllipticCurve(PseudoRandomGenerator*);
    ~EllipticCurve();
    void InitAsWeierstrass(GaloisField* GF, const BigInt cardinality, const GFElement a, const GFElement b, const EcPoint* BP = NULL);
    void InitAsEdwards(GaloisField* GF, const BigInt cardinality, const GFElement d, const EcPoint* BP = NULL);
    void SetPseudoRandomProvider(PseudoRandomGenerator* p);

    bool CheckSupersingularity();
    void GenerateBasePoint();
    void GetJInvariant(GFElement J);

    bool IsPointOnCurve(const EcPoint* P);
    bool CheckPointTorsion(const EcPoint* P, const BigInt order); 
    
    bool PointEqual(const EcPoint* A, const EcPoint* B);
    bool PointEqual(const EcPointProj* A, const EcPointProj* B);

    void PointCopy(EcPoint* dst, const EcPoint* src);
    void PointCopy(EcPointProj* dst, const EcPointProj* B);

    void ToProjective(const EcPoint* src, EcPointProj* dst);
    void ToAffine(const EcPointProj* src, EcPoint* dst);

    void Add(const EcPoint* A, const EcPoint* B, EcPoint* C);
    void Add(const EcPointProj* A, const EcPointProj* B, EcPointProj* C);

    void Dbl(const EcPoint* A, EcPoint* C);
    void Dbl(const EcPointProj* A, EcPointProj* C);

    void SetNaiveScalarMulEngine(); // DoubleAndAdd algorithm
    void SetMontgomeryScalarMulEngine(); // Suitable for cryptologic usage
    //void SetScalarMulWindowedEngine(const EcPoint* A, int windowSize); // Very fast and still suitable for cryptology in Edwards form

    void ScalarMul(const EcPoint* P, const BigInt k, EcPoint* Q, int bitLen=0);
    void ScalarMul(const EcPointProj* P, const BigInt k, EcPointProj* Q, int bitLen=0);
    
};


#endif /* EC_H */