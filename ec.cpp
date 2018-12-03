#include "ec.h"
#include "gf.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sstream>
#include <iostream>

PseudoRandomGenerator::PseudoRandomGenerator(unsigned char* seed, int len) {
    srand(time(0));
    memcpy(this->state, seed, len);
}

PseudoRandomGenerator::PseudoRandomGenerator() {

}

void PseudoRandomGenerator::Run() {
    int a = rand();
    *((int*)state) = a;
}
void PseudoRandomGenerator::DeriveRandomFromState(unsigned char* ptr, int byteCnt) {
    memcpy(ptr, state, byteCnt);
}

void PseudoRandomGenerator::GenerateSequence(int bitLen, unsigned char* dest) {
    int byteCnt = bitLen / 8;
    int restBits = bitLen % 8;
    for (int i=0; i<byteCnt; i++) {
        Run();
        DeriveRandomFromState(dest + i, 1);
    }
}

void PseudoRandomGenerator::GenerateBaseEl(int bit_len, BaseEl a) {
    a[bit_len/ARCH] = 0;
    GenerateSequence(bit_len, (unsigned char*)a);
}

EllipticCurveException::EllipticCurveException(int c): code(c) {}
std::string EllipticCurveException::What() {
    std::stringstream ss;
    ss << "EllipticCurveException #";
    ss << code;
    return ss.str();
}

EllipticCurve::EllipticCurve() : isBasePointPresent(false), scalarMulEngine(&EllipticCurve::ScalarMulMontgomery) {}

EllipticCurve::EllipticCurve(PseudoRandomGenerator* p): prng(p), isBasePointPresent(false), scalarMulEngine(&EllipticCurve::ScalarMulMontgomery) {}

void EllipticCurve::InitAsWeierstrass(GaloisField* GF, const BigInt cardinality, const GFElement a, const GFElement b, const EcPoint* BP) {
    this->GF = GF;
    this->form = WEIERSTRASS;
    memcpy(this->n, cardinality, 2*BYTES_IN_WORD*GF->GetWordSize());
    this->GF->Copy(this->a, a);
    this->GF->Copy(this->b, b);
    if (BP != NULL) {
        PointCopy(&(this->BasePoint), BP);
    } else {
        GenerateBasePoint();
    }
    this->isSupersingular = CheckSupersingularity();

    this->UnityPoint = { {{0}, {0}}, {{0}, {0}} };
    this->UnityPointProj = { {{0}, {0}}, {{1}, {0}}, {{0}, {0}} };
}

void EllipticCurve::PointCopy(EcPoint* dst, const EcPoint* src) {
    this->GF->Copy(dst->X, src->X);
    this->GF->Copy(dst->Y, src->Y);
}

void EllipticCurve::PointCopy(EcPointProj* dst, const EcPointProj* src) {
    this->GF->Copy(dst->X, src->X);
    this->GF->Copy(dst->Y, src->Y);
    this->GF->Copy(dst->Z, src->Z);
}

bool EllipticCurve::PointEqual(const EcPoint* A, const EcPoint* B) {
    return GF->Equal(A->X, B->X) && GF->Equal(A->Y, B->Y);
}

void EllipticCurve::ToProjective(const EcPoint* src, EcPointProj* dst) {
    GF->Copy(dst->Z, GF->Zero);
    prng->GenerateSequence(GF->GetBitSize()-1, (unsigned char*)dst->Z[0]); // Z from GF(p)
    GF->MulByBase(src->X, dst->Z[0], dst->X);
    GF->MulByBase(src->Y, dst->Z[0], dst->Y);
}

#define UNDEFINED_POINT -45

void EllipticCurve::ToAffine(const EcPointProj* src, EcPoint* dst) {
    GFElement z_inv;
    #ifdef DEBUG
    if (GF->Equal(src->Z, GF->Zero)) {  
        std::cout << "UNDEFINED POINT: \n" << PointDump(src);
    }
    #endif
    GF->Inv(src->Z, z_inv);
    GF->Mul(src->X, z_inv, dst->X);
    GF->Mul(src->Y, z_inv, dst->Y);
}

void EllipticCurve::InitAsEdwards(GaloisField* GF, const BigInt cardinality, const GFElement d, const EcPoint* BP) {
    this->GF = GF;
    this->form = EDWARDS;
    memcpy(this->n, cardinality, 2*BYTES_IN_WORD*GF->GetWordSize());
    this->GF->Copy(this->d, d);
    if (BP != NULL) {
        PointCopy(&(this->BasePoint), BP);
    } else {
        GenerateBasePoint();
    }
    this->isSupersingular = CheckSupersingularity();

    this->UnityPoint = { {{0}, {0}}, {{1}, {0}} };
    this->UnityPointProj = { {{0}, {0}}, {{1}, {0}}, {{1}, {0}} };
}

void EllipticCurve::SetPseudoRandomProvider(PseudoRandomGenerator* p) {
    this->prng = p;
}

bool EllipticCurve::CheckSupersingularity() {
    BigInt r;
    divide(GF->GetWordSize(), this->n, *(this->GF->GetChar()), NULL, r); // n = 1 (mod p)
    return GF->BaseCmp(r, GF->Unity[0]) == 0;
}

void EllipticCurve::GetJInvariant(GFElement J) {
    GFElement t,g,c,e;
    switch (form) {
        case EDWARDS: // 16(1 + 14*d + d^2)^3 / d(1-d)^4
            
            this->GF->Copy(t, GF->Zero);
            t[0][0] = 14;
            GF->Mul(t, d, t);
            GF->Add(t, GF->Unity, t);
            GF->Sqr(d, g);
            GF->Add(t, g, t);
            GF->Sqr(t, e);
            GF->Mul(t, e, t);
            GF->MulBy2Power(t, 4, t);

            GF->Sub(GF->Unity, d, c);
            GF->Sqr(c, c);
            GF->Sqr(c, c);
            GF->Mul(c, d, c);
            GF->Inv(c, e);

            GF->Mul(e, t, J);
            break;
        case WEIERSTRASS: // 6912 a^3 / (4 a^3 + 27 b^2)

            GF->Sqr(a, c);
            GF->Mul(c, a, c); // a ^ 3
            e[0][0] = 6912;
            GF->Mul(e, c, e);
            
            GF->MulBy2Power(c, 2, c);
            t[0][0] = 27;
            GF->Sqr(b, g);
            GF->Mul(g, t, g);
            GF->Add(g, c, g);
            GF->Inv(g, t);

            GF->Mul(e, t, J);
            break;
    }
}

bool EllipticCurve::IsPointOnCurve(const EcPoint* P) {
    GFElement r,t,f,g;
    switch (form) {
        case EDWARDS:
            
            GF->Sqr(P->X, r);
            GF->Sqr(P->Y, t);
            GF->Add(r, t, f);
            GF->Mul(r, t, g);
            GF->Mul(g, d, g);
            GF->Sub(f, g, f);

            return GF->Equal(f, GF->Unity);
        case WEIERSTRASS:

            GF->Sqr(P->Y, r);
            GF->Sqr(P->X, t);
            GF->Mul(t, P->X, t);
            GF->Mul(P->X, a, f);
            GF->Add(t, f, t);
            GF->Sub(r, t, r);

            return GF->Equal(r, b);
    }
}

bool EllipticCurve::CheckPointTorsion(const EcPoint* P, const BigInt order) {
    EcPoint Z;
    ScalarMul(P, order, &Z);
    return PointEqual(&Z, &UnityPoint);
}

void EllipticCurve::AcquireEdwardsForm() {
    if (form != EDWARDS) throw EllipticCurveException(NOT_SUPPORTED);
}

// In my world (0,1) is an Edwards unity point
void EllipticCurve::Add(const EcPoint* A, const EcPoint* B, EcPoint* C) {
    AcquireEdwardsForm();
    GFElement z1, z2, z3, z4, z5, z6, z7;
    GF->Mul(A->X, B->Y, z1); // z1 = x1 * x2
    GF->Mul(A->Y, B->Y, z2); // z2 = y1 * y2
    GF->Mul(z1, z2, z3); 
    GF->Mul(z3, d, z3); // z3 = d * x1 * x2 * y1 * y2

    GF->Neg(z3, z4); // z4 = - z3
    GF->Add(z3, GF->Unity, z3);

    GF->Inv(z3, z3);

    GF->Add(z4, GF->Unity, z4);

    GF->Inv(z4, z4);

    GF->Mul(A->X, B->Y, z5); // z5 = x1 * y2
    GF->Mul(A->Y, B->X, z6); // z6 = x2 * y1
    GF->Add(z5, z6, z5);
    GF->Sub(z2, z1, z2);

    GF->Mul(z5, z3, C->X);
    GF->Mul(z2, z4, C->Y);
}

void EllipticCurve::Dbl(const EcPoint* A, EcPoint* B) {
    AcquireEdwardsForm();
    GFElement z1, z2, z3, z4, z5;
    GF->Sqr(A->X, z1);
    GF->Sqr(A->Y, z2);
    GF->Mul(A->X, A->Y, z3);
    GF->Add(z3, z3, z3);
    GF->Mul(z1, z2, z4);
    GF->Mul(z4, d, z4);
    GF->Neg(z4, z5);
    GF->Add(z4, GF->Unity, z4);
    GF->Inv(z4, z4);
    GF->Add(z5, GF->Unity, z5);
    GF->Inv(z5, z5);
    GF->Sub(z2, z1, z2);

    GF->Mul(z4, z3, B->X);
    GF->Mul(z2, z5, B->Y);
}

void EllipticCurve::Add(const EcPointProj* P1, const EcPointProj* P2, EcPointProj* P3) {
    AcquireEdwardsForm();
    #ifdef DEBUG
    EcPointProj P1_copy, P2_copy;
    PointCopy(&P1_copy, P1);
    PointCopy(&P2_copy, P2);
    #endif

    /* 10M + 1S */
    GFElement A, B, C, D, E, F, G, T;
    GF->Mul(P1->Z, P2->Z, A); // A = Z1Z2
    GF->Sqr(A, B); // B = A^2
    GF->Mul(P1->Y, P2->Y, C); // C = Y1Y2
    GF->Mul(P1->X, P2->X, D); // D = X1X2
    GF->Mul(C, D, E); 
    GF->Mul(E, d, E); // E = dCD
    GF->Sub(B, E, F); // F = B-E
    GF->Add(B, E, G); // G = B+E
    
    GF->Add(P1->Y, P1->X, T);
    GF->Add(P2->Y, P2->X, P3->X);
    GF->Mul(P3->X, T, P3->X);
    GF->Sub(P3->X, C, P3->X);
    GF->Sub(P3->X, D, P3->X);
    GF->Mul(P3->X, A, P3->X);
    GF->Mul(P3->X, F, P3->X); // X3 = AF((X1+Y1)(X2+Y2)-C-D)

    GF->Sub(C, D, P3->Y);
    GF->Mul(P3->Y, A, P3->Y);
    GF->Mul(P3->Y, G, P3->Y); // Y3 = AG(C-D) 

    GF->Mul(F, G, P3->Z); // Z3 = FG

    #ifdef DEBUG
    if (GF->Equal(P3->Z, GF->Zero)) {
        std::cout << "!ADD, Z=0: \n" << PointDump(&P1_copy) << "\n" << PointDump(&P2_copy) << "\n" << PointDump(P3);
    }
    #endif
}

void EllipticCurve::Dbl(const EcPointProj* P, EcPointProj* P2) {
    AcquireEdwardsForm();
    #ifdef DEBUG
    EcPointProj P_copy;
    PointCopy(&P_copy, P);
    #endif
    /* 3M + 4S */
    GFElement A,B,C,D,E,F,G;

    GF->Sqr(P->Y, A); // A = Y^2
    GF->Sqr(P->X, B); // B = X^2
    GF->Sqr(P->Z, C); // C = Z^2
    GF->Add(A, B, D); // D = A+B
    GF->Sub(A, B, E); // E = A-B
    GF->Add(C, C, F); 
    GF->Sub(F, D, F); // F = 2C - A - B
    GF->Add(P->Y, P->X, G);
    GF->Sqr(G, G);
    GF->Sub(G, D, G);
    
    GF->Mul(F, G, P2->X); // X2 = FG 
    GF->Mul(D, E, P2->Y); // Y2 = DE
    GF->Mul(D, F, P2->Z); // Z3 = DF 
    #ifdef DEBUG
    if (GF->Equal(P2->Z, GF->Zero)) {
        std::cout << "!DBL, Z=0: \n" << PointDump(&P_copy) << "\n" << PointDump(P2);
        std::cout << "!D: " << GF->Dump(D) << "\n";
        std::cout << "!F: " << GF->Dump(F) << "\n";
    }
    #endif
}

void EllipticCurve::ScalarMulNaive(const EcPointProj* P, const BigInt k, EcPointProj* Q, int bitLen) {
    EcPointProj H;
    PointCopy(&H, P); // H := A
    PointCopy(Q, &UnityPointProj);
    bitLen = (bitLen != 0) ? bitLen : GF->GetBitSize()*GF->GetExtension();
    for (u32 i=0; i<bitLen; i++) {
        if (get_bit(k, i)) {
            Add(Q, &H, Q);
        }
        Dbl(&H, &H); 
    }
}
void EllipticCurve::ScalarMul(const EcPoint* P, const BigInt k, EcPoint* Q, int bitLen) {
    EcPointProj H;
    ToProjective(P, &H);
    ScalarMul(&H, k, &H, bitLen);
    ToAffine(&H, Q);
}

void EllipticCurve::ScalarMul(const EcPointProj* P, const BigInt k, EcPointProj* Q, int bitLen) {
    (this->*scalarMulEngine)(P, k, Q, bitLen);
}

void EllipticCurve::ScalarMulMontgomery(const EcPointProj* P, const BigInt k, EcPointProj* Q, int bitLen) {
    EcPointProj H;
    PointCopy(&H, P); // H := A, H = P1
    PointCopy(Q, &UnityPointProj);
    bitLen = (bitLen != 0) ? bitLen : GF->GetBitSize()*GF->GetExtension();

    for (int i=bitLen-1; i>=0; i--) {
        if (get_bit(k, i) == 0) {
            Add(Q, &H, &H);
            Dbl(Q, Q); 
        }
        else {
            Add(Q, &H, Q);
            Dbl(&H, &H); 
        }
    }
}

// ksi(x,y) = (xi, y^{-1})
void EllipticCurve::ApplyDistortionMap(const EcPoint* P, EcPoint* Q) {
    GF->Mul(P->X, GF->I, Q->X);
    GF->Inv(P->X, Q->Y);
}

EllipticCurve::~EllipticCurve() {

}

#define NOT_QUAD_RESIDUE -44

void EllipticCurve::GenerateBasePoint() {
    BaseEl x2;
    GFElement y;
    GF->Copy(y, GF->Zero);

    prng->GenerateBaseEl(GF->GetBitSize()-1, BasePoint.X[0]);
    GF->BaseCopy(BasePoint.X[1], GF->Zero[1]);

    GF->BaseSqr(BasePoint.X[0], x2);
    GF->BaseSub(x2, GF->Unity[0], y[0]);

    GF->BaseMul(x2, d[0], x2);
    GF->BaseSub(x2, GF->Unity[0], x2);
    GF->BaseInv(x2, x2);
    GF->BaseMul(x2, y[0], y[0]);
    #ifdef DEBUG
    if (!GF->IsQuadraticResidue(y)) throw GaloisFieldException(NOT_QUAD_RESIDUE);
    #endif
    GF->Sqrt(y, BasePoint.Y);

}

std::string EllipticCurve::PointDump(const EcPoint* X) {
    return  "x: " + GF->Dump(X->X) + "\n" + "y: " + GF->Dump(X->Y) + "\n";
}

std::string EllipticCurve::PointDump(const EcPointProj* X) {
    return  "X: " + GF->Dump(X->X) + "\n" + "Y: " + GF->Dump(X->Y) + "\n" + "Z: " + GF->Dump(X->Z) + "\n";
}
