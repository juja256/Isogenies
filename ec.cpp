#include "ec.h"
#include "gf.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sstream>

PseudoRandomGenerator::PseudoRandomGenerator(unsigned char* seed, int len) {
    srand(time(0));
    memcpy(this->state, seed, len);
}

void PseudoRandomGenerator::Run() {
    int a = rand();
    *((int*)state) = a;
}
void PseudoRandomGenerator::DeriveRandomFromState(char* ptr, int byteCnt) {
    memcpy(ptr, state, byteCnt);
}

void PseudoRandomGenerator::GenerateSequence(int bitLen, char* dest) {
    int byteCnt = bitLen / 8;
    int restBits = bitLen % 8;
    for (int i=0; i<byteCnt; i++) {
        Run();
        DeriveRandomFromState(dest + i, 1);
    }
}

EllipticCurveException::EllipticCurveException(int c): code(c) {}
std::string EllipticCurveException::What() {
    std::stringstream ss;
    ss << "EllipticCurveException #";
    ss << code;
    return ss.str();
}

EllipticCurve::EllipticCurve() {}

EllipticCurve::EllipticCurve(PseudoRandomGenerator* p): prng(p), isBasePointPresent(false) {}

void EllipticCurve::InitAsWeierstrass(GaloisField* GF, const BigInt cardinality, const GFElement a, const GFElement b, const EcPoint* BP = NULL) {
    this->GF = GF;
    this->form = WEIERSTRASS;
    memcpy(this->n, cardinality, 2*BYTES_IN_WORD*GF->GetWordSize());
    this->GF->GFCopy(this->a, a);
    this->GF->GFCopy(this->b, b);
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
    this->GF->GFCopy(dst->X, src->X);
    this->GF->GFCopy(dst->Y, src->Y);
}

void EllipticCurve::PointCopy(EcPointProj* dst, const EcPointProj* src) {
    this->GF->GFCopy(dst->X, src->X);
    this->GF->GFCopy(dst->Y, src->Y);
    this->GF->GFCopy(dst->Z, src->Z);
}

void EllipticCurve::InitAsEdwards(GaloisField* GF, const BigInt cardinality, const GFElement d, const EcPoint* BP = NULL) {
    this->GF = GF;
    this->form = EDWARDS;
    memcpy(this->n, cardinality, 2*BYTES_IN_WORD*GF->GetWordSize());
    this->GF->GFCopy(this->d, a);
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
    divide(2*GF->GetWordSize(), this->n, *(this->GF->GetChar()), NULL, r); // n = 1 (mod p)
    return GF->GFBaseCmp(r, GF->Unity[0]) == 0;
}

void EllipticCurve::GenerateBasePoint() {
    // ?
}

void EllipticCurve::GetJInvariant(GFElement J) {
    
    switch (form) {
        case EDWARDS: // 16(1 + 14*d + d^2)^3 / d(1-d)^4
            GFElement t,g,c,e;
            this->GF->GFCopy(t, GF->Zero);
            t[0][0] = 14;
            GF->GFMul(t, d, t);
            GF->GFAdd(t, GF->Unity, t);
            GF->GFSqr(d, g);
            GF->GFAdd(t, g, t);
            GF->GFSqr(t, e);
            GF->GFMul(t, e, t);
            GF->GFMulBy2Power(t, 4, t);

            GF->GFSub(GF->Unity, d, c);
            GF->GFSqr(c, c);
            GF->GFSqr(c, c);
            GF->GFMul(c, d, c);
            GF->GFInv(c, e);

            GF->GFMul(e, t, J);
            break;
        case WEIERSTRASS: // 6912 a^3 / (4 a^3 + 27 b^2)
            GFElement c,e,t,g;
            GF->GFSqr(a, c);
            GF->GFMul(c, a, c); // a ^ 3
            e[0][0] = 6912;
            GF->GFMul(e, c, e);
            
            GF->GFMulBy2Power(c, 2, c);
            t[0][0] = 27;
            GF->GFSqr(b, g);
            GF->GFMul(g, t, g);
            GF->GFAdd(g, c, g);
            GF->GFInv(g, t);

            GF->GFMul(e, t, J);
            break;
    }
}

bool EllipticCurve::IsPointOnCurve(const EcPoint* P) {
    switch (form) {
        case EDWARDS:
            GFElement r,t,f,g;
            GF->GFSqr(P->X, r);
            GF->GFSqr(P->Y, t);
            GF->GFAdd(r, t, f);
            GF->GFMul(r, t, g);
            GF->GFMul(g, d, g);
            GF->GFSub(f, g, f);

            return GF->GFCmp(f, GF->Unity);
        case WEIERSTRASS:
            GFElement r,t,f,g;
            GF->GFSqr(P->Y, r);
            GF->GFSqr(P->X, t);
            GF->GFMul(t, P->X, t);
            GF->GFMul(P->X, a, f);
            GF->GFAdd(t, f, t);
            GF->GFSub(r, t, r);

            return GF->GFCmp(r, b);
    }
}