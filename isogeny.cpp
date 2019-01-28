#include "isogeny.h"

SIDHEngine::SIDHEngine(int l2, int l3, int f) : l2(l2), l3(l3), f(f) {
}

SIDHEngine::SIDHEngine(int l2, int l3, int f, BigInt p, int bit_size) : l2(l2), l3(l3), f(f) {
    GF = new GaloisField(p, 2, bit_size);
    unsigned char seed[4] = { 0x05, 0x05, 0x05, 0x05 };
    prng = new PseudoRandomGenerator((unsigned char*)seed, 4);
    BaseCurve = new EllipticCurve(prng);
    BigInt cardinality;
    add(GF->GetWordSize(), *(GF->GetChar()), GF->Unity[0], cardinality);
    sqr(2*GF->GetWordSize(), cardinality, cardinality);
    GFElement d;
    GF->Neg(GF->Unity, d);
    BaseCurve->InitAsEdwards(GF, cardinality, d, NULL);

    shl(GF->GetWordSize(), GF->Unity[0], L2, l2);

    BaseEl a3;
    GF->BaseCopy(L3, GF->Unity[0]);
    for (int i=0; i<l3; i++) {
        GF->BaseCopy(a3, L3);
        shl(GF->GetWordSize(), a3, L3, 1);
        add(GF->GetWordSize(), L3, a3, L3);
    }
}

SIDHEngine::~SIDHEngine() {
    delete prng;
    delete GF;
    delete BaseCurve;
}

// 1 - x^2y^2 = x^2 + y^2
void SIDHEngine::Generate4LTorsionPoint(EcPoint* P4L) {
    BaseEl x, x2, y, y2;
    GF->Copy(P4L->X, GF->Unity);
    P4L->X[0][0] = 2; // X = 2
    GF->Copy(P4L->Y, GF->Zero); // Y = 0
    for (;;) {
        GF->BaseSqr(P4L->X[0], x2);
        GF->BaseSub(x2, GF->Unity[0], y);
        GF->BaseMul(x2, BaseCurve->d[0], x2);
        GF->BaseSub(x2, GF->Unity[0], x2);
        GF->BaseInv(x2, x2);
        GF->BaseMul(x2, y, y);
        if (GF->IsQuadraticResidueBase(y)) {
            EcPoint J;
            GF->BaseSqrt(y, P4L->Y[0]);
            BaseCurve->ScalarMul(P4L, L3, &J);
            if (BaseCurve->PointEqual(&J, &BaseCurve->UnityPoint)) {

            }
        }
        P4L->X[0][0] ++;
    }
}

void SIDHEngine::Generate3LTorsionPoint(EcPoint* P3L) {

}

void SIDHEngine::GenerateBasePoints() {
    Generate3LTorsionPoint(&P_A);
    Generate4LTorsionPoint(&P_B);
    BaseCurve->ApplyDistortionMap(&P_A, &Q_A);
    BaseCurve->ApplyDistortionMap(&P_B, &Q_B);
}

void SIDHEngine::SetBasePoints(const EcPoint* p_A, const EcPoint* q_A, const EcPoint* p_B, const EcPoint* q_B) {
    BaseCurve->PointCopy(&P_A, p_A);
    BaseCurve->PointCopy(&Q_A, q_A);
    BaseCurve->PointCopy(&P_B, p_B);
    BaseCurve->PointCopy(&Q_B, q_B);   
}

void SIDHEngine::Compute3Isogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, GFElement D) {
    GFElement c0, c1, c2, c3, t0, t1, C;
    E->GF->Add(kernelPoint->Y, kernelPoint->Z, c0);
    E->GF->Sqr(c0, c0);
    E->GF->Sqr(kernelPoint->Y, c1);
    E->GF->Sqr(kernelPoint->Z, c2);
    E->GF->Sub(c0, c1, c3);
    E->GF->Sub(c3, c2, c3);
    E->GF->Add(c1, c1, t0);
    E->GF->Add(t0, c1, t0);
    E->GF->Add(c0, c3, t1);
    E->GF->Add(t1, t0, t1);
    E->GF->Add(c2, c3, t0);
    E->GF->Mul(t0, t1, D);

    E->GF->Sub(c1, c2, t0);
    E->GF->Add(c3, c3, t1);
    E->GF->Sub(c0, t1, t1);
    E->GF->Mul(t0, t1, t0);
    E->GF->Add(t0, D, C);

    E->GF->Inv(C, C);
    E->GF->Mul(C, D, D);
}

void SIDHEngine::Evaluate3Isogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, const EcPointProj* P, EcPointProj* PImage) {
    GFElement t0, t1, t2, t3;
    E->GF->Mul(kernelPoint->Y, P->Z, t0);
    E->GF->Mul(kernelPoint->Z, P->Y, t1);
    E->GF->Add(t0, t1, t2);
    E->GF->Add(P->Y, P->Z, t3);
    E->GF->Sqr(t2, t2);
    E->GF->Mul(t2, t3, t2);
    E->GF->Sub(t0, t1, t0);
    E->GF->Sqr(t0, t0);
    E->GF->Sub(P->Y, P->Z, t1);
    E->GF->Mul(t0, t1, t0);
    E->GF->Add(t0, t2, PImage->Y);
    E->GF->Sub(t2, t0, PImage->Z);
}

void SIDHEngine::Compute4Isogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, GFElement D) {
    GFElement c0, C;
    E->GF->Add(kernelPoint->Y, kernelPoint->Z, c0);
    E->GF->Sqr(c0, c0);
    E->GF->Sqr(c0, C);
    E->GF->Sub(kernelPoint->Y, kernelPoint->Z, c0);
    E->GF->Sqr(c0, c0);
    E->GF->Sqr(c0, c0);
    E->GF->Sub(C, c0, D);

    E->GF->Inv(C, C);
    E->GF->Mul(C, D, D);
}

void SIDHEngine::Evaluate4Isogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, const EcPointProj* P, EcPointProj* PImage) {
    GFElement t0, t1, t2, t3, t4, t5, c0;
    E->GF->Mul(P->Z, kernelPoint->Y, t0);
    E->GF->Mul(P->Y, kernelPoint->Z, t1);
    E->GF->Mul(t0, t1, t2);
    E->GF->Add(t2, t2, t2);
    E->GF->Add(t0, t1, t3);
    E->GF->Sqr(t3, t3);
    E->GF->Sub(t3, t2, t5);
    E->GF->Mul(P->Y, P->Z, t4);
    E->GF->Add(kernelPoint->Y, kernelPoint->Z, c0);
    E->GF->Sqr(c0, c0);
    E->GF->Mul(t4, c0, t4);
    E->GF->Sub(t4, t2, t4);
    E->GF->Add(t4, t5, PImage->Y);
    E->GF->Sub(t4, t5, PImage->Z);
    E->GF->Mul(PImage->Y, t3, t0);
    E->GF->Sub(t5, t2, t3);
    E->GF->Mul(PImage->Z, t3, t1);
    E->GF->Add(t0, t1, PImage->Y);
    E->GF->Sub(t0, t1, PImage->Z);
}

void ComputeAndEvaluate3LIsogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, const EcPointProj* P, const EcPointProj* Q, EcPointProj* PImage, EcPointProj* QImage, EllipticCurve* F) {

}

void ComputeAndEvaluate4LIsogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, const EcPointProj* P, const EcPointProj* Q, EcPointProj* PImage, EcPointProj* QImage, EllipticCurve* F) {

}