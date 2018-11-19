#include "isogeny.h"

void IsogenyEngine::Compute3Isogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, EllipticCurve* F) {
    GFElement c0, c1, c2, c3, t0, t1, C, D;
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
    F->InitAsEdwards(E->GF, E->n, D, NULL);
}

void IsogenyEngine::Evaluate3Isogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, const EcPointProj* P, EcPointProj* PImage) {
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

void IsogenyEngine::Compute4Isogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, EllipticCurve* F) {
    GFElement c0, C, D;
    E->GF->Add(kernelPoint->Y, kernelPoint->Z, c0);
    E->GF->Sqr(c0, c0);
    E->GF->Sqr(c0, C);
    E->GF->Sub(kernelPoint->Y, kernelPoint->Z, c0);
    E->GF->Sqr(c0, c0);
    E->GF->Sqr(c0, c0);
    E->GF->Sub(C, c0, D);

    E->GF->Inv(C, C);
    E->GF->Mul(C, D, D);
    F->InitAsEdwards(E->GF, E->n, D, NULL);
}

void IsogenyEngine::Evaluate4Isogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, const EcPointProj* P, EcPointProj* PImage) {
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