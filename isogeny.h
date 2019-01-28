#ifndef ISOGENY_H
#define ISOGENY_H

#include "ec.h"

#define ALICE_KEY 3
#define BOB_KEY 4

typedef struct {
    BigInt d; // param of Edwards curve
    EcPoint P_img;
    EcPoint Q_img;
} SIDHPublicKey;

typedef struct {
    BigInt r;
    EcPoint R;
    int type; // Alice or Bob
} SIDHPrivateKey;

class SIDHEngine {
    /* 
        Parameters identifying group structure of E over field GF(p^2), p = (2^l2 * 3^l3 * f)-1 prime; 
        E - supersingular |E| = (2^l2 * 3^l3 * f)^2 
        Alice must choose secret subgroup E_A = E[2^l2] from E as E[2^l2] = <[m_A]P_A + Q_A>, m_A - secret, P_A, Q_A - public
        Bob does mutatis mutandis.
    */
    int l2, l3, f; 
    
    void Generate4LTorsionPoint(EcPoint* P4L);
    void Generate3LTorsionPoint(EcPoint* P3L);
    
    GaloisField* GF;
    PseudoRandomGenerator* prng;
public:
    EllipticCurve* BaseCurve; // think usual form 'd be y^2 = x^3 + x, need to check Edwards equiv.
    // Edwards equiv. would be x^2 + y^2 = 1 -x^2y^2 (d = -1) or isomorphic a=2,d=-2
    EcPoint P_A, Q_A, P_B, Q_B;
    BigInt L3, L2;

    SIDHEngine();
    SIDHEngine(int l2, int l3, int f);
    SIDHEngine(int l2, int l3, int f, BigInt p, int bit_size);
    ~SIDHEngine();

    void GenerateBasePoints();
    void SetBasePoints(const EcPoint* p_A, const EcPoint* q_A, const EcPoint* p_B, const EcPoint* q_B);

    void GeneratePrivateKey(int type, SIDHPrivateKey* priv); // Alice or Bob
    void GeneratePublicKey(const SIDHPrivateKey* priv, SIDHPublicKey* pub );
    void DeriveSharedSecret(const SIDHPrivateKey* priv, const SIDHPublicKey* pub_other, GFElement* sharedSecret);

    /* via Velu formulas, small isogenies of degree 3 and 4 */
    void Compute3Isogeny(const EllipticCurve* E, const EcPointProj* R, GFElement dImg); // E_dImg <- E_d/<R>
    void Compute4Isogeny(const EllipticCurve* E, const EcPointProj* R, GFElement dImg);

    void Evaluate3Isogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, const EcPointProj* P, EcPointProj* PImage);
    void Evaluate4Isogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, const EcPointProj* P, EcPointProj* PImage);

    /*  
        Natural isogeny from E to E/<P>, P belongs to E_A or E_B
        via smooth degree isogeny algorithm(perhaps multiplication strategy), isogenies of degree 3^l, 4^l, phi(E) = F
        Need to check optimal strategies(I suppose they might be vulnarable to timing attacks)
    */
    //void Compute3LIsogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, EllipticCurve* F, int l);
    //void Compute4LIsogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, EllipticCurve* F, int l);


    void ComputeAndEvaluate3LIsogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, const EcPointProj* P, const EcPointProj* Q, EcPointProj* PImage, EcPointProj* QImage, EllipticCurve* F);
    void ComputeAndEvaluate4LIsogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, const EcPointProj* P, const EcPointProj* Q, EcPointProj* PImage, EcPointProj* QImage, EllipticCurve* F);
};

#endif