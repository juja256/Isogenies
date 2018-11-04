#ifndef ISOGENY_H
#define ISOGENY_H

#include "ec.h"

class IsogenyEngine {
    /* 
        Parameters identifying group structure of E over field GF(p^2), p = (2^l2 * 3^l3 * f)-1 prime; 
        E - supersingular |E| = (2^l2 * 3^l3 * f)^2 
        Alice must choose secret subgroup E_A = E[2^l2] from E as E[2^l2] = <[m_A]P_A + Q_A>, m_A - secret, P_A, Q_A - public
        Bob does mutatis mutandis.
    */
    int l2, l3, f; 
    EllipticCurve* BaseCurve; // think usual form 'd be y^2 = x^3 + x, need to check Edwards equiv.

public:

    IsogenyEngine(int l2, int l3, int f);
    ~IsogenyEngine();
    EllipticCurve* GetBaseCurve();

    /* via Velu formulas, small isogenies of degree 3 and 4 */
    void Compute3Isogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, EllipticCurve* F, EcPointProj* image);
    void Compute4Isogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, EllipticCurve* F, EcPointProj* image);

    /*  
        Natural isogeny from E to E/<P>, P belongs to E_A or E_B
        via smooth degree isogeny algorithm(perhaps multiplication strategy), isogenies of degree 3^l, 4^l, phi(E) = F
        Need to check optimal strategies(I suppose they might be vulnarable to timing attacks)
    */
    void Compute3LIsogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, int l, EllipticCurve* F, EcPointProj* image);
    void Compute4LIsogeny(const EllipticCurve* E, const EcPointProj* kernelPoint, int l, EllipticCurve* F, EcPointProj* image);
};



#endif