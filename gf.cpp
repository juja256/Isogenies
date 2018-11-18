#include "gf.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <iomanip>

#define MAX_U64    0xFFFFFFFFFFFFFFFF
#define MSB_M      0x8000000000000000
#define HEX_FORMAT "%.16llX"

static BigInt zero = { 0 };

inline u64 get_bit(const u64* a, u64 num) {
    return a[num/64] & ((u64)1 << (num % 64));
}

inline void copy(u64* a, const u64* b, int len) {
    for (int i=0;i<len;i++)
        a[i] = b[i];
}

inline u64 add(u64 n, const u64* a, const u64* b, u64* c) {
    u64 msb_a, msb_b, carry = 0;

    for (u32 i=0;i < n; i++) {
        msb_a = a[i] & MSB_M;
        msb_b = b[i] & MSB_M;
        c[i] = a[i] + b[i] + carry;
        carry = ( (msb_a && msb_b) || ((msb_a ^ msb_b) && !(MSB_M & c[i])) );
    }
    return carry;
}

inline u64 sub(u64 n, const u64* a, const u64* b, u64* c) {
    u64 borrow = 0;
    for (int i=0; i<n; i++) {
        u64 t_a = a[i];
        u64 t_b = b[i];
        c[i] = t_a - t_b - borrow;
        borrow = ( (~t_a) & (c[i] | t_b) | (c[i] & t_b) ) >> (63);
    }
    return borrow;
}

inline void _mul_raw(u64 a, u64 b, u64* low, u64* high) {
#ifdef _WIN64
    *low = _umul128(a, b, high);
#else
    __int128 r = (__int128)a * (__int128)b;
    *low = (unsigned long long)r;
    *high = r >> 64;
#endif // _WIN64
}

inline u64 _add_raw(u64 a, u64 b, u64* c) {
#ifdef _WIN64
    u64 msb_a, msb_b;
    msb_a = a & MSB_M;
    msb_b = b & MSB_M;
    *c = a + b;
    return ((msb_a && msb_b) || ((msb_a ^ msb_b) && !(MSB_M & *c)));
#else 
    __int128 r = (__int128)a + (__int128)b;
    *c = (u64)r;
    return r >> 64;
#endif 
}

inline void mul_by_word(u64 n, const u64* a, u64 d, u64* c) {
    u64 carry = 0, carry_tmp;
    for (int i=0; i < n; i++) {
        carry_tmp = carry;
        _mul_raw(d, a[i], &(c[i]), &carry);
        carry += _add_raw(c[i], carry_tmp, &(c[i]));
    }
    c[n] = carry;
}

inline void mul(u64 n, const u64* a, const u64* b, u64* c) {
    BigInt tmp;
    memset(c, 0, 2*8*n);
    for (u64 i=0; i < n; i++) {
        mul_by_word(n, a, b[i], (u64*)tmp);
        add(n + 1, c+i, tmp, c+i);
    }
}

/*
    squaring in O(n*logn)
    (a[0] + B a[1] + B^2 a[2] + ... + B^{n-1} a[n-1])^2 = a[0]^2 + B^2 a[1]^2 + ... + B^{2(n-1)} a[n-1]^2 + 2 \Sum{i,j} a[i]a[j] B^{i+j}
    (FF FF)^2 = FE 01 00 00 + 01 FC 02 00 + FE 01 = FF FE 00 01
*/
inline void sqr(u64 n, const u64* a, u64* res) {
    u64 c, c2;
    memset(res, 0, 2*8*n);
    u64 mult [2] = {0, 0};
    u64 carry[2] = {0, 0}; 
    for (int i=0; i<n; i++) {
        _mul_raw(a[i], a[i], &c, &carry[0]);
        carry[0] += _add_raw(res[2*i], c, &res[2*i]);
        carry[1] = 0;
        /* carry[0] - actual 64 bit carry after squaring, carry[1] - still 0 */
        for (int j=i+1; j<n; j++) {
            _mul_raw(a[i], a[j], &mult[0], &mult[1]);

            c2 = _add_raw(mult[0], mult[0], &mult[0]); 
            c = _add_raw(mult[1], mult[1], &mult[1]); 

            mult[1] |= c2; /* mult = 2*a[i]*a[j] */
            mult[1] += _add_raw(mult[0], carry[0], &mult[0]);
            c += _add_raw(mult[1], carry[1], &mult[1]); /* mult = mult + carry */

            mult[1] += _add_raw(mult[0], res[i+j], &res[i+j]); /* res[i+j] += mult[0] */

            carry[0] = mult[1];
            carry[1] = c;
        }
        res[i+n+carry[1]] += _add_raw(res[i+n], carry[0], &res[i+n]) + carry[1]; /* some sort of hack, due to carry[1] might be only 0 or 1 */
    }
}

inline int word_bit_len(u64 n) {
    int c = 64;
    while (c) {
        if ((((u64)1 << (64-1)) & n) >> (64-1))
            return c;
        n <<= 1;
        --c;
    }
    return 0;
}

inline int bigint_bit_len(u64 nWords, const u64* a) {
    int bit_len = nWords * 64;
    int i=nWords-1;
    do {
        bit_len-=64;
    } while ((i>=0) && (a[i--] == 0));

    bit_len += word_bit_len(a[i+1]);
    return bit_len;
}

inline void shl(u64 n, const u64* a, u64* res, u64 bits) {
    u64 buf = 0;
    int chk = bits / 64;
    bits = bits % 64;
    u64 cur;

    if (bits) {
        for (int i = 0; i < n; i++) {
            cur = a[i];
            res[i+chk] = (cur << bits) ^ buf;
            buf = (cur & (MAX_U64 << ( 64-bits ))) >> ( 64-bits );
        }
    }
    else {
        for (int i = 0; i < n; i++) {
            cur = a[i];
            res[i+chk] = cur;
        }
    }

    for (int i = 0; i < chk; i++) {
        res[i] = 0;
    }
    res[n+chk] = buf;
}

inline void shr(u64 n, const u64* a, u64* res, u64 bits) {
    u64 buf = 0;
    int chk = bits / 64;
    bits = bits % 64;
    u64 cur;

    if (bits) {
        for (int i = n-1+chk; i >= chk; i--) {
            cur = a[i];
            res[i-chk] = (cur >> bits) ^ buf;
            buf = (cur & (MAX_U64 >> ( 64-bits ))) << (64 - bits);
        }
    }
    else {
        for (int i = n-1+chk; i >= chk; i--) {
            cur = a[i];
            res[i-chk] = cur;
        }
    }
    for (int i = n; i<n+chk; i++) {
        res[i] = 0;
    }
}

static inline void dump(u64 n, const u64* a) {
    for (int i=n-1; i>=0; i--)
        printf(HEX_FORMAT, a[i]);
    printf("\n");
}

inline int cmp(u64 n, const u64* a, const u64* b) {
    for (int i = n - 1; i >= 0; i--)
    {
        if (a[i] > b[i]) return 1;
        else if (a[i] < b[i]) return -1;
    }
    return 0;
}

void zero_int(u64 n, u64* a) {
    for (int i=0; i<n; i++) a[i] = 0;
}

void add_mod(u64 n, const u64* a, const u64* b, const u64* m, u64* res) {
    res[n] = add(n, a, b, res);
    if (cmp(n+1, res, m) == 1) sub(n+1, res, m, res);
}

void mul_mod(u64 n, const u64* a, const u64* b, const u64* m, u64* res) {
    /* new multiplication with reduction by division */
    BigInt d;
    mul(n, a, b, d);
    divide(n, d, m, NULL, res);
    
    /* old multiplication using only additions */
    /*
    u64 b_len = bigint_bit_len(n, b);
    BigInt mm, r;
    zero_int(n, r);
    
    copy(mm, a, n);
    for (int i=0; i<b_len; i++) {
        if (get_bit(b, i)) add_mod(n, r, mm, m, r);
        add_mod(n, mm, mm, m, mm);
    }
    copy(res, r, n);
    */    
}

void sqr_mod(u64 n, const u64* a, const u64* m, u64* res) {
    BigInt d;
    sqr(n, a, d);
    divide(n, d, m, NULL, res);
}

void exp_mod(u64 n, const u64* a, const u64* p, const u64* m, u64* res) {
    u64 b_len = bigint_bit_len(n, p);
    BigInt mm, r;
    zero_int(n+1, r);
    r[0] = 1;
    copy(mm, a, n);
    for (int i=0; i<b_len; i++) {
        if (get_bit(p, i)) mul_mod(n, r, mm, m, r);
        mul_mod(n, mm, mm, m, mm);
    }
    copy(res, r, n);
}

void imul(u64 n, const u64* a, const u64* b, u64* c) {
    /* c := a * b, where b could be signed value */
    int b_isneg = b[n] & MSB_M;
    BigInt bb;
    if (b_isneg) {
        sub(n+1, zero, b, bb);
    }
    else {
        copy(bb, b, n+1);
    }
    
    mul(n+1, a, bb, c);
    
    if (b_isneg) {
        sub(2*n, zero, c, c);
    }
}

void inv_mod(u64 n, const BigInt a, const BigInt m, BigInt res) {
    /* Old realization of inversion using powering to p-2 */
    /*
    BigInt mm;
    copy(mm, m, n);
    mm[0] -= 2; 
    exp_mod(n, a, mm, m, res);
    */

    /* Compute a^{-1} mod m with Extended Euclidean Algorithm */
    BigInt q, tmp, r;
    BigInt t, newt, newr;
    zero_int(2*n, r);
    zero_int(n+1, t); // t := 0
    zero_int(n+1, newt); newt[0] = 1; // newt := 1
    copy(r, m, n); // r := m
    copy(newr, a, n); // newr := a

    while(cmp(n, newr, zero) != 0) {
        divide(n, r, newr, q, tmp); // q := r div newr, tmp := r mod newr
        copy(r, newr, n); // r := newr
        copy(newr, tmp, n); // newr := tmp 

        imul(n, q, newt, tmp); // tmp := q*newt

        sub(n+1, t, tmp, tmp); // tmp := t - tmp
        copy(t, newt, n+1);
        copy(newt, tmp, n+1);
    }
    
    if (t[n] & MSB_M) {
        add(n, t, m, res);
    }
    else {
        copy(res, t, n);
    }
}

#define div2(n, a) shr((n), (a), (a), 1)

void divide(u64 n, const u64* a, const u64* b, u64* quotient, u64* reminder) {
    BigInt q;
    BigInt tmp;
    BigInt r;

    copy(r, a, 2*n);
    zero_int(2*n, q);
    zero_int(2*n, tmp);

    int k = bigint_bit_len(2*n, a);
    int t = bigint_bit_len(n, b);

    if (k < t) {
        if (quotient) copy(quotient, q, 2*n);
        copy(reminder, r, n);
        return;
    }

    k = k-t;
    shl(n, b, tmp, k);
    while (k >= 0) {
        if (sub(2*n, r, tmp, r) == 0) {
            q[ k/64 ] |= (u64)1 << (k % 64);
        }
        else {
            add(2*n, r, tmp, r);
        }
        
        div2(2*n, tmp);
        k--;
    }

    if (quotient) copy(quotient, q, 2*n);
    copy(reminder, r, n);
}


GaloisFieldException::GaloisFieldException(int c): code(c) {}
std::string GaloisFieldException::What() {
    std::stringstream ss;
    ss << "GaloisFieldException #";
    ss << code;
    return ss.str();
}

GaloisField::GaloisField() {

}

GaloisField::GaloisField(const BigInt characteristic, int ext, int bitSize) {
    this->bitSize = bitSize;
    this->wordSize = (bitSize % ARCH == 0) ? bitSize / ARCH : bitSize / ARCH + 1;
    copy( this->characteristic, characteristic, 2*wordSize);

    if (ext > 2) {
        throw GaloisFieldException(UNSUPPORTED_PARAM);
    }
    this->extension = ext;
    sqr(wordSize, characteristic, this->size);
    shr(2*wordSize, this->size, this->halfSize, 1);
}

const BigInt* GaloisField::GetChar() {
    return &(this->characteristic);
}

int GaloisField::GetWordSize() {
    return wordSize;
}

int GaloisField::GetBitSize() {
    return bitSize;
}

int GaloisField::GetExtension() {
    return extension;
}

const GFElement GaloisField::Zero = {{0}, {0}};
const GFElement GaloisField::Unity = {{1}, {0}};
const GFElement GaloisField::I = {{0}, {1}};

void GaloisField::BaseAdd(const BaseEl a, const BaseEl b, BaseEl c) {
    add_mod(wordSize, a, b, characteristic, c);
}

void GaloisField::BaseSub(const BaseEl a, const BaseEl b, BaseEl c) {
    int borrow  = sub(wordSize, a, b, c);
    if (bitSize % ARCH != 0) {
        borrow = c[wordSize - 1] & ((u64)1 << (wordSize % ARCH) );
    }
    if (borrow != 0) {
        add(wordSize, c, characteristic, c);
    }
}

void GaloisField::Add(const GFElement a, const GFElement b, GFElement c) {
    for (int i=0; i<extension; i++) {
        BaseAdd(a[i], b[i], c[i]);
    }
}

int GaloisField::BaseCmp(const BaseEl a, const BaseEl b) {
    return cmp(wordSize, a, b);
}

#define div2(n, a) shr((n), (a), (a), 1)

// Tonelli-Shanks
void GaloisField::BaseSqrt(const BaseEl a, BaseEl r) {
    // p = Q*2^s
    BaseEl Q, pp, z, tmp, c, t, R, b;
    BaseCopy(z, Unity[0]);
    z[0]++;
    BaseCopy(Q, characteristic);
    u64 M, S = 1;
    div2(wordSize, Q);

    BaseCopy(pp, Q);
    //p-1 = 2^s*Q

    while ((Q[0] & 1) == 0) {
        div2(wordSize, Q);
        S++;
    }
    //finding z - non quadratic residue
    while(1) {
        BasePow(z, pp, tmp); 
        if (BaseCmp(tmp, Unity[0]))
            break;
        else 
            BaseAdd(z, Unity[0], z);
    }
    
    M = S;
    BasePow(z, Q, c);
    BasePow(a, Q, t);
    BaseAdd(Q, Unity[0], tmp);
    div2(wordSize, tmp);
    u64 i=1;
    BasePow(a, tmp, R); // tmp = (Q+1)/2

    while (1) {
        if (!BaseCmp(t, Unity[0])) {
            BaseCopy(r, R);
            break;
        } 
        BaseCopy(tmp, t); //copying value of t, we'll need it later
        for (i=1; i<M; i++) {
            BaseSqr(t, t);
            if (!BaseCmp(t, Unity[0]))
                break;
        }
        BaseCopy(b, c);
        for (u64 j=0; j<M-i-1; j++) {
            BaseSqr(b, b);
        }
        M = i;
        BaseSqr(b, c);
        BaseMul(tmp, c, t); //original value of t multiplied by c^2 
        BaseMul(R, b, R);
    }
}

bool GaloisField::IsQuadraticResidue(const GFElement a) {
    GFElement b;
    Pow(a, halfSize, 2*wordSize, b);
    return BaseCmp(b[0], Unity[0]) == 0;
}

// Algorithm 9 in https://eprint.iacr.org/2012/685.pdf for p = 3 (mod 4)
bool GaloisField::Sqrt(const GFElement a, GFElement r) {
    BigInt p, p2;
    GFElement a1, x0, a0, alpha, b;
    BigInt c;
    BaseCopy(c, characteristic);
    c[0] &= (u64)(-1) << 1;

    shr(wordSize, characteristic, p, 2);
    shr(wordSize, characteristic, p2, 1);

    Pow(a, p, wordSize, a1);
    Mul(a1, a, x0);
    Mul(x0, a1, alpha);
    Pow(alpha, characteristic, wordSize, a0);
    Mul(a0, alpha, a0);
    if (BaseCmp(a0[0], c) == 0) {
        return false;
    }
    if (BaseCmp(alpha[0], c) == 0) { 
        Mul(x0, I, r);
    } else {
        Add(alpha, Unity, b);
        Pow(b, p2, wordSize, b);
        Mul(x0, b, r);
    }
    return true;
} 

void GaloisField::InitFromString(BaseEl a, const char* str) {
    u64 s_len = strlen(str);
    u64 tmp;

    memset(a, 0, sizeof(BaseEl));

    for (int i = s_len-1; i >= 0; i--) {
        if ((str[i] >= '0') && (str[i] <= '9')) { // 0,1,2,...
            tmp = str[i] - 48;
        }
        else if ((str[i] >= 'A') && (str[i] <= 'F')) { // A,B,...
            tmp = str[i] - 55;
        }
        else {
            return;
        }
        a[(s_len - 1 - i) / 16] ^= ( tmp << (((s_len - 1 - i) % 16)*4) );
    }
}

void GaloisField::InitFromString(GFElement a, const char* str) {
    std::string s = str;
    size_t pos = s.find_first_of(";,+");
    
    if (pos != std::string::npos) {
        std::string s1 = s.substr(0, pos);
        std::string s2 = s.substr(pos+1);
        InitFromString(a[0], s1.c_str());
        InitFromString(a[1], s2.c_str());
    } else {
        throw GaloisFieldException(INVALID_DATA);
    }    
}

std::string GaloisField::Dump(const GFElement a) {
    std::stringstream out;
    auto oldFlags = out.flags();
    for (int i=wordSize-1; i>=0; i--) {
        out << std::uppercase << std::setfill('0') << std::setw(16) << std::hex << a[0][i];
    }
    out.flags(oldFlags);
    out << " + ";
    for (int i=wordSize-1; i>=0; i--) {
        out << std::uppercase << std::setfill('0') << std::setw(16) << std::hex << a[1][i];
    }
    out.flags(oldFlags);
    out << "j";
    return out.str();
}

void GaloisField::Sub(const GFElement a, const GFElement b, GFElement c) {
    for (int i=0; i<extension; i++) {
        BaseSub(a[i], b[i], c[i]);
    }
}

void GaloisField::BaseNeg(const BaseEl a, BaseEl c) {
    sub(wordSize, characteristic, a, c);
}

void GaloisField::Neg(const GFElement a, GFElement c) {
    for (int i=0; i<extension; i++) {
        BaseNeg(a[i], c[i]);
    }
}

void GaloisField::BasePow(const BaseEl a, const BigInt n, int nlen, BaseEl b) {
    BaseEl tmp;
    BaseCopy(tmp, a);
    BaseCopy(b, Unity[0]);
    for (u64 i=0; i<nlen; i++) {
        if (get_bit(n, i))
            BaseMul(b, tmp, b);
        BaseSqr(tmp, tmp);
    }
}

void GaloisField::BasePow(const BaseEl a, const BigInt n, BaseEl b) {
    BasePow(a, n, wordSize, b);
}

void GaloisField::Pow(const GFElement a, const BigInt n, int nlen, GFElement b) {
    GFElement tmp;
    Copy(tmp, a);
    Copy(b, Unity);
    for (u64 i=0; i<nlen; i++) {
        if (get_bit(n, i))
            Mul(b, tmp, b);
        Sqr(tmp, tmp);
    }
}

void GaloisField::Pow(const GFElement a, const BigInt n, GFElement b) {
    Pow(a, n, 2*wordSize, b);
}

void GaloisField::BaseInv(const BaseEl a, BaseEl b) {
    // via extended Euclidean
    inv_mod(wordSize, a, characteristic, b);
}

void GaloisField::Inv(const GFElement a, GFElement b) {
    // Thanks to YaSV for reminding basic complex analysis: 
    // 1/(x+yi) = (x-yi)/(x^2 + y^2)
    BaseEl c, d;
    BaseSqr(a[0], c);
    BaseSqr(a[1], d);
    BaseAdd(c, d, c);
    BaseInv(c, d);
    Copy(b, a);
    BaseNeg(b[1], b[1]);
    MulByBase(b, d, b);
}

void GaloisField::BaseMul(const BaseEl a, const BaseEl b, BaseEl c) {
    mul_mod(wordSize, a, b, characteristic, c);
}

void GaloisField::Mul(const GFElement a, const GFElement b, GFElement c) {
    // (x1 + y1*i) * (x2 + y2*i) = x1*x2 - y1*y2 + (y1*x2 + y2*x1)*i
    BaseEl d1, d2;
    BaseMul(a[0], b[0], d1);
    BaseMul(a[0], b[0], d2);

    BaseSub(d1, d2, c[0]);

    BaseMul(a[1], b[0], d1);
    BaseMul(a[0], b[1], d2);

    BaseAdd(d1, d2, c[1]);
}

void GaloisField::BaseSqr(const BaseEl a, BaseEl c) {
    sqr_mod(wordSize, a, characteristic, c);
}

void GaloisField::Sqr(const GFElement a, GFElement c) {
    BaseEl d, b;
    BaseSqr(a[0], d);
    BaseSqr(a[1], b);
    BaseSub(d, b, c[0]);
    BaseMul(a[0], a[1], c[1]);
    BaseAdd(c[1], c[1], c[1]);
}

void GaloisField::BaseCopy(BaseEl a, const BaseEl b) {
    copy(a, b, wordSize);
}

void GaloisField::Copy(GFElement a, const GFElement b) {
    for (int i=0; i<extension; i++) {
        BaseCopy(a[i], b[i]);
    }
}

void GaloisField::MulByBase(const GFElement a, const BaseEl e, GFElement c) {
    for (int i=0; i<extension; i++) {
        BaseMul(a[i], e, c[i]);
    }
}

void GaloisField::MulBy2Power(const GFElement a, int pp, GFElement b) {
    for (int i=0; i<extension; i++) {
        BigInt d;
        BaseCopy(d, Zero[0]);
        shl(wordSize, a[i], d, pp);
        divide(wordSize, d, characteristic, NULL, b[i]);
    }
}

bool GaloisField::Equal(const GFElement a, const GFElement b) {
    return (BaseCmp(a[0], b[0]) == 0) && (BaseCmp(a[1], b[1]) == 0);
}