#ifndef GF_H
#define GF_H

#define BASE_FIELD_LEN (751/64 + 1 + 1) // 13 64bit words
#define FIELD_EXTENSION 2

#define ARCH 64
#define BYTES_IN_WORD (ARCH/8)

#define UNSUPPORTED_PARAM -1
#define INVALID_DATA -2


#include <string>

typedef unsigned long long u64;
typedef unsigned u32;
typedef unsigned char u8;

typedef u64 BaseEl[BASE_FIELD_LEN];
typedef u64 BigInt[2*BASE_FIELD_LEN];

typedef BaseEl GFElement[FIELD_EXTENSION];

class GaloisFieldException {
    int code;
public:
    GaloisFieldException(int c);
    std::string What();
};

class GaloisField {
    int bitSize;
    int wordSize;
    BigInt characteristic;
    int extension;
    
    BigInt size;
    BigInt halfSize;
public:
   

    GaloisField();
    GaloisField(const BigInt characteristic, int ext, int bitSize);
    const BigInt* GetChar();
    int GetWordSize();
    int GetBitSize();
    int GetExtension();

    static const GFElement Zero;
    static const GFElement Unity;
    static const GFElement I;

    void GFInitFromString(BaseEl a, const char* str);
    void GFInitFromString(GFElement a, const char* str);
    std::string GFDump(const GFElement a);

    bool IsQuadraticResidue(const GFElement a);
    void GFCopy(GFElement a, const GFElement b);
    bool GFCmp(const GFElement a, const GFElement b);
    void GFAdd(const GFElement a, const GFElement b, GFElement c);
    void GFSub(const GFElement a, const GFElement b, GFElement c);
    void GFNeg(const GFElement a, GFElement c);
    void GFPow(const GFElement a, const BigInt n, int nlen, GFElement b);
    void GFPow(const GFElement a, const BigInt n, GFElement b);
    void GFInv(const GFElement a, GFElement b);
    void GFMul(const GFElement a, const GFElement b, GFElement c);
    void GFSqr(const GFElement a, GFElement c);
    bool GFSqrt(const GFElement a, GFElement r); 
    void GFMulBy2Power(const GFElement a, int pp, GFElement b);
    void GFMulByBase(const GFElement a, const BaseEl e, GFElement c);

    /* base field operations */
    void GFBaseCopy(BaseEl a, const BaseEl b);
    void GFBaseSqr(const BaseEl a, BaseEl c);
    void GFBaseMul(const BaseEl a, const BaseEl b, BaseEl c);
    void GFBaseAdd(const BaseEl a, const BaseEl b, BaseEl c);
    void GFBaseSub(const BaseEl a, const BaseEl b, BaseEl c);
    void GFBasePow(const BaseEl a, const BigInt n, int nlen, BaseEl b);
    void GFBasePow(const BaseEl a, const BigInt n, BaseEl b);
    void GFBaseInv(const BaseEl a, BaseEl b);
    void GFBaseNeg(const BaseEl a, BaseEl c);
    void GFBaseSqrt(const BaseEl a, BaseEl b);
    int GFBaseCmp(const BaseEl a, const BaseEl b);
    
};

/* Common Arithmetics */

void shr(u64 n, const u64* a, u64* res, u64 bits);
void shl(u64 n, const u64* a, u64* res, u64 bits);
u64 get_bit(const u64* a, u64 num);
void copy(u64* a, const u64* b, int len);
u64 add(u64 n, const u64* a, const u64* b, u64* c);
u64 sub(u64 n, const u64* a, const u64* b, u64* c);
void _mul_raw(u64 a, u64 b, u64* low, u64* high);
u64 _add_raw(u64 a, u64 b, u64* c);
void mul_by_word(u64 n, const u64* a, u64 d, u64* c);
void mul(u64 n, const u64* a, const u64* b, u64* c);
void imul(u64 n, const u64* a, const u64* b, u64* c);
void sqr(u64 n, const u64* a, u64* res);
int word_bit_len(u64 w);
int bigint_bit_len(u64 n, const u64* a);
u64 sub_word(u64 n, const u64* a, const u64 b, u64* c);
u64 add_word(u64 n, const u64* a, const u64 b, u64* c);

void divide(u64 n, const u64* a, const u64* b, u64* quotient, u64* reminder);
int cmp(u64 n, const u64* a, const u64* b);


void add_mod(u64 n, const u64* a, const u64* b, const u64* m, u64* res);
void mul_mod(u64 n, const u64* a, const u64* b, const u64* m, u64* res);
void exp_mod(u64 n, const u64* a, const u64* b, const u64* m, u64* res);
//void inv_mod(u64 n, const u64* a, const u64* m, u64* res); 


#endif /* GF_H */
