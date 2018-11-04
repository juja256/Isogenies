#ifndef GF_H
#define GF_H

#define BASE_FIELD_LEN (751/64 + 1 + 1) // 13 64bit words
#define FIELD_EXTENSION 2

typedef unsigned long long u64;
typedef unsigned u32;
typedef unsigned char u8;

typedef u64 BaseEl[BASE_FIELD_LEN];
typedef u64 BigInt[2*BASE_FIELD_LEN];



typedef struct {
    BaseEl a[FIELD_EXTENSION];
} GFElement;

class GaloisField {
    BigInt characteristic;
    int bitSize;
    int extension;

public:
    static GFElement Zero;
    static GFElement Unity;

    void GFAdd(const GFElement, const GFElement, GFElement);
    void GFSqrt(const GFElement a, GFElement r); // via Tonelli-Shanks
    void GFInitFromString(GFElement a, const char* str);
    void GFDump(const GFElement a);
    void GFAdd(const GFElement a, const GFElement b, GFElement c);
    void GFSub(const GFElement a, const GFElement b, GFElement c);
    void GFNeg(const GFElement a, GFElement c);
    void GFPow(const GFElement a, const BigInt n, GFElement b);
    void GFInv(const GFElement a, GFElement b);
    int  GFCmp(const GFElement a, const GFElement b);
    void GFMul(const GFElement a, const GFElement b, GFElement c);
    void GFSqr(const GFElement a, GFElement c);
    void GFMulBy2Power(const GFElement a, int pp, GFElement b);
};

/* Common Arithmetics */

void shr(u64 n, const u64* a, u64* res, u64 bits);
#define div2(n, a) shr((n), (a), (a), 1)
void shl(u64 n, const u64* a, u64* res, u64 bits);
#define mul2(n, a) shl((n), (a), (a), 1)
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


void add_mod(u64 n, const BigInt a, const BigInt b, const BigInt m, BigInt res);
void mul_mod(u64 n, const BigInt a, const BigInt b, const BigInt m, BigInt res);
void exp_mod(u64 n, const BigInt a, const BigInt b, const BigInt m, BigInt res);
void inv_mod(u64 n, const BigInt a, const BigInt m, BigInt res);


#endif /* GF_H */
