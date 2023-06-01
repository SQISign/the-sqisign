#include <protocols.h>
#include <string.h>

int sqisign_keypair(unsigned char *pk, unsigned char *sk) { 
    int ret;
    secret_key_t skt;
    public_key_t pkt = { 0 };
    secret_key_init(&skt);

    ret = !protocols_keygen(&pkt, &skt);

    secret_key_encode(sk, &skt, &pkt);
    public_key_encode(pk, &pkt);
    secret_key_finalize(&skt);
    return ret;
}

int sqisign_sign(unsigned char *sm,
              unsigned long long *smlen, const unsigned char *m,
              unsigned long long mlen, const unsigned char *sk) {
    int ret = 0;
    secret_key_t skt;
    public_key_t pkt = { 0 };
    signature_t sigt;
    secret_key_init(&skt);
    signature_init(&sigt);
    secret_key_decode(&skt, &pkt, sk);

    ret = !protocols_sign(&sigt, &pkt, &skt, m, mlen);
    signature_encode(sm, &sigt);

    memcpy(sm + SIGNATURE_LEN, m, mlen);
    *smlen = SIGNATURE_LEN + mlen;

    secret_key_finalize(&skt);
    signature_finalize(&sigt);
    return ret;
}

int sqisign_open(unsigned char *m,
              unsigned long long *mlen, const unsigned char *sm,
              unsigned long long smlen, const unsigned char *pk) { 
    int ret = 0;
    public_key_t pkt = { 0 };
    signature_t sigt;
    signature_init(&sigt);

    public_key_decode(&pkt, pk);
    signature_decode(&sigt, sm);

    ret = !protocols_verif(&sigt, &pkt, sm + SIGNATURE_LEN, smlen - SIGNATURE_LEN);

    if (!ret) {
        *mlen = smlen - SIGNATURE_LEN;
        memmove(m, sm + SIGNATURE_LEN, *mlen);
    }

    signature_finalize(&sigt);
    return ret;
}

int sqisign_verify(const unsigned char *m,
                unsigned long long mlen, const unsigned char *sig,
                unsigned long long siglen, const unsigned char *pk) {

    int ret = 0;
    public_key_t pkt = { 0 };
    signature_t sigt;
    signature_init(&sigt);

    public_key_decode(&pkt, pk);
    signature_decode(&sigt, sig);

    ret = !protocols_verif(&sigt, &pkt, m, mlen);

    signature_finalize(&sigt);
    return ret;
}

