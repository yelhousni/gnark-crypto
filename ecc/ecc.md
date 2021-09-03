## Supported curves

* BLS12-381 (Zcash)
* BN254 (Ethereum)
* BLS12-377 (ZEXE)
* BW6-761 (EC supporting pairing on BLS12-377 field of definition)
* BLS12-379 (alternative to BLS12-377 with higher 2-adicity)
* BW6-672 (EC supporting pairing on BLS12-379 field of definition)
* BW6-672 (EC supporting pairing on BLS12-379 field of definition with conservative security)
* BLS24-315 (optmized for KZG-based SNARKs)
* BW6-633 (EC supporting pairing on BLS24-315 field of definition)

### Twisted edwards curves

Each of these curve has a `twistededwards` sub-package with its companion curve. In particular, BLS12-381 comapnion curve is known as [Jubjub](https://z.cash/technology/jubjub/) and BN254's [Baby-Jubjub](https://iden3-docs.readthedocs.io/en/latest/_downloads/33717d75ab84e11313cc0d8a090b636f/Baby-Jubjub.pdf).

They are of particular interest as they allow efficient elliptic curve cryptography inside zkSNARK circuits.
