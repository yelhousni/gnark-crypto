<a name="v0.5.0"></a>
## [v0.5.0] - 2021-08-20

### Breaking changes
- twisted Edwards BN-companion in reduced form (a=-1): this affect `eddsa`. `v0.4.0` and `v0.5.0` keys and signatures are not compatible. 

### Feat
- adds new curve bls24-315
- adds new curve bw6-633
- adds kzg polynomial commitment scheme
- adds fiat shamir
- Element.SetInterface returns an error instead of panicking if unsupported type
- MultiExp now takes a nbTasks parameter and splits until we have nbTasks <= nbChunks
- MultiExp returns error if len(points) != len(scalars)
- ecc encoder now handles []Element so gnark don't have to
- ecc encoders uses binary.Write and binary.Read to support basic types
- added ecc.Implemented() that returns list of curve fully implemented
- added Reference bencharks for continuous benchmarking. fixes [#54](https://github.com/ConsenSys/gnark-crypto/issues/54)
- added curve level go-fuzz fuzz functions
- **all curves:** faster GT memebership
- **twisted Edwards:** tests use gopter, no more hardcoded values
- **bls12-377:** change G2 generator (+Fp QNR) to match other libs
- **bls12-377:** change G1 generator to match other libs
- **bw6:** Pairing according to ABLR 2013/722 with Fp6/Fp3

### Fix
- use crypto/rand instead of math/rand in ecc/../utils.go
- fixes [#51](https://github.com/ConsenSys/gnark-crypto/issues/51)
- e2 x86 asm incorrect offset when x is 0
- fixes [#49](https://github.com/ConsenSys/gnark-crypto/issues/49)
- **twisted Edwards:** fixed Neg(), and fixes [#57](https://github.com/ConsenSys/gnark-crypto/issues/57)

### Perf
- **all curves:** twisted Edwards companions arithmetic with a=-1
- **bls12:** faster G2 clear cofactor
- **bls12:** faster G2 subgroup checks --> psi^2=phi+1
- **bls12:** faster G2 subgroup checks
- **bls12-377:** remove one add, one sub in e2.Square
- **bn:** optimize Expt (no conditional branching)
- **bn254:** Expt in 2-NAF
- **bw6:** replace Inverse and FrobeniusCube by conjugate
- **bw6:** new optimized final exp (hard part)
- **bw6-633:** divide G1 cofactor formula by 4
- **bw6-633:** optimized hard part in final exp
- **fft:** introduced flatten kernel for n==8 and asm impl for butterfly to minimze memory writes

### Refactor
- ported accumulator/ and polynomial/ from gnark
- moved fr/polynomial/kzg to fr/kzg
- removed deprecated MulAssign, AddAssign and SubAssign apis
- removed hash functions recorded in transcript.go
- moved crypto/* under /
- **kzg:** Proof -> OpeningProof. BatchProofsSinglePoint -> BatchOpeningProof
- **kzg:** removed Scheme, package level methods with SRS and domain as parameter

### Test
- added mulGeneric vs mul assembly on E2
- **curves:** use IsInSubGroup instead IsOnCurve MapToCurveG1Svdw test
- added e2.Neg test in code generation





<a name="v0.4.0"></a>
## [v0.4.0] - 2021-03-31

### Refactor
- gurvy -> gnark-crypto
- moved interop tests under github.com/consensys/gnark-tests
- bls381 -> bls12-381
- bls377 -> bls12-377
- bn256 -> bn254
- migrated MiMC and EdDSA from gnark into gnark-crypto
- migrated gnark/backend/fft into gnark-crypto
- migrated goff packages into ./field/...
- cleaning internal/generator pattern

### Ci
- testing with go 1.15, go 1.16 on Windows, MacOS, Linux (+arch=32bits)

### Docs
- added ecc/ecc.md and field/field.md

### Feat
- multiExp in full extended jacobian coordinates

### Fix
- handle case where numCPU < 4 in precomputeExpTable
- incorrect comment and size returned in twistededwards SetBytes fixes [#34](https://github.com/ConsenSys/gnark-crypto/issues/34)
- point.SetBytes can now be called concurently with same byte slice input



<a name="v0.3.8"></a>
## [v0.3.8] - 2021-02-01

### Bls377
- final exp hard part eprint 2020/875
- ML entirely on the twist (ABLR)

### Bls381
- final exp hard part eprint 2020/875
- ML entirely on the twist (ABLR)
- change G1 and G2 generators for interop

### Bn256
- inline lineEval() in MilleLoop
- ML entirely on the twist (ABLR)
- change G1 and G2 generators for interop

### Bw6
- add E6 and pairing tests
- correct comments in FinalExp
- fix bw6 pairing API to take slices of points and mutualize squares
- change G1 and G2 generators for interop

### Pull Requests
- Merge pull request [#29](https://github.com/ConsenSys/gnark-crypto/issues/29) from ConsenSys/youssef/bls12-finalExp
- Merge pull request [#27](https://github.com/ConsenSys/gnark-crypto/issues/27) from ConsenSys/experimental/pairing
- Merge pull request [#26](https://github.com/ConsenSys/gnark-crypto/issues/26) from ConsenSys/youssef/ML-ABLR
- Merge pull request [#25](https://github.com/ConsenSys/gnark-crypto/issues/25) from ConsenSys/csquare
- Merge pull request [#23](https://github.com/ConsenSys/gnark-crypto/issues/23) from ConsenSys/youssef/bw6-API-pairing

