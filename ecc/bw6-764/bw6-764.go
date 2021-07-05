// Copyright 2020 ConsenSys AG
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package bw6764

import (
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/consensys/gnark-crypto/ecc/bw6-764/fp"
	"github.com/consensys/gnark-crypto/ecc/bw6-764/fr"
)

// https://eprint.iacr.org/2020/351.pdf

// E: y**2=x**3+1
// Etwist: y**2 = x**3+2
// Tower: Fp->Fp6, u**6=2
// Generator (same as BLS379): x=11170052975785672705
// optimal Ate loops: x+1, x**2-x-1
// Fp: p=0xb445f0691766d95a1ba736fcc2d419ddb11f4639db77f4fec6ee87d3b4d1280a904786c4f97be0e2f12a8f4577b62835d84ff875d577945e4d8b84d76f9ef9ec47c22f24fbf3b3e63c2f132d64000585dd86710000000abf870000000000085
// Fr: r=0x434e417092cfe75084c660572fb575ee95ce40bf8a2570025e48630b58000001f49f2b0000000009b04000000000001

// ID BW6_764 ID
const ID = ecc.BW6_764

// bCurveCoeff b coeff of the curve
var bCurveCoeff fp.Element

// bTwistCurveCoeff b coeff of the twist (defined over Fp) curve
var bTwistCurveCoeff fp.Element

// twoInv 1/2 mod p (needed for DoubleStep in Miller loop)
var twoInv fp.Element

// generators of the r-torsion group, resp. in ker(pi-id), ker(Tr)
var g1Gen G1Jac
var g2Gen G2Jac

var g1GenAff G1Affine
var g2GenAff G2Affine

// point at infinity
var g1Infinity G1Jac
var g2Infinity G2Jac

// optimal Ate loop counters
// Miller loop 1: f(P), div(f) = (x+1)(Q)-([x+1]Q)-x(O)
// Miller loop 2: f(P), div(f) = (x**3-x**2-x)(Q) -([x**3-x**2-x]Q)-(x**3-x**2-x-1)(O)
var loopCounter1 [64]int8
var loopCounter2 [191]int8

// Parameters useful for the GLV scalar multiplication. The third roots define the
//  endomorphisms phi1 and phi2 for <G1Affine> and <G2Affine>. lambda is such that <r, phi-lambda> lies above
// <r> in the ring Z[phi]. More concretely it's the associated eigenvalue
// of phi1 (resp phi2) restricted to <G1Affine> (resp <G2Affine>)
// cf https://www.cosic.esat.kuleuven.be/nessie/reports/phase2/GLV.pdf
var thirdRootOneG1 fp.Element
var thirdRootOneG2 fp.Element
var lambdaGLV big.Int

// glvBasis stores R-linearly independant vectors (a,b), (c,d)
// in ker((u,v)->u+vlambda[r]), and their determinant
var glvBasis ecc.Lattice

// generator of the curve
var xGen big.Int

func init() {

	bCurveCoeff.SetOne()
	bTwistCurveCoeff.SetUint64(2)

	twoInv.SetOne().Double(&twoInv).Inverse(&twoInv)

	g1Gen.X.SetString("44031393204798314061162301011099048410666125502551608273807339655348007258308439555989129477416050025082158725561903091504130709478789292781757614002787101455791814151682932718971253146149898579793180882066624371549243207336115417")
	g1Gen.Y.SetString("48501644543217168960833655345577303547571873451429994737312152557265291217152891120792008214005929577974077078556592160946426634130594721147752394732659881102518901729284240600820967888508198036478650854754237946510651144879101841")
	g1Gen.Z.SetString("1")

	g2Gen.X.SetString("47681727899945315873713685533093641276713506398361651846416031161365545349674593454585382891255597279216567457830083595832686502063956383419129647840479540512677250396440012527207006600158651592613071368553421400948346660108850388")
	g2Gen.Y.SetString("50448368292571048916194308036311359878125593899436824920346893655914800862589673667005835100405280622843727689201669776611292055612381731192826342137074796676758316455104383817907603182334857958413856510528811358426958309253880221")
	g2Gen.Z.SetString("1")

	g1GenAff.FromJacobian(&g1Gen)
	g2GenAff.FromJacobian(&g2Gen)

	//binary decomposition of 11170052975785672706, little endian
	// xGen+1
	loopCounter1 = [64]int8{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1}

	// xGen^3-xGen^2-xGen
	T, _ := new(big.Int).SetString("1393688442285558804941521504723062763365228464444220112895", 10)
	ecc.NafDecomposition(T, loopCounter2[:])

	g1Infinity.X.SetOne()
	g1Infinity.Y.SetOne()
	g2Infinity.X.SetOne()
	g2Infinity.Y.SetOne()

	thirdRootOneG1.SetString("12058137234552858436467312031288370469153666779562229030530168593599403961809626412081945255294308728010875871177842428330178991050270518836672998307525535915184491141469668374776448025924800720761988128920570186787148628915049193") // (-48-127*x-115*x^2+162*x^3+109*x^4-500*x^5+1119*x^6-481*x^7-1168*x^8+1668*x^9-848*x^10+163*x^11)/51
	thirdRootOneG2.Square(&thirdRootOneG1)
	lambdaGLV.SetString("173890623291670311836512360316100953774879339360001038672113943180607164901767459447515297873921", 10) // (1-x+3*x^3-3*x^4+x^5)
	_r := fr.Modulus()
	ecc.PrecomputeLattice(_r, &lambdaGLV, &glvBasis)

	xGen.SetString("11170052975785672705", 10)

}

// Generators return the generators of the r-torsion group, resp. in ker(pi-id), ker(Tr)
func Generators() (g1Jac G1Jac, g2Jac G2Jac, g1Aff G1Affine, g2Aff G2Affine) {
	g1Aff = g1GenAff
	g2Aff = g2GenAff
	g1Jac = g1Gen
	g2Jac = g2Gen
	return
}
