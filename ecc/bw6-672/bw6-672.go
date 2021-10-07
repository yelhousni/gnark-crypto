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

package bw6672

import (
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/consensys/gnark-crypto/ecc/bw6-672/fp"
	"github.com/consensys/gnark-crypto/ecc/bw6-672/fr"
)

// E: y**2=x**3-4
// Etwist: y**2 = x**3-4/3
// Tower: Fp->Fp6, u**6=3
// Generator (same as BLS24-315): x=-3218079743
// Fp: p=10298672283669831003294162666372200112722413308864335408449041451781560030635870917520584453632890853794509234527426523862901218215709508411167619988692340017052754441213777163725504913151409615827435527
// Fr: r=39705142709513438335025689890408969744933502416914749335064285505637884093126342347073617133569

// ID BW6_672 ID
const ID = ecc.BW6_672

// bCurveCoeff b coeff of the curve
var bCurveCoeff fp.Element

// twist
var twist fp.Element

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
var loopCounterAte0 [33]int8
var loopCounterAte1 [159]int8
var loopCounterOptTate0 [158]int8
var loopCounterOptTate1 [158]int8

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

	bCurveCoeff.SetUint64(4).Neg(&bCurveCoeff)
	twist.SetUint64(3)
	bTwistCurveCoeff.Inverse(&twist).Mul(&bTwistCurveCoeff, &bCurveCoeff) // D-twist

	twoInv.SetOne().Double(&twoInv).Inverse(&twoInv)

	// E1(1,y)*cofactor
	g1Gen.X.SetString("8839124564855791793291047299347745152935735884399481072595505567430987344292088438819284921669192377581484035322844603120596848898640559427651237481367056532887348973571742747864903336241016425676648922")
	g1Gen.Y.SetString("7030649000584214638488653735360672247062444076266380763893123242676220279584815192649273277833411485723399215055414665529503749372485469042045766616632147119532832921197569331464959943179273541348355402")
	g1Gen.Z.SetString("1")

	// E2(1,y))*cofactor
	g2Gen.X.SetString("6864525156298809514087621905098536023812176437678895007490281362985910323101181274535591103134655916073644759268162288059811105140364213622145788959502942353195723241371518155372106541540935111983324661")
	g2Gen.Y.SetString("4310630220448298840815264308633992670665637623470083878928169211767078121165763030830443894380267956933199856451252613613846744817142124346847488950342334298598067007454754698671283277762399410956207422")
	g2Gen.Z.SetString("1")

	g1GenAff.FromJacobian(&g1Gen)
	g2GenAff.FromJacobian(&g2Gen)

	//binary decomposition of 3218079742, little endian
	// xGen+1 (negative)
	T, _ := new(big.Int).SetString("3218079742", 10)
	ecc.NafDecomposition(T, loopCounterAte0[:])

	// xGen^5-xGen^4-xGen (negative)
	T, _ = new(big.Int).SetString("345131030376204096837580131803633448876874137601", 10)
	ecc.NafDecomposition(T, loopCounterAte1[:])

	// binary decomposition of u+1 (negative)
	loopCounterOptTate0 = [158]int8{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
	// binary decomposition of u^5-u^4-u (negative)
	loopCounterOptTate1 = [158]int8{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1}

	g1Infinity.X.SetOne()
	g1Infinity.Y.SetOne()
	g2Infinity.X.SetOne()
	g2Infinity.Y.SetOne()

	thirdRootOneG1.SetString("3018874115714335214819350859895104554342871833755237358128964215684826616475473951096168379451772131547146048038621211116857857462584990239008506640664973892715664806650352079858457256024893835939479141")
	thirdRootOneG2.Square(&thirdRootOneG1)
	lambdaGLV.SetString("39705142672498995661671850106945620852186608752525090699191017895721506694646055668218723303426", 10) // 1-x+2*x^2-2*x^3+3*x^5-4*x^6+4*x^7-3*x^8+x^9
	_r := fr.Modulus()
	ecc.PrecomputeLattice(_r, &lambdaGLV, &glvBasis)

	xGen.SetString("3218079743", 10) // negative

}

// Generators return the generators of the r-torsion group, resp. in ker(pi-id), ker(Tr)
func Generators() (g1Jac G1Jac, g2Jac G2Jac, g1Aff G1Affine, g2Aff G2Affine) {
	g1Aff = g1GenAff
	g2Aff = g2GenAff
	g1Jac = g1Gen
	g2Jac = g2Gen
	return
}
