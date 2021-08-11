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

package bw6767

import (
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/consensys/gnark-crypto/ecc/bw6-767/fp"
	"github.com/consensys/gnark-crypto/ecc/bw6-767/fr"
)

// E: y**2=x**3+1
// Etwist: y**2 = x**3+3
// Tower: Fp->Fp6, u**6=3
// Generator (same as BLS381): x=-15132376222941642752
// optimal Ate loops: x+1, x**2-x-1
// Fp: p=496597749679620867773432037469214230242402307330180853437434581099336634619713640485778675608223760166307530047354464605410050411581079376994803852937842168733702867087556948851016246640584660942486895230518034810309227309966899431
// Fr: r=4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787

// ID BW6_767 ID
const ID = ecc.BW6_767

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
var loopCounter1 [65]int8
var loopCounter2 [192]int8

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
	bTwistCurveCoeff.SetUint64(3)

	twoInv.SetOne().Double(&twoInv).Inverse(&twoInv)

	// E(1,y)*cofactor
	g1Gen.X.SetString("127687253511432941835499154999732953539969793860764514205013635996439242747457934431893570832266740963864950713809357287070846939000367049554519743864924323440810949629217677483481194663331926309250818003412838087592587472550707218")
	g1Gen.Y.SetString("415570529523170147223250223671601071129165798689804006717876771297003017718159840368703823786319144396618898691682149260290217115399107531975419658973137909698922937988511368601419289861827304905241655385035120916874417442125721204")
	g1Gen.Z.SetString("1")

	// Etwist(1,2)*cofactor
	g2Gen.X.SetString("370611171465172359348863648443534520144617072349884185652206813771489664034831143983178049920510836078361116088420840622225267322852644540540617123958979924966938307707664543525950567252218300954395355151658118858470703533448342222")
	g2Gen.Y.SetString("455144308204607096185992716699045373884508292978508084510087807751472279103896568109582325400258900176330927780121791269969939391813736974371796892558810828460226121428602798229282770695472612961143258458821149661074127679136388603")
	g2Gen.Z.SetString("1")

	g1GenAff.FromJacobian(&g1Gen)
	g2GenAff.FromJacobian(&g2Gen)

	//binary decomposition of -15132376222941642751, little endian
	// xGen+1
	T1, _ := new(big.Int).SetString("15132376222941642751", 10) // negative
	ecc.NafDecomposition(T1, loopCounter1[:])
	// xGen^3-xGen^2-xGen
	T2, _ := new(big.Int).SetString("3465144826073652319005258340840392356319973669502814453760", 10) // negative
	ecc.NafDecomposition(T2, loopCounter2[:])

	g1Infinity.X.SetOne()
	g1Infinity.Y.SetOne()
	g2Infinity.X.SetOne()
	g2Infinity.Y.SetOne()

	thirdRootOneG1.SetString("451452499708746243421442696394275804592767119751118962106882058158528025766103643615697202253207413006991058800455542766924935899310685166148099708594514571753800103096705086912881023032622324847956780035251378028187894066092550170")
	thirdRootOneG2.Square(&thirdRootOneG1)
	lambdaGLV.SetString("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436", 10) // (x**5-3*x**4+3*x**3-x+1)
	_r := fr.Modulus()
	ecc.PrecomputeLattice(_r, &lambdaGLV, &glvBasis)

	xGen.SetString("15132376222941642752", 10)

}

// Generators return the generators of the r-torsion group, resp. in ker(pi-id), ker(Tr)
func Generators() (g1Jac G1Jac, g2Jac G2Jac, g1Aff G1Affine, g2Aff G2Affine) {
	g1Aff = g1GenAff
	g2Aff = g2GenAff
	g1Jac = g1Gen
	g2Jac = g2Gen
	return
}
