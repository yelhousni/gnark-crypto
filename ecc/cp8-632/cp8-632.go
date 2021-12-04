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

package cp8632

import (
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/consensys/gnark-crypto/ecc/cp8-632/fp"
	"github.com/consensys/gnark-crypto/ecc/cp8-632/internal/fptower"
)

// E: y**2=x**3-x
// Etwist: y**2=x**3-x/u
// Tower: Fp->Fp8 (w^2=2)
// Fp: p=16857842227199225999646786835636637546980511455593181936643515867446449681250963223102065881329922732942164707714798074010252040453159496268796306973476256192757039177873952554508334906059909 (632 bits)
// Fr: r=39705142709513438335025689890408969744933502416914749335064285505637884093126342347073617133569 (= Fp BLS12-315, 315 bits)

// ID CP8_632 ID
const ID = ecc.CP8_632

// aCurveCoeff a coeff of the curve
var aCurveCoeff fp.Element

// twist
var twist fptower.E2

// cofactor of the curve
var cofactor big.Int

// aTwistCurveCoeff b coeff of the twist (defined over Fp2) curve
var aTwistCurveCoeff fptower.E2

// cofactorTwist of the curve
var cofactorTwist big.Int

// twoInv 1/2 mod p (needed for DoubleStep in Miller loop)
var twoInv fp.Element

// generators of the r-torsion group, resp. in ker(pi-id), ker(Tr)
var g1Gen G1W12
var g2Gen G2Jac

var g1GenAff G1Affine
var g2GenAff G2Affine

// point at infinity
var g1Infinity G1W12
var g2Infinity G2Jac

// Miller loop counters
var loopCounter [33]int8 // ???

/*
// Parameters useful for the GLV scalar multiplication. The third roots define the
//  endomorphisms phi1 and phi2 for <G1Affine> and <G2Affine>. lambda is such that <r, phi-lambda> lies above
// <r> in the ring Z[phi]. More concretely it's the associated eigenvalue
// of phi1 (resp phi2) restricted to <G1Affine> (resp <G2Affine>)
// cf https://www.cosic.esat.kuleuven.be/nessie/reports/phase2/GLV.pdf
var thirdRootOneG1 fp.Element
var thirdRootOneG2 fp.Element
var lambdaGLV big.Int
*/

// glvBasis stores R-linearly independant vectors (a,b), (c,d)
// in ker((u,v)->u+vlambda[r]), and their determinant
var glvBasis ecc.Lattice

// psi o pi o psi**-1, where psi:E->E' is the degree 6 iso defined over Fp8
var endo struct {
	u fptower.E2
	v fptower.E2
}

func init() {

	aCurveCoeff.SetUint64(1).Neg(&aCurveCoeff)
	twist.A1.SetUint64(1)
	aTwistCurveCoeff.Inverse(&twist).Mul(&aTwistCurveCoeff, &twist) // D-twist

	twoInv.SetUint64(2).Inverse(&twoInv)

    cofactor.SetString("424575787336486525919881433149779109065117535008715395416548814899802759612501346622389123445160", 10)
    cofactorTwist.SetString("7157431636407380950126640349841351740878278991117729287857214041004658272402879684508872603937823124718589641314727702909291567620337487154057402119860855486765390772491717966442838511543044147992594893712246699516873984803821277362429144406100389878217469407870650419787619681584093298", 10)

	// random, not deterministic
	g1Gen.X.SetString("16669241569793144689576358184015273697172473599133409527765597244456059063645632082294888602687297549313408226230020161845797893266328116331226196965713808831871989472075519327658967992069490")
	g1Gen.Y.SetString("1248052385589995465912547602370558205558212431448411201256289955347018540821720119347513967254758856831331395724637858977285509974129230385564546422828003505266952129796380389964953474092787")
	g1Gen.Z.SetString("1")

	// random, not deterministic
    g2Gen.X.SetString("14354547420196287803704182366476680632365830220216070176151071330154402298831834623870632909972194516703719909648086691714742669516707712338402515313317968011902252305693941493080807819072390",
		"2395863994431045294058697368227003439644839809767851608618176893080329310698605614360965199725467782498440632025216054915211615636225384587391692730438957640608095745163547903910786493577272")
	g2Gen.Y.SetString("9776954322322743732313986496143600875624250284803189422766662481205754219164790571343833769147920817555784777036310013265777572781104130018069386569670888186047120343466911735495135013382678",
		"16145243289410111745998998990108372459202724164009495870131936753931166691184115416492282110780957479364462561638319281671437384573226876857329949788153639802828393721896123566433408057468194")
	g2Gen.Z.SetString("1",
		"0")


	g1GenAff.FromW12(&g1Gen)
	g2GenAff.FromJacobian(&g2Gen)
	g1Infinity.X.SetOne()
	g1Infinity.Y.SetOne()
	g2Infinity.X.SetOne()
	g2Infinity.Y.SetOne()

    /*
	// r
	_r := fr.Modulus()

	thirdRootOneG1.SetString("4098895725012429242072311240482566844345873033931481129362557724405008256668293241245050359832461015092695507587185678086043587575438449040313411246717257958467499181450742260777082884928318")
	thirdRootOneG2.Square(&thirdRootOneG1)
	lambdaGLV.SetString("39705142672498995661671850106945620852186608752525090699191017895721506694646055668218723303426", 10)
	ecc.PrecomputeLattice(_r, &lambdaGLV, &glvBasis)
    */

    // TODO: correct values
	endo.u.A0.SetString("1")
	endo.v.A0.SetString("1")
}

// Generators return the generators of the r-torsion group, resp. in ker(pi-id), ker(Tr)
func Generators() (g1Jac G1W12, g2Jac G2Jac, g1Aff G1Affine, g2Aff G2Affine) {
	g1Aff = g1GenAff
	g2Aff = g2GenAff
	g1Jac = g1Gen
	g2Jac = g2Gen
	return
}
