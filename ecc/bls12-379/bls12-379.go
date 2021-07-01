package bls12379

import (
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/consensys/gnark-crypto/ecc/bls12-379/fp"
	"github.com/consensys/gnark-crypto/ecc/bls12-379/fr"
	"github.com/consensys/gnark-crypto/ecc/bls12-379/internal/fptower"
)

// E: y**2=x**3+1
// Etwist: y**2 = x**3+(u+5)**-1
// Tower: Fp->Fp2, u**2=5 -> Fp12, v**6=u+5
// Generator (BLS12 family): x=11170052975785672705
// optimal Ate loop: trace(frob)-1=x
// trace of pi: x+1
// Fp: p=647455824720115791999401948377863964948475498868568187799869012112522291430466516542182303552703940519176779595777
// Fr: r=15567573732069904898445906795858855143537221806202697669674256244108007833601

// ID bls379 ID
const ID = ecc.BLS12_379

// bCurveCoeff b coeff of the curve
var bCurveCoeff fp.Element

// twist
var twist fptower.E2

// bTwistCurveCoeff b coeff of the twist (defined over Fp2) curve
var bTwistCurveCoeff fptower.E2

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

// optimal Ate loop counter (=trace-1 = x in BLS family)
var loopCounter [64]int8

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

// psi o pi o psi**-1, where psi:E->E' is the degree 6 iso defined over Fp12
var endo struct {
	u fptower.E2
	v fptower.E2
}

// generator of the curve
var xGen big.Int

// expose the tower -- github.com/consensys/gnark uses it in a gnark circuit

// E2 is a degree two finite field extension of fp.Element
type E2 = fptower.E2

// E6 is a degree three finite field extension of fp2
type E6 = fptower.E6

// E12 is a degree two finite field extension of fp6
type E12 = fptower.E12

func init() {

	bCurveCoeff.SetUint64(1)
	twist.A1.SetUint64(1)
	twist.A0.SetUint64(5)
	bTwistCurveCoeff.Inverse(&twist)

	twoInv.SetOne().Double(&twoInv).Inverse(&twoInv)

	g1Gen.X.SetString("132523964661110384907332505045337948264398196315391712221656804725644620418869173495434061941186407932288361189427")
	g1Gen.Y.SetString("489540456653765585587483238625648106453748411558834855892046639881836253922937210183166545003353426582699245889888")
	g1Gen.Z.SetString("1")

	g2Gen.X.SetString("63536250141944042999329934468565286750304783693126395492168700859056213268331021509996747553200477568653273815834",
		"543911846289311625681475737404806809310985050399821208075813618556502455423364557877980547865120514429405372728749")
	g2Gen.Y.SetString("218437454626769505098282817930895835135474259385784458258055745507127879418351714828772994086735819028645374589000",
		"190125814215625101662261605981382591306631856523959927924837144028510735770756689374464298081356205177232360494197")
	g2Gen.Z.SetString("1",
		"0")

	g1GenAff.FromJacobian(&g1Gen)
	g2GenAff.FromJacobian(&g2Gen)

	g1Infinity.X.SetOne()
	g1Infinity.Y.SetOne()
	g2Infinity.X.SetOne()
	g2Infinity.Y.SetOne()

	thirdRootOneG1.SetString("173890623291670311836512360316100953774879339360001038672113943180607164901767459447515297873921") // x**5-3*x**4+3*x**3-x+1
	thirdRootOneG2.Square(&thirdRootOneG1)
	lambdaGLV.SetString("124770083481858362097340376349382017024", 10) //(x**2-1)
	_r := fr.Modulus()
	ecc.PrecomputeLattice(_r, &lambdaGLV, &glvBasis)

	endo.u.A0.SetString("222745979508143237255341693926393623218761351828393355617642805276543242213727214918790809594387987052224624366927")
	endo.u.A1.SetString("53608708391251269430337404591046336897908586198581852257006640780158581836659927408814995415886520684560049511630")
	endo.v.A0.SetString("205343422135153100463306508509910492825903335462897260204806706089537866940279563177678964270966750592224060664130")
	endo.v.A1.SetString("205343422135153100463306508509910492825903335462897260204806706089537866940279563177678964270966750592224060664130")

	// binary decomposition of 11170052975785672705 little endian
	loopCounter = [64]int8{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1}

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
