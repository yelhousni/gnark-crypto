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
	"errors"

	"github.com/consensys/gnark-crypto/ecc/cp8-632/fp"
	"github.com/consensys/gnark-crypto/ecc/cp8-632/internal/fptower"
)

// GT target group of the pairing
type GT = fptower.E8

type lineEvaluation struct {
	r0 fptower.E2
	r1 fptower.E2
	r2 fptower.E2
}

// Pair calculates the reduced pairing for a set of points
func Pair(P []G1Affine, Q []G2Affine) (GT, error) {
	f, err := MillerLoop(P, Q)
	if err != nil {
		return GT{}, err
	}
	return FinalExponentiation(&f), nil
}

// FinalExponentiation computes the final expo x**(c*(p**3-1)(p+1)(p**2-p+1)/r)
func FinalExponentiation(z *GT, _z ...*GT) GT {

	var result GT
	result.Set(z)

	for _, e := range _z {
		result.Mul(&result, e)
	}

    result.Exp(&result, finalExponent)

	return result
}

// PairingCheck calculates the reduced pairing for a set of points and returns True if the result is One
func PairingCheck(P []G1Affine, Q []G2Affine) (bool, error) {
	f, err := Pair(P, Q)
	if err != nil {
		return false, err
	}
	var one GT
	one.SetOne()
	return f.Equal(&one), nil
}

// MillerLoop Miller loop
func MillerLoop(P []G1Affine, Q []G2Affine) (GT, error) {
	return MillerLoopTate(P, Q)
}

// MillerLoop Tate
func MillerLoopTate(P []G1Affine, Q []G2Affine) (GT, error) {
	// check input size match
	n := len(P)
	if n == 0 || n != len(Q) {
		return GT{}, errors.New("invalid inputs sizes")
	}

	// filter infinity points
	p := make([]G1Affine, 0, n)
	q := make([]G2Affine, 0, n)

	for k := 0; k < n; k++ {
		if P[k].IsInfinity() || Q[k].IsInfinity() {
			continue
		}
		p = append(p, P[k])
		q = append(q, Q[k])
	}

	n = len(p)

	// projective points for P
	pW12 := make([]g1W12, n)
	pNeg := make([]G1Affine, n)
	for k := 0; k < n; k++ {
		pW12[k].FromAffine(&p[k])
		pNeg[k].Neg(&p[k])
	}

	// f_{r,P}(Q)
	var result GT
	result.SetOne()
	var l lineEvaluation

	for i := 313; i >= 0; i-- {
		result.Square(&result)

		for k := 0; k < n; k++ {
			pW12[k].DoubleStep(&l)
			// line evaluation
			l.r1.Mul(&l.r1, &q[k].Y)
			l.r2.Mul(&l.r2, &q[k].X)
			result.MulBy023(&l.r0, &l.r1, &l.r2)

			if loopCounterTate[i] == 1 {
				pW12[k].AddMixedStep(&l, &p[k])
				// line evaluation
				l.r1.Mul(&l.r1, &q[k].Y)
				l.r2.Mul(&l.r2, &q[k].X)
				result.MulBy023(&l.r0, &l.r1, &l.r2)

			} else if loopCounterTate[i] == -1 {
				pW12[k].AddMixedStep(&l, &pNeg[k])
				// line evaluation
				l.r1.Mul(&l.r1, &q[k].Y)
				l.r2.Mul(&l.r2, &q[k].X)
				result.MulBy023(&l.r0, &l.r1, &l.r2)
			}
		}
	}

	return result, nil
}

// For Tate pairing
// DoubleStep doubles a point in W12 coordinates, and evaluates the line in Miller loop
// https://eprint.iacr.org/2009/615.pdf (section 4)
func (p *g1W12) DoubleStep(evaluations *lineEvaluation) {

	var A, B, C, D, E, F fp.Element
	A.Square(&p.x)
	B.Square(&p.y)
	C.Square(&p.z)
    D.Mul(&C, &aCurveCoeff)

    evaluations.r0.A0.Add(&p.x, &A).
        Sub(&evaluations.r0.A0, &D).
        Square(&evaluations.r0.A0).
        Sub(&evaluations.r0.A0, &A)

	p.x.Sub(&A, &D).
		Square(&p.x)

    evaluations.r0.A0.Sub(&evaluations.r0.A0, &p.x)

    evaluations.r1.A0.Add(&p.y, &p.z).
        Square(&evaluations.r1.A0).
        Sub(&evaluations.r1.A0, &B).
        Sub(&evaluations.r1.A0, &C)

    evaluations.r2.A0.Double(&A).
        Add(&evaluations.r2.A0, &A).
        Add(&evaluations.r2.A0, &D).
        Mul(&evaluations.r2.A0, &p.z).
        Double(&evaluations.r2.A0).
        Neg(&evaluations.r2.A0)

	E.Add(&A, &D).
        Square(&E).
        Double(&E).
        Sub(&E, &p.x)
    F.Sub(&A, &D).
        Add(&F, &p.y).
        Square(&F).
        Sub(&F, &B).
        Sub(&F, &p.x)

	p.y.Mul(&E, &F)

	p.z.Double(&B).
        Double(&p.z)
}

// AddMixedStep point addition in Mixed Homogenous projective and Affine coordinates
// https://eprint.iacr.org/2009/615.pdf (section 4)
func (p *g1W12) AddMixedStep(evaluations *lineEvaluation, a *G1Affine) {

	// get some Element from our pool
	var tmp, A, one, C, D, E, F, G, H, I, II, J, K fp.Element
    one.SetOne()
    A.Square(&p.z)
    C.Add(&p.z, &one).
        Square(&C).
        Sub(&C, &A).
        Sub(&C, &one)
    D.Set(&p.x)
    E.Mul(&p.z, &a.X)
    F.Set(&p.y)
    G.Mul(&A, &a.Y)
    H.Sub(&D, &E)
    I.Sub(&F, &G).
        Double(&I)
    II.Square(&I)
    J.Mul(&C, &H)
    K.Mul(&J, &H).
        Double(&K).
        Double(&K)
    tmp.Add(&D, &E).
        Mul(&tmp, &K)

    p.x.Double(&II).
        Sub(&p.x, &tmp)

    p.z.Square(&J)

    tmp.Mul(&D, &K).
        Sub(&tmp, &p.x)

    p.y.Add(&J, &I).
        Square(&p.y).
        Sub(&p.y, &p.z).
        Sub(&p.y, &II).
        Mul(&p.y, &tmp)

    tmp.Square(&K).
        Mul(&tmp, &F)

    p.y.Sub(&p.y, &tmp)

    p.z.Double(&p.z)

    tmp.Mul(&J, &a.Y)

    evaluations.r0.A0.Mul(&I, &a.X).
        Sub(&evaluations.r0.A0, &tmp)
    evaluations.r1.A0.Set(&J)
    evaluations.r2.A0.Neg(&I)
}
