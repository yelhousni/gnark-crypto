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

// FinalExponentiation computes the final expo x**((p**4-1)(p**4+1)/r)
func FinalExponentiation(z *GT, _z ...*GT) GT {

	var result GT
	result.Set(z)

	for _, e := range _z {
		result.Mul(&result, e)
	}

	var t GT

	// easy part
	t.Conjugate(&result)
	result.Inverse(&result)
	result.Mul(&t, &result)

	// hard part (up to permutation)
	// TODO: use cyclo squares
	// TODO: use a short addition chain
	result.Exp(&result, finalExponentHardPart)

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
	// return MillerLoopAte(P, Q)
	return MillerLoopOptAte(P, Q)
}

func MillerLoopOptAte(P []G1Affine, Q []G2Affine) (GT, error) {
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

	var f GT
	f.SetOne()
	for k := 0; k < n; k++ {
		m, err := MillerLoopOptAteSingle(p[k], q[k])
		if err != nil {
			return GT{}, err
		}
		f.Mul(&f, &m)
	}

	return f, nil
}

// MillerLoopOptAteSingle optimal ate
func MillerLoopOptAteSingle(p G1Affine, q G2Affine) (GT, error) {

	// W12 points for Q
	var qW12, qW12_a0, qFrobSquare g2W12
	var qNeg, qFrob, qFrobCube G2Affine
	qW12.FromAffine(&q)
	qNeg.Neg(&q)

	var result [4]GT
	var l, l1, l2 lineEvaluation

	// f_{a0,Q}(P)
	result[0].SetOne()

	for i := 77; i >= 0; i-- {
		result[0].Square(&result[0])

		qW12.DoubleStep(&l)
		// line evaluation
		l.r1.MulByElement(&l.r1, &p.Y)
		l.r2.MulByElement(&l.r2, &p.X)
		result[0].MulBy023(&l.r1, &l.r2, &l.r0)

		if loopCounterOptAte0[i] == 1 {
			qW12.AddMixedStep(&l, &q)
			// line evaluation
			l.r1.MulByElement(&l.r1, &p.Y)
			l.r2.MulByElement(&l.r2, &p.X)
			result[0].MulBy023(&l.r1, &l.r2, &l.r0)

		} else if loopCounterOptAte0[i] == -1 {
			qW12.AddMixedStep(&l, &qNeg)
			// line evaluation
			l.r1.MulByElement(&l.r1, &p.Y)
			l.r2.MulByElement(&l.r2, &p.X)
			result[0].MulBy023(&l.r1, &l.r2, &l.r0)
		}
	}

	result[0].Conjugate(&result[0])

	qW12_a0.Neg(&qW12)

	qW12.FromAffine(&q)

	// f_{a1,Q}(P)
	result[1].SetOne()

	for i := 76; i >= 0; i-- {
		result[1].Square(&result[1])

		qW12.DoubleStep(&l)
		// line evaluation
		l.r1.MulByElement(&l.r1, &p.Y)
		l.r2.MulByElement(&l.r2, &p.X)
		result[1].MulBy023(&l.r1, &l.r2, &l.r0)

		if loopCounterOptAte1[i] == 1 {
			qW12.AddMixedStep(&l, &q)
			// line evaluation
			l.r1.MulByElement(&l.r1, &p.Y)
			l.r2.MulByElement(&l.r2, &p.X)
			result[1].MulBy023(&l.r1, &l.r2, &l.r0)

		} else if loopCounterOptAte1[i] == -1 {
			qW12.AddMixedStep(&l, &qNeg)
			// line evaluation
			l.r1.MulByElement(&l.r1, &p.Y)
			l.r2.MulByElement(&l.r2, &p.X)
			result[1].MulBy023(&l.r1, &l.r2, &l.r0)
		}
	}

	result[1].Frobenius(&result[1])

	qFrob.FromW12(&qW12)
	qFrob.Frobenius(&qFrob)

	qW12.FromAffine(&q)

	// f_{a2,Q}(P)
	result[2].SetOne()

	for i := 76; i >= 0; i-- {
		result[2].Square(&result[2])

		qW12.DoubleStep(&l)
		// line evaluation
		l.r1.MulByElement(&l.r1, &p.Y)
		l.r2.MulByElement(&l.r2, &p.X)
		result[2].MulBy023(&l.r1, &l.r2, &l.r0)

		if loopCounterOptAte2[i] == 1 {
			qW12.AddMixedStep(&l, &q)
			// line evaluation
			l.r1.MulByElement(&l.r1, &p.Y)
			l.r2.MulByElement(&l.r2, &p.X)
			result[2].MulBy023(&l.r1, &l.r2, &l.r0)

		} else if loopCounterOptAte2[i] == -1 {
			qW12.AddMixedStep(&l, &qNeg)
			// line evaluation
			l.r1.MulByElement(&l.r1, &p.Y)
			l.r2.MulByElement(&l.r2, &p.X)
			result[2].MulBy023(&l.r1, &l.r2, &l.r0)
		}
	}

	result[2].FrobeniusSquare(&result[2]).
		Conjugate(&result[2])

	qFrobSquare.Neg(&qW12)
	qFrobSquare.FrobeniusSquare(&qFrobSquare)

	qW12.FromAffine(&q)

	// f_{a3,Q}(P)
	result[3].SetOne()

	for i := 77; i >= 0; i-- {
		result[3].Square(&result[3])

		qW12.DoubleStep(&l)
		// line evaluation
		l.r1.MulByElement(&l.r1, &p.Y)
		l.r2.MulByElement(&l.r2, &p.X)
		result[3].MulBy023(&l.r1, &l.r2, &l.r0)

		if loopCounterOptAte3[i] == 1 {
			qW12.AddMixedStep(&l, &q)
			// line evaluation
			l.r1.MulByElement(&l.r1, &p.Y)
			l.r2.MulByElement(&l.r2, &p.X)
			result[3].MulBy023(&l.r1, &l.r2, &l.r0)

		} else if loopCounterOptAte3[i] == -1 {
			qW12.AddMixedStep(&l, &qNeg)
			// line evaluation
			l.r1.MulByElement(&l.r1, &p.Y)
			l.r2.MulByElement(&l.r2, &p.X)
			result[3].MulBy023(&l.r1, &l.r2, &l.r0)
		}
	}

	result[3].FrobeniusCube(&result[3]).
		Conjugate(&result[3])

	qW12.Neg(&qW12)
	qFrobCube.FromW12(&qW12)
	qFrobCube.FrobeniusCube(&qFrobCube)

	// l_{a0*Q, a1*pi(Q)}(P)
	qW12_a0.AddMixedStep(&l1, &qFrob)
	l1.r1.MulByElement(&l1.r1, &p.Y)
	l1.r2.MulByElement(&l1.r2, &p.X)
	result[3].MulBy023(&l1.r1, &l1.r2, &l1.r0)

	// l_{a2*pi^2(Q), a3*pi^3(Q)}(P)
	qFrobSquare.AddMixedStep(&l2, &qFrobCube)
	l2.r1.MulByElement(&l2.r1, &p.Y)
	l2.r2.MulByElement(&l2.r2, &p.X)
	result[3].MulBy023(&l2.r1, &l2.r2, &l2.r0)

	// f0 * f1^q * f2^q2 * f3^q3 * l1 * l2
	result[3].Mul(&result[3], &result[2]).
		Mul(&result[3], &result[1]).
		Mul(&result[3], &result[0])

	return result[3], nil
}

// MillerLoopAte ate
func MillerLoopAte(P []G1Affine, Q []G2Affine) (GT, error) {
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

	// W12 points for Q
	qW12 := make([]g2W12, n)
	qNeg := make([]G2Affine, n)
	for k := 0; k < n; k++ {
		qW12[k].FromAffine(&q[k])
		qNeg[k].Neg(&q[k])
	}

	// f_{t-1,Q}(P)
	var result GT
	result.SetOne()
	var l lineEvaluation

	for i := 316; i >= 0; i-- {
		result.Square(&result)

		for k := 0; k < n; k++ {
			qW12[k].DoubleStep(&l)
			// line evaluation
			l.r1.MulByElement(&l.r1, &p[k].Y)
			l.r2.MulByElement(&l.r2, &p[k].X)
			result.MulBy023(&l.r1, &l.r2, &l.r0)

			if loopCounterAte[i] == 1 {
				qW12[k].AddMixedStep(&l, &q[k])
				// line evaluation
				l.r1.MulByElement(&l.r1, &p[k].Y)
				l.r2.MulByElement(&l.r2, &p[k].X)
				result.MulBy023(&l.r1, &l.r2, &l.r0)

			} else if loopCounterAte[i] == -1 {
				qW12[k].AddMixedStep(&l, &qNeg[k])
				// line evaluation
				l.r1.MulByElement(&l.r1, &p[k].Y)
				l.r2.MulByElement(&l.r2, &p[k].X)
				result.MulBy023(&l.r1, &l.r2, &l.r0)
			}
		}
	}

	return result, nil
}

// DoubleStep doubles a point in W12 coordinates, and evaluates the line in Miller loop
// https://eprint.iacr.org/2009/615.pdf (section 4)
func (p *g2W12) DoubleStep(evaluations *lineEvaluation) {

	var A, B, C, D, E, F fptower.E2
	A.Square(&p.x)
	B.Square(&p.y)
	C.Square(&p.z)
	D.Mul(&C, &aTwistCurveCoeff)

	evaluations.r0.Add(&p.x, &A).
		Sub(&evaluations.r0, &D).
		Square(&evaluations.r0).
		Sub(&evaluations.r0, &A)

	p.x.Sub(&A, &D).
		Square(&p.x)

	evaluations.r0.Sub(&evaluations.r0, &p.x)

	evaluations.r1.Add(&p.y, &p.z).
		Square(&evaluations.r1).
		Sub(&evaluations.r1, &B).
		Sub(&evaluations.r1, &C).
		Double(&evaluations.r1)

	evaluations.r2.Double(&A).
		Add(&evaluations.r2, &A).
		Add(&evaluations.r2, &D).
		Mul(&evaluations.r2, &p.z).
		Double(&evaluations.r2).
		Neg(&evaluations.r2)

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

// AddMixedStep point addition in Mixed W12 and Affine coordinates
// https://eprint.iacr.org/2009/615.pdf (section 4)
func (p *g2W12) AddMixedStep(evaluations *lineEvaluation, a *G2Affine) {

	// get some Element from our pool
	var tmp, A, one, C, D, E, F, G, H, I, II, J, K fptower.E2
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

	evaluations.r0.Mul(&I, &a.X).
		Sub(&evaluations.r0, &tmp)
	evaluations.r1.Set(&J)
	evaluations.r2.Neg(&I)
}
