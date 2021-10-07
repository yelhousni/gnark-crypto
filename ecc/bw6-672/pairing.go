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
	"errors"

	"github.com/consensys/gnark-crypto/ecc/bw6-672/fp"
	"github.com/consensys/gnark-crypto/ecc/bw6-672/internal/fptower"
)

// GT target group of the pairing
type GT = fptower.E6

type lineEvaluation struct {
	r0 fp.Element
	r1 fp.Element
	r2 fp.Element
}

// Pair calculates the reduced pairing for a set of points
func Pair(P []G1Affine, Q []G2Affine) (GT, error) {
	f, err := MillerLoop(P, Q)
	if err != nil {
		return GT{}, err
	}
	return FinalExponentiation(&f), nil
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

// FinalExponentiation computes the final expo x**(c*(p**3-1)(p+1)(p**2-p+1)/r)
func FinalExponentiation(z *GT, _z ...*GT) GT {

	var result GT
	result.Set(z)

	for _, e := range _z {
		result.Mul(&result, e)
	}

	var buf GT

	// easy part exponent: (p**3 - 1)*(p+1)
	buf.Conjugate(&result)
	result.Inverse(&result)
	buf.Mul(&buf, &result)
	result.Frobenius(&buf).
		Mul(&result, &buf)

	// hard part exponent: 3(u+1)(p^2-p+1)/r
	var m [11]GT
	var f10, _m1, _m3, _m4, _m5, _m7, _m8, _m8m5m0, _m6, f11, f11f10, f12, f1, f1u, f1q, f1a GT
	m[0].Set(&result)
	for i := 1; i < 11; i++ {
		m[i].Expt(&m[i-1])
	}
	result.Mul(&m[3], &m[1]).
		Conjugate(&result).
		Mul(&result, &m[2]).
		CyclotomicSquare(&result).
		Mul(&result, &m[0]).
		Mul(&result, &m[4])
	buf.Frobenius(&m[0])
	result.Mul(&result, &buf)
	buf.CyclotomicSquare(&result).
		CyclotomicSquare(&buf).
		CyclotomicSquare(&buf)
	result.Mul(&result, &buf)
	_m1.Conjugate(&m[1])
	_m3.Conjugate(&m[3])
	_m4.Conjugate(&m[4])
	_m5.Conjugate(&m[5])
	_m7.Conjugate(&m[7])
	f10.Mul(&m[4], &_m3).
		CyclotomicSquare(&f10).
		Mul(&f10, &m[2]).
		Mul(&f10, &m[6]).
		Mul(&f10, &_m5).
		CyclotomicSquare(&f10).
		Mul(&f10, &_m1).
		Mul(&f10, &_m5).
		Mul(&f10, &_m7).
		CyclotomicSquare(&f10).
		Mul(&f10, &m[0]).
		Mul(&f10, &m[2]).
		Mul(&f10, &m[3]).
		Mul(&f10, &_m1).
		CyclotomicSquare(&f10).
		Mul(&f10, &m[0]).
		Mul(&f10, &m[8]).
		Mul(&f10, &_m4)
	_m8.Conjugate(&m[8])
	_m6.Conjugate(&m[6])
	_m8m5m0.Mul(&m[5], &_m8).
		Mul(&_m8m5m0, &m[0])
	f11.Mul(&m[7], &_m6).
		CyclotomicSquare(&f11).
		Mul(&f11, &m[2]).
		Mul(&f11, &_m3).
		Mul(&f11, &_m8m5m0).
		CyclotomicSquare(&f11).
		Mul(&f11, &_m8m5m0).
		Mul(&f11, &m[9]).
		Mul(&f11, &_m1)
	f11.Expc1(&f11) // f11^(ht/2)
	f11f10.Mul(&f11, &f10)
	buf.CyclotomicSquare(&f11f10)
	f11f10.Mul(&f11f10, &buf)
	f12.Mul(&m[0], &m[1]).
		Mul(&f12, &m[2]).
		Mul(&f12, &m[8]).
		Mul(&f12, &m[10])
	buf.CyclotomicSquare(&m[5])
	f12.Mul(&f12, &buf)
	buf.CyclotomicSquare(&m[9]).
		Mul(&buf, &m[6]).
		Mul(&buf, &m[4]).
		Conjugate(&buf)
	f12.Mul(&f12, &buf)
	f12.Expc2(&f12) // f12^(ht/2)^2
	f1.Mul(&f11f10, &f12)
	f1u.Expt(&f1)
	f1q.Mul(&f1u, &f1).
		Frobenius(&f1q)
	buf.Conjugate(&f1u)
	f1a.Expt(&f1u).
		Mul(&f1a, &buf).
		Expt(&f1a).
		Expt(&f1a).
		Expt(&f1a).
		Mul(&f1a, &buf)

	result.Mul(&result, &f1a).
		Mul(&result, &f1q)

	return result
}

// MillerLoop Miller loop
func MillerLoop(P []G1Affine, Q []G2Affine) (GT, error) {
	// return MillerLoopOptAte(P, Q)
	return MillerLoopOptTate(P, Q)
}

// MillerLoop Optimal Tate (or twisted ate or Eta revisited)
// Alg.2 in ...
func MillerLoopOptTate(P []G1Affine, Q []G2Affine) (GT, error) {
	// check input size match
	n := len(P)
	if n == 0 || n != len(Q) {
		return GT{}, errors.New("invalid inputs sizes")
	}

	// filter infinity points
	p0 := make([]G1Affine, 0, n)
	q := make([]G2Affine, 0, n)

	for k := 0; k < n; k++ {
		if P[k].IsInfinity() || Q[k].IsInfinity() {
			continue
		}
		p0 = append(p0, P[k])
		q = append(q, Q[k])
	}

	n = len(q)

	// precomputations
	pProj1 := make([]g1Proj, n)
	p1 := make([]G1Affine, n)
	p01 := make([]G1Affine, n)
	pProj01 := make([]g1Proj, n)
	l01 := make([]lineEvaluation, n)
	for k := 0; k < n; k++ {
		p0[k].Neg(&p0[k])
		p1[k].phi(&p0[k]).Add(&p1[k], &p0[k])
		pProj1[k].FromAffine(&p1[k])
		// l_{p0,p1}(q)
		pProj01[k].Set(&pProj1[k])
		pProj01[k].AddMixedStep(&l01[k], &p0[k])
		l01[k].r1.Mul(&l01[k].r1, &q[k].X)
		l01[k].r2.Mul(&l01[k].r2, &q[k].Y)
	}
	BatchProjectiveToAffineG1(pProj01, p01)

	// f_{a0+lambda*a1,P}(Q)
	var result, ss GT
	result.SetOne()
	var l, l0 lineEvaluation

	for i := 156; i >= 0; i-- {
		result.Square(&result)

		for k := 0; k < n; k++ {
			pProj1[k].DoubleStep(&l0)
			// line evaluation
			l0.r1.Mul(&l0.r1, &q[k].X)
			l0.r2.Mul(&l0.r2, &q[k].Y)

			if loopCounterOptTate0[i] == 0 && loopCounterOptTate1[i] == 0 {
				result.MulBy014(&l0.r0, &l0.r1, &l0.r2)
			} else if loopCounterOptTate0[i] == 1 && loopCounterOptTate1[i] == 1 {
				pProj1[k].AddMixedStep(&l, &p01[k])
				l.r1.Mul(&l.r1, &q[k].X)
				l.r2.Mul(&l.r2, &q[k].Y)
				ss.Mul014By014(&l.r0, &l.r1, &l.r2, &l01[k].r0, &l01[k].r1, &l01[k].r2)
				result.MulBy014(&l0.r0, &l0.r1, &l0.r2).
					Mul(&result, &ss)
			} else if loopCounterOptTate0[i] == 1 && loopCounterOptTate1[i] == 0 {
				pProj1[k].AddMixedStep(&l, &p0[k])
				l.r1.Mul(&l.r1, &q[k].X)
				l.r2.Mul(&l.r2, &q[k].Y)
				ss.Mul014By014(&l.r0, &l.r1, &l.r2, &l0.r0, &l0.r1, &l0.r2)
				result.Mul(&result, &ss)
			} else {
				pProj1[k].AddMixedStep(&l, &p1[k])
				l.r1.Mul(&l.r1, &q[k].X)
				l.r2.Mul(&l.r2, &q[k].Y)
				ss.Mul014By014(&l.r0, &l.r1, &l.r2, &l0.r0, &l0.r1, &l0.r2)
				result.Mul(&result, &ss)
			}
		}
	}

	return result, nil
}

// MillerLoop Optimal ate
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

	// projective points for Q
	qProj1 := make([]g2Proj, n)
	qProj2 := make([]g2Proj, n)
	qNeg := make([]G2Affine, n)
	for k := 0; k < n; k++ {
		qProj1[k].FromAffine(&q[k])
		qProj2[k].FromAffine(&q[k])
		qNeg[k].Neg(&q[k])
	}

	// f_{u+1,Q}(P)
	var result1 GT
	result1.SetOne()
	var l lineEvaluation

	for i := 31; i >= 0; i-- {
		result1.Square(&result1)

		for k := 0; k < n; k++ {
			qProj1[k].DoubleStep(&l)
			// line evaluation
			l.r0.Mul(&l.r0, &p[k].Y)
			l.r1.Mul(&l.r1, &p[k].X)
			result1.MulBy034(&l.r0, &l.r1, &l.r2)

			if loopCounterAte0[i] == 1 {
				qProj1[k].AddMixedStep(&l, &q[k])
				// line evaluation
				l.r0.Mul(&l.r0, &p[k].Y)
				l.r1.Mul(&l.r1, &p[k].X)
				result1.MulBy034(&l.r0, &l.r1, &l.r2)

			} else if loopCounterAte0[i] == -1 {
				qProj1[k].AddMixedStep(&l, &qNeg[k])
				// line evaluation
				l.r0.Mul(&l.r0, &p[k].Y)
				l.r1.Mul(&l.r1, &p[k].X)
				result1.MulBy034(&l.r0, &l.r1, &l.r2)
			}
		}
	}

	result1.Conjugate(&result1)

	// f_{u^5-u^4-u,Q}(P)
	var result2 GT
	result2.SetOne()

	for i := 157; i >= 0; i-- {
		result2.Square(&result2)

		for k := 0; k < n; k++ {
			qProj2[k].DoubleStep(&l)
			// line evaluation
			l.r0.Mul(&l.r0, &p[k].Y)
			l.r1.Mul(&l.r1, &p[k].X)
			result2.MulBy034(&l.r0, &l.r1, &l.r2)

			if loopCounterAte1[i] == 1 {
				qProj2[k].AddMixedStep(&l, &q[k])
				// line evaluation
				l.r0.Mul(&l.r0, &p[k].Y)
				l.r1.Mul(&l.r1, &p[k].X)
				result2.MulBy034(&l.r0, &l.r1, &l.r2)

			} else if loopCounterAte1[i] == -1 {
				qProj2[k].AddMixedStep(&l, &qNeg[k])
				// line evaluation
				l.r0.Mul(&l.r0, &p[k].Y)
				l.r1.Mul(&l.r1, &p[k].X)
				result2.MulBy034(&l.r0, &l.r1, &l.r2)
			}
		}
	}

	result2.Conjugate(&result2)

	result2.Frobenius(&result2).
		Mul(&result2, &result1)

	return result2, nil
}

// For Tate pairing
// DoubleStep doubles a point in Homogenous projective coordinates, and evaluates the line in Miller loop
// https://eprint.iacr.org/2013/722.pdf (Section 4.3)
func (p *g1Proj) DoubleStep(evaluations *lineEvaluation) {

	// get some Element from our pool
	var t0, t1, A, B, C, D, E, EE, F, G, H, I, J, K fp.Element
	t0.Mul(&p.x, &p.y)
	A.Mul(&t0, &twoInv)
	B.Square(&p.y)
	C.Square(&p.z)
	D.Double(&C).
		Add(&D, &C)
	E.Mul(&D, &bCurveCoeff)
	F.Double(&E).
		Add(&F, &E)
	G.Add(&B, &F)
	G.Mul(&G, &twoInv)
	H.Add(&p.y, &p.z).
		Square(&H)
	t1.Add(&B, &C)
	H.Sub(&H, &t1)
	I.Sub(&E, &B)
	J.Square(&p.x)
	EE.Square(&E)
	K.Double(&EE).
		Add(&K, &EE)

	// X, Y, Z
	p.x.Sub(&B, &F).
		Mul(&p.x, &A)
	p.y.Square(&G).
		Sub(&p.y, &K)
	p.z.Mul(&B, &H)

	// Line evaluation
	evaluations.r0.Set(&I)
	evaluations.r1.Double(&J).
		Add(&evaluations.r1, &J)
	evaluations.r2.Neg(&H)
}

// AddMixedStep point addition in Mixed Homogenous projective and Affine coordinates
// https://eprint.iacr.org/2013/722.pdf (Section 4.3)
func (p *g1Proj) AddMixedStep(evaluations *lineEvaluation, a *G1Affine) {

	// get some Element from our pool
	var Y2Z1, X2Z1, O, L, C, D, E, F, G, H, t0, t1, t2, J fp.Element
	Y2Z1.Mul(&a.Y, &p.z)
	O.Sub(&p.y, &Y2Z1)
	X2Z1.Mul(&a.X, &p.z)
	L.Sub(&p.x, &X2Z1)
	C.Square(&O)
	D.Square(&L)
	E.Mul(&L, &D)
	F.Mul(&p.z, &C)
	G.Mul(&p.x, &D)
	t0.Double(&G)
	H.Add(&E, &F).
		Sub(&H, &t0)
	t1.Mul(&p.y, &E)

	// X, Y, Z
	p.x.Mul(&L, &H)
	p.y.Sub(&G, &H).
		Mul(&p.y, &O).
		Sub(&p.y, &t1)
	p.z.Mul(&E, &p.z)

	t2.Mul(&L, &a.Y)
	J.Mul(&a.X, &O).
		Sub(&J, &t2)

	// Line evaluation
	evaluations.r0.Set(&J)
	evaluations.r1.Neg(&O)
	evaluations.r2.Set(&L)
}

// For ate pairing
// DoubleStep doubles a point in Homogenous projective coordinates, and evaluates the line in Miller loop
// https://eprint.iacr.org/2013/722.pdf (Section 4.3)
func (p *g2Proj) DoubleStep(evaluations *lineEvaluation) {

	// get some Element from our pool
	var t0, t1, A, B, C, D, E, EE, F, G, H, I, J, K fp.Element
	t0.Mul(&p.x, &p.y)
	A.Mul(&t0, &twoInv)
	B.Square(&p.y)
	C.Square(&p.z)
	D.Double(&C).
		Add(&D, &C)
	E.Mul(&D, &bTwistCurveCoeff)
	F.Double(&E).
		Add(&F, &E)
	G.Add(&B, &F)
	G.Mul(&G, &twoInv)
	H.Add(&p.y, &p.z).
		Square(&H)
	t1.Add(&B, &C)
	H.Sub(&H, &t1)
	I.Sub(&E, &B)
	J.Square(&p.x)
	EE.Square(&E)
	K.Double(&EE).
		Add(&K, &EE)

	// X, Y, Z
	p.x.Sub(&B, &F).
		Mul(&p.x, &A)
	p.y.Square(&G).
		Sub(&p.y, &K)
	p.z.Mul(&B, &H)

	// Line evaluation
	evaluations.r0.Neg(&H)
	evaluations.r1.Double(&J).
		Add(&evaluations.r1, &J)
	evaluations.r2.Set(&I)
}

// AddMixedStep point addition in Mixed Homogenous projective and Affine coordinates
// https://eprint.iacr.org/2013/722.pdf (Section 4.3)
func (p *g2Proj) AddMixedStep(evaluations *lineEvaluation, a *G2Affine) {

	// get some Element from our pool
	var Y2Z1, X2Z1, O, L, C, D, E, F, G, H, t0, t1, t2, J fp.Element
	Y2Z1.Mul(&a.Y, &p.z)
	O.Sub(&p.y, &Y2Z1)
	X2Z1.Mul(&a.X, &p.z)
	L.Sub(&p.x, &X2Z1)
	C.Square(&O)
	D.Square(&L)
	E.Mul(&L, &D)
	F.Mul(&p.z, &C)
	G.Mul(&p.x, &D)
	t0.Double(&G)
	H.Add(&E, &F).
		Sub(&H, &t0)
	t1.Mul(&p.y, &E)

	// X, Y, Z
	p.x.Mul(&L, &H)
	p.y.Sub(&G, &H).
		Mul(&p.y, &O).
		Sub(&p.y, &t1)
	p.z.Mul(&E, &p.z)

	t2.Mul(&L, &a.Y)
	J.Mul(&a.X, &O).
		Sub(&J, &t2)

	// Line evaluation
	evaluations.r0.Set(&L)
	evaluations.r1.Neg(&O)
	evaluations.r2.Set(&J)
}
