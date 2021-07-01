package fptower

import (
	"math/bits"
)

// Expt set z to x^t in E12 and return z
// TODO: optimized Expt
func (z *E12) Expt(x *E12) *E12 {
	const tAbsVal uint64 = 11170052975785672705

	var result E12
	result.Set(x)

	l := bits.Len64(tAbsVal) - 2
	for i := l; i >= 0; i-- {
		result.CyclotomicSquare(&result)
		if tAbsVal&(1<<uint(i)) != 0 {
			result.Mul(&result, x)
		}
	}

	z.Set(&result)
	return z
}

// MulBy034 multiplication by sparse element (c0,0,0,c3,c4,0)
func (z *E12) MulBy034(c0, c3, c4 *E2) *E12 {

	var a, b, d E6

	a.MulByE2(&z.C0, c0)

	b.Set(&z.C1)
	b.MulBy01(c3, c4)

	c0.Add(c0, c3)
	d.Add(&z.C0, &z.C1)
	d.MulBy01(c0, c4)

	z.C1.Add(&a, &b).Neg(&z.C1).Add(&z.C1, &d)
	z.C0.MulByNonResidue(&b).Add(&z.C0, &a)

	return z
}
