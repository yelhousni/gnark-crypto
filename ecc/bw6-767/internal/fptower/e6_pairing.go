package fptower

import "github.com/consensys/gnark-crypto/ecc/bw6-767/fp"

func (z *E6) nSquare(n int) {
	for i := 0; i < n; i++ {
		z.CyclotomicSquare(z)
	}
}

// ExptHalf set z to x^(t/2) in E6 and return z
// const t/2 uint64 = 7566188111470821376 // negative
func (z *E6) ExptHalf(x *E6) *E6 {
	var result E6
	result.CyclotomicSquare(x)
	result.Mul(&result, x)
	result.nSquare(2)
	result.Mul(&result, x)
	result.nSquare(3)
	result.Mul(&result, x)
	result.nSquare(9)
	result.Mul(&result, x)
	result.nSquare(32)
	result.Mul(&result, x)
	result.nSquare(15)
	return z.Conjugate(&result) // because tAbsVal is negative
}

// Expt set z to x^t in E6 and return z
// const t uint64 = 15132376222941642752 // negative
func (z *E6) Expt(x *E6) *E6 {
	var result E6
	result.ExptHalf(x)
	return z.CyclotomicSquare(&result)
}

// Expc1 set z to x^c1 in E6 and return z
// ht, hy = -4, -6
// c1 = ht**2+3*hy**2 = 31 (11111)
func (z *E6) Expc1(x *E6) *E6 {

	var result E6

	result.CyclotomicSquare(x)
	result.Mul(&result, x)
	result.CyclotomicSquare(&result)
	result.Mul(&result, x)
	result.CyclotomicSquare(&result)
	result.Mul(&result, x)
	result.CyclotomicSquare(&result)
	result.Mul(&result, x)

	z.Set(&result)

	return z
}

// MulBy014 multiplication by sparse element (c0,c1,0,0,c4,0)
func (z *E6) MulBy014(c0, c1, c4 *fp.Element) *E6 {

	var a, b E3
	var d fp.Element

	a.Set(&z.B0)
	a.MulBy01(c0, c1)

	b.Set(&z.B1)
	b.MulBy1(c4)
	d.Add(c1, c4)

	z.B1.Add(&z.B1, &z.B0)
	z.B1.MulBy01(c0, &d)
	z.B1.Sub(&z.B1, &a)
	z.B1.Sub(&z.B1, &b)
	z.B0.MulByNonResidue(&b)
	z.B0.Add(&z.B0, &a)

	return z
}
