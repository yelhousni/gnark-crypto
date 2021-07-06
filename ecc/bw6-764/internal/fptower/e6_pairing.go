package fptower

import (
	"github.com/consensys/gnark-crypto/ecc/bw6-764/fp"
)

// Expt set z to x^t in E6 and return z
func (z *E6) Expt(x *E6) *E6 {
	// const tAbsVal uint64 = 11170052975785672705
	// tAbsVal in binary: 1001101100000100000000000000000000000000000000000000000000000001
	// drop the low 50 bits (all 0 except the least significant bit): 10011011000001 = 9921
	// Shortest addition chains can be found at https://wwwhomes.uni-bielefeld.de/achim/addition_chain.html
	var result, x4, x5 E6
	result.Set(x)                    //  0                1
	result.CyclotomicSquare(&result) //  1( 0)            2
	result.CyclotomicSquare(&result) //  2( 1)            4
	x4.Set(&result)                  //  save x4 for step 4
	result.Mul(&result, x)           //  3( 2, 0)         5
	x5.Set(&result)                  //  save x5 for step 8
	result.Mul(&result, &x4)         //  4( 3, 2)         9
	result.CyclotomicSquare(&result) //  5( 4)           18
	result.CyclotomicSquare(&result) //  6( 5)           36
	result.CyclotomicSquare(&result) //  7( 6)           72
	result.Mul(&result, &x5)         //  8( 7, 3)        77
	result.CyclotomicSquare(&result) //  9( 8)          154
	result.Mul(&result, x)           // 10( 9, 0)       155
	result.CyclotomicSquare(&result) // 11(10)          310
	result.CyclotomicSquare(&result) // 12(11)          620
	result.CyclotomicSquare(&result) // 13(12)         1240
	result.CyclotomicSquare(&result) // 14(13)         2480
	result.CyclotomicSquare(&result) // 15(14)         4960
	result.CyclotomicSquare(&result) // 16(15)         9920
	result.Mul(&result, x)           // 17(16, 0)      9921

	// the remaining 50 bits
	for i := 0; i < 50; i++ {
		result.CyclotomicSquare(&result)
	}
	result.Mul(&result, x)

	z.Set(&result)
	return z
}

// Expc2 set z to x^c2 in E6 and return z
// ht, hy = -25, 3
// c2 = ht+hy = -22 (10110)
func (z *E6) Expc2(x *E6) *E6 {

	var result E6

	result.CyclotomicSquare(x)
	result.CyclotomicSquare(&result)
	result.Mul(&result, x)
	result.CyclotomicSquare(&result)
	result.Mul(&result, x)
	result.CyclotomicSquare(&result)

	z.Conjugate(&result)

	return z
}

// Expc1 set z to x^c1 in E6 and return z
// ht, hy = -25, 3
// c1 = ht**2+3*hy**2 = 652 (1010001100)
func (z *E6) Expc1(x *E6) *E6 {

	var result E6

	result.CyclotomicSquare(x)
	result.CyclotomicSquare(&result)
	result.Mul(&result, x)
	result.CyclotomicSquare(&result)
	result.CyclotomicSquare(&result)
	result.CyclotomicSquare(&result)
	result.CyclotomicSquare(&result)
	result.Mul(&result, x)
	result.CyclotomicSquare(&result)
	result.Mul(&result, x)
	result.CyclotomicSquare(&result)
	result.CyclotomicSquare(&result)

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
