package fptower

// Expt set z to x^t in E12 and return z
func (z *E12) Expt(x *E12) *E12 {
	// const tAbsVal uint64 = 11170052975785672705
	// tAbsVal in binary: 1001101100000100000000000000000000000000000000000000000000000001
	// drop the low 50 bits (all 0 except the least significant bit): 10011011000001 = 9921
	// Shortest addition chains can be found at https://wwwhomes.uni-bielefeld.de/achim/addition_chain.html
	var result, x4, x5 E12
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
