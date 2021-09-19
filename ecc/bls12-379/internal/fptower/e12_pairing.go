package fptower

func (z *E12) nSquare(n int) {
	for i := 0; i < n; i++ {
		z.CyclotomicSquare(z)
	}
}

func (z *E12) nSquareCompressed(n int) {
	for i := 0; i < n; i++ {
		z.CyclotomicSquareCompressed(z)
	}
}

// Expt set z to x^t in E12 and return z
func (z *E12) Expt(x *E12) *E12 {
	// const tAbsVal uint64 = 11170052975785672705
	// tAbsVal in binary: 1001101100000100000000000000000000000000000000000000000000000001
	// drop the low 50 bits (all 0 except the least significant bit): 10011011000001 = 9921
	// Shortest addition chains can be found at https://wwwhomes.uni-bielefeld.de/achim/addition_chain.html
	var result, x4, x5 E12
	result.Set(x)
	result.nSquare(2)
	x4.Set(&result)
	result.Mul(&result, x)
	x5.Set(&result)
	result.Mul(&result, &x4)
	result.nSquare(3)
	result.Mul(&result, &x5)
	result.CyclotomicSquare(&result)
	result.Mul(&result, x)
	result.nSquare(6)
	result.Mul(&result, x)

	// the remaining 50 bits
    result.nSquareCompressed(50)
    result.Decompress(&result)
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
