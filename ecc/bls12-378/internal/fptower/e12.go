// Copyright 2020 ConsenSys Software Inc.
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

// Code generated by consensys/gnark-crypto DO NOT EDIT

package fptower

import (
	"encoding/binary"
	"errors"
	"github.com/consensys/gnark-crypto/ecc"
	"github.com/consensys/gnark-crypto/ecc/bls12-378/fp"
	"github.com/consensys/gnark-crypto/ecc/bls12-378/fr"
	"math/big"
)

// E12 is a degree two finite field extension of fp6
type E12 struct {
	C0, C1 E6
}

// Equal returns true if z equals x, fasle otherwise
func (z *E12) Equal(x *E12) bool {
	return z.C0.Equal(&x.C0) && z.C1.Equal(&x.C1)
}

// String puts E12 in string form
func (z *E12) String() string {
	return (z.C0.String() + "+(" + z.C1.String() + ")*w")
}

// SetString sets a E12 from string
func (z *E12) SetString(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11 string) *E12 {
	z.C0.SetString(s0, s1, s2, s3, s4, s5)
	z.C1.SetString(s6, s7, s8, s9, s10, s11)
	return z
}

// Set copies x into z and returns z
func (z *E12) Set(x *E12) *E12 {
	z.C0 = x.C0
	z.C1 = x.C1
	return z
}

// SetOne sets z to 1 in Montgomery form and returns z
func (z *E12) SetOne() *E12 {
	*z = E12{}
	z.C0.B0.A0.SetOne()
	return z
}

// ToMont converts to Mont form
func (z *E12) ToMont() *E12 {
	z.C0.ToMont()
	z.C1.ToMont()
	return z
}

// FromMont converts from Mont form
func (z *E12) FromMont() *E12 {
	z.C0.FromMont()
	z.C1.FromMont()
	return z
}

// Add set z=x+y in E12 and return z
func (z *E12) Add(x, y *E12) *E12 {
	z.C0.Add(&x.C0, &y.C0)
	z.C1.Add(&x.C1, &y.C1)
	return z
}

// Sub sets z to x sub y and return z
func (z *E12) Sub(x, y *E12) *E12 {
	z.C0.Sub(&x.C0, &y.C0)
	z.C1.Sub(&x.C1, &y.C1)
	return z
}

// Double sets z=2*x and returns z
func (z *E12) Double(x *E12) *E12 {
	z.C0.Double(&x.C0)
	z.C1.Double(&x.C1)
	return z
}

// SetRandom used only in tests
func (z *E12) SetRandom() (*E12, error) {
	if _, err := z.C0.SetRandom(); err != nil {
		return nil, err
	}
	if _, err := z.C1.SetRandom(); err != nil {
		return nil, err
	}
	return z, nil
}

// Mul set z=x*y in E12 and return z
func (z *E12) Mul(x, y *E12) *E12 {
	var a, b, c E6
	a.Add(&x.C0, &x.C1)
	b.Add(&y.C0, &y.C1)
	a.Mul(&a, &b)
	b.Mul(&x.C0, &y.C0)
	c.Mul(&x.C1, &y.C1)
	z.C1.Sub(&a, &b).Sub(&z.C1, &c)
	z.C0.MulByNonResidue(&c).Add(&z.C0, &b)
	return z
}

// Square set z=x*x in E12 and return z
func (z *E12) Square(x *E12) *E12 {

	//Algorithm 22 from https://eprint.iacr.org/2010/354.pdf
	var c0, c2, c3 E6
	c0.Sub(&x.C0, &x.C1)
	c3.MulByNonResidue(&x.C1).Neg(&c3).Add(&x.C0, &c3)
	c2.Mul(&x.C0, &x.C1)
	c0.Mul(&c0, &c3).Add(&c0, &c2)
	z.C1.Double(&c2)
	c2.MulByNonResidue(&c2)
	z.C0.Add(&c0, &c2)

	return z
}

// Karabina's compressed cyclotomic square
// https://eprint.iacr.org/2010/542.pdf
// Th. 3.2 with minor modifications to fit our tower
func (z *E12) CyclotomicSquareCompressed(x *E12) *E12 {

	var t [7]E2

	// t0 = g1^2
	t[0].Square(&x.C0.B1)
	// t1 = g5^2
	t[1].Square(&x.C1.B2)
	// t5 = g1 + g5
	t[5].Add(&x.C0.B1, &x.C1.B2)
	// t2 = (g1 + g5)^2
	t[2].Square(&t[5])

	// t3 = g1^2 + g5^2
	t[3].Add(&t[0], &t[1])
	// t5 = 2 * g1 * g5
	t[5].Sub(&t[2], &t[3])

	// t6 = g3 + g2
	t[6].Add(&x.C1.B0, &x.C0.B2)
	// t3 = (g3 + g2)^2
	t[3].Square(&t[6])
	// t2 = g3^2
	t[2].Square(&x.C1.B0)

	// t6 = 2 * nr * g1 * g5
	t[6].MulByNonResidue(&t[5])
	// t5 = 4 * nr * g1 * g5 + 2 * g3
	t[5].Add(&t[6], &x.C1.B0).
		Double(&t[5])
	// z3 = 6 * nr * g1 * g5 + 2 * g3
	z.C1.B0.Add(&t[5], &t[6])

	// t4 = nr * g5^2
	t[4].MulByNonResidue(&t[1])
	// t5 = nr * g5^2 + g1^2
	t[5].Add(&t[0], &t[4])
	// t6 = nr * g5^2 + g1^2 - g2
	t[6].Sub(&t[5], &x.C0.B2)

	// t1 = g2^2
	t[1].Square(&x.C0.B2)

	// t6 = 2 * nr * g5^2 + 2 * g1^2 - 2*g2
	t[6].Double(&t[6])
	// z2 = 3 * nr * g5^2 + 3 * g1^2 - 2*g2
	z.C0.B2.Add(&t[6], &t[5])

	// t4 = nr * g2^2
	t[4].MulByNonResidue(&t[1])
	// t5 = g3^2 + nr * g2^2
	t[5].Add(&t[2], &t[4])
	// t6 = g3^2 + nr * g2^2 - g1
	t[6].Sub(&t[5], &x.C0.B1)
	// t6 = 2 * g3^2 + 2 * nr * g2^2 - 2 * g1
	t[6].Double(&t[6])
	// z1 = 3 * g3^2 + 3 * nr * g2^2 - 2 * g1
	z.C0.B1.Add(&t[6], &t[5])

	// t0 = g2^2 + g3^2
	t[0].Add(&t[2], &t[1])
	// t5 = 2 * g3 * g2
	t[5].Sub(&t[3], &t[0])
	// t6 = 2 * g3 * g2 + g5
	t[6].Add(&t[5], &x.C1.B2)
	// t6 = 4 * g3 * g2 + 2 * g5
	t[6].Double(&t[6])
	// z5 = 6 * g3 * g2 + 2 * g5
	z.C1.B2.Add(&t[5], &t[6])

	return z
}

// Decompress Karabina's cyclotomic square result
func (z *E12) Decompress(x *E12) *E12 {

	var t [3]E2
	var one E2
	one.SetOne()

	// t0 = g1^2
	t[0].Square(&x.C0.B1)
	// t1 = 3 * g1^2 - 2 * g2
	t[1].Sub(&t[0], &x.C0.B2).
		Double(&t[1]).
		Add(&t[1], &t[0])
		// t0 = E * g5^2 + t1
	t[2].Square(&x.C1.B2)
	t[0].MulByNonResidue(&t[2]).
		Add(&t[0], &t[1])
	// t1 = 1/(4 * g3)
	t[1].Double(&x.C1.B0).
		Double(&t[1]).
		Inverse(&t[1]) // costly
	// z4 = g4
	z.C1.B1.Mul(&t[0], &t[1])

	// t1 = g2 * g1
	t[1].Mul(&x.C0.B2, &x.C0.B1)
	// t2 = 2 * g4^2 - 3 * g2 * g1
	t[2].Square(&z.C1.B1).
		Sub(&t[2], &t[1]).
		Double(&t[2]).
		Sub(&t[2], &t[1])
	// t1 = g3 * g5
	t[1].Mul(&x.C1.B0, &x.C1.B2)
	// c_0 = E * (2 * g4^2 + g3 * g5 - 3 * g2 * g1) + 1
	t[2].Add(&t[2], &t[1])
	z.C0.B0.MulByNonResidue(&t[2]).
		Add(&z.C0.B0, &one)

	z.C0.B1.Set(&x.C0.B1)
	z.C0.B2.Set(&x.C0.B2)
	z.C1.B0.Set(&x.C1.B0)
	z.C1.B2.Set(&x.C1.B2)

	return z
}

// BatchDecompress multiple Karabina's cyclotomic square results
func BatchDecompress(x []E12) []E12 {

	n := len(x)
	if n == 0 {
		return x
	}

	t0 := make([]E2, n)
	t1 := make([]E2, n)
	t2 := make([]E2, n)

	var one E2
	one.SetOne()

	for i := 0; i < n; i++ {
		// t0 = g1^2
		t0[i].Square(&x[i].C0.B1)
		// t1 = 3 * g1^2 - 2 * g2
		t1[i].Sub(&t0[i], &x[i].C0.B2).
			Double(&t1[i]).
			Add(&t1[i], &t0[i])
			// t0 = E * g5^2 + t1
		t2[i].Square(&x[i].C1.B2)
		t0[i].MulByNonResidue(&t2[i]).
			Add(&t0[i], &t1[i])
		// t1 = 4 * g3
		t1[i].Double(&x[i].C1.B0).
			Double(&t1[i])
	}

	t1 = BatchInvert(t1) // costs 1 inverse

	for i := 0; i < n; i++ {
		// z4 = g4
		x[i].C1.B1.Mul(&t0[i], &t1[i])

		// t1 = g2 * g1
		t1[i].Mul(&x[i].C0.B2, &x[i].C0.B1)
		// t2 = 2 * g4^2 - 3 * g2 * g1
		t2[i].Square(&x[i].C1.B1)
		t2[i].Sub(&t2[i], &t1[i])
		t2[i].Double(&t2[i])
		t2[i].Sub(&t2[i], &t1[i])

		// t1 = g3 * g5
		t1[i].Mul(&x[i].C1.B0, &x[i].C1.B2)
		// z0 = E * (2 * g4^2 + g3 * g5 - 3 * g2 * g1) + 1
		t2[i].Add(&t2[i], &t1[i])
		x[i].C0.B0.MulByNonResidue(&t2[i]).
			Add(&x[i].C0.B0, &one)
	}

	return x
}

// Granger-Scott's cyclotomic square
// https://eprint.iacr.org/2009/565.pdf, 3.2
func (z *E12) CyclotomicSquare(x *E12) *E12 {

	// x=(x0,x1,x2,x3,x4,x5,x6,x7) in E2^6
	// cyclosquare(x)=(3*x4^2*u + 3*x0^2 - 2*x0,
	//					3*x2^2*u + 3*x3^2 - 2*x1,
	//					3*x5^2*u + 3*x1^2 - 2*x2,
	//					6*x1*x5*u + 2*x3,
	//					6*x0*x4 + 2*x4,
	//					6*x2*x3 + 2*x5)

	var t [9]E2

	t[0].Square(&x.C1.B1)
	t[1].Square(&x.C0.B0)
	t[6].Add(&x.C1.B1, &x.C0.B0).Square(&t[6]).Sub(&t[6], &t[0]).Sub(&t[6], &t[1]) // 2*x4*x0
	t[2].Square(&x.C0.B2)
	t[3].Square(&x.C1.B0)
	t[7].Add(&x.C0.B2, &x.C1.B0).Square(&t[7]).Sub(&t[7], &t[2]).Sub(&t[7], &t[3]) // 2*x2*x3
	t[4].Square(&x.C1.B2)
	t[5].Square(&x.C0.B1)
	t[8].Add(&x.C1.B2, &x.C0.B1).Square(&t[8]).Sub(&t[8], &t[4]).Sub(&t[8], &t[5]).MulByNonResidue(&t[8]) // 2*x5*x1*u

	t[0].MulByNonResidue(&t[0]).Add(&t[0], &t[1]) // x4^2*u + x0^2
	t[2].MulByNonResidue(&t[2]).Add(&t[2], &t[3]) // x2^2*u + x3^2
	t[4].MulByNonResidue(&t[4]).Add(&t[4], &t[5]) // x5^2*u + x1^2

	z.C0.B0.Sub(&t[0], &x.C0.B0).Double(&z.C0.B0).Add(&z.C0.B0, &t[0])
	z.C0.B1.Sub(&t[2], &x.C0.B1).Double(&z.C0.B1).Add(&z.C0.B1, &t[2])
	z.C0.B2.Sub(&t[4], &x.C0.B2).Double(&z.C0.B2).Add(&z.C0.B2, &t[4])

	z.C1.B0.Add(&t[8], &x.C1.B0).Double(&z.C1.B0).Add(&z.C1.B0, &t[8])
	z.C1.B1.Add(&t[6], &x.C1.B1).Double(&z.C1.B1).Add(&z.C1.B1, &t[6])
	z.C1.B2.Add(&t[7], &x.C1.B2).Double(&z.C1.B2).Add(&z.C1.B2, &t[7])

	return z
}

// Inverse set z to the inverse of x in E12 and return z
func (z *E12) Inverse(x *E12) *E12 {
	// Algorithm 23 from https://eprint.iacr.org/2010/354.pdf

	var t0, t1, tmp E6
	t0.Square(&x.C0)
	t1.Square(&x.C1)
	tmp.MulByNonResidue(&t1)
	t0.Sub(&t0, &tmp)
	t1.Inverse(&t0)
	z.C0.Mul(&x.C0, &t1)
	z.C1.Mul(&x.C1, &t1).Neg(&z.C1)

	return z
}

// Exp sets z=x**e and returns it
// uses 2-bits windowed method
func (z *E12) Exp(x *E12, e big.Int) *E12 {

	var res E12
	var ops [3]E12

	res.SetOne()
	ops[0].Set(x)
	ops[1].Square(&ops[0])
	ops[2].Set(&ops[0]).Mul(&ops[2], &ops[1])

	b := e.Bytes()
	for i := range b {
		w := b[i]
		mask := byte(0xc0)
		for j := 0; j < 4; j++ {
			res.Square(&res).Square(&res)
			c := (w & mask) >> (6 - 2*j)
			if c != 0 {
				res.Mul(&res, &ops[c-1])
			}
			mask = mask >> 2
		}
	}
	z.Set(&res)

	return z
}

// CyclotomicExp sets z=x**e and returns it
// uses 2-NAF decomposition
// x must be in the cyclotomic subgroup
// TODO: use a windowed method
func (z *E12) CyclotomicExp(x *E12, e big.Int) *E12 {
	var res, xInv E12
	xInv.InverseUnitary(x)
	res.SetOne()
	eNAF := make([]int8, e.BitLen()+3)
	n := ecc.NafDecomposition(&e, eNAF[:])
	for i := n - 1; i >= 0; i-- {
		res.CyclotomicSquare(&res)
		if eNAF[i] == 1 {
			res.Mul(&res, x)
		} else if eNAF[i] == -1 {
			res.Mul(&res, &xInv)
		}
	}
	z.Set(&res)
	return z
}

// ExpGLV sets z=x**e and returns it
// uses 2-dimensional GLV with 2-bits windowed method
// x must be in GT
// TODO: use 2-NAF
// TODO: use higher dimensional decomposition
func (p *E12) ExpGLV(a *E12, s *big.Int) *E12 {

	var table [15]E12
	var res E12
	var k1, k2 fr.Element

	res.SetOne()

	// table[b3b2b1b0-1] = b3b2*Frobinius(a) + b1b0*a
	table[0].Set(a)
	table[3].Frobenius(a)

	// split the scalar, modifies +-x, Frob(x) accordingly
	k := ecc.SplitScalar(s, &glvBasis)

	if k[0].Sign() == -1 {
		k[0].Neg(&k[0])
		table[0].InverseUnitary(&table[0])
	}
	if k[1].Sign() == -1 {
		k[1].Neg(&k[1])
		table[3].InverseUnitary(&table[3])
	}

	// precompute table (2 bits sliding window)
	// table[b3b2b1b0-1] = b3b2*Frobenius(a) + b1b0*a if b3b2b1b0 != 0
	table[1].CyclotomicSquare(&table[0])
	table[2].Mul(&table[1], &table[0])
	table[4].Mul(&table[3], &table[0])
	table[5].Mul(&table[3], &table[1])
	table[6].Mul(&table[3], &table[2])
	table[7].CyclotomicSquare(&table[3])
	table[8].Mul(&table[7], &table[0])
	table[9].Mul(&table[7], &table[1])
	table[10].Mul(&table[7], &table[2])
	table[11].Mul(&table[7], &table[3])
	table[12].Mul(&table[11], &table[0])
	table[13].Mul(&table[11], &table[1])
	table[14].Mul(&table[11], &table[2])

	// bounds on the lattice base vectors guarantee that k1, k2 are len(r)/2 bits long max
	k1.SetBigInt(&k[0]).FromMont()
	k2.SetBigInt(&k[1]).FromMont()

	// loop starts from len(k1)/2 due to the bounds
	for i := len(k1) / 2; i >= 0; i-- {
		mask := uint64(3) << 62
		for j := 0; j < 32; j++ {
			res.CyclotomicSquare(&res).CyclotomicSquare(&res)
			b1 := (k1[i] & mask) >> (62 - 2*j)
			b2 := (k2[i] & mask) >> (62 - 2*j)
			if b1|b2 != 0 {
				s := (b2<<2 | b1)
				res.Mul(&res, &table[s-1])
			}
			mask = mask >> 2
		}
	}

	p.Set(&res)
	return p
}

// InverseUnitary inverse a unitary element
func (z *E12) InverseUnitary(x *E12) *E12 {
	return z.Conjugate(x)
}

// Conjugate set z to x conjugated and return z
func (z *E12) Conjugate(x *E12) *E12 {
	*z = *x
	z.C1.Neg(&z.C1)
	return z
}

// SizeOfGT represents the size in bytes that a GT element need in binary form
const SizeOfGT = 48 * 12

// Marshal converts z to a byte slice
func (z *E12) Marshal() []byte {
	b := z.Bytes()
	return b[:]
}

// Unmarshal is an allias to SetBytes()
func (z *E12) Unmarshal(buf []byte) error {
	return z.SetBytes(buf)
}

// Bytes returns the regular (non montgomery) value
// of z as a big-endian byte array.
// z.C1.B2.A1 | z.C1.B2.A0 | z.C1.B1.A1 | ...
func (z *E12) Bytes() (r [SizeOfGT]byte) {
	_z := *z
	_z.FromMont()
	binary.BigEndian.PutUint64(r[568:576], _z.C0.B0.A0[0])
	binary.BigEndian.PutUint64(r[560:568], _z.C0.B0.A0[1])
	binary.BigEndian.PutUint64(r[552:560], _z.C0.B0.A0[2])
	binary.BigEndian.PutUint64(r[544:552], _z.C0.B0.A0[3])
	binary.BigEndian.PutUint64(r[536:544], _z.C0.B0.A0[4])
	binary.BigEndian.PutUint64(r[528:536], _z.C0.B0.A0[5])

	binary.BigEndian.PutUint64(r[520:528], _z.C0.B0.A1[0])
	binary.BigEndian.PutUint64(r[512:520], _z.C0.B0.A1[1])
	binary.BigEndian.PutUint64(r[504:512], _z.C0.B0.A1[2])
	binary.BigEndian.PutUint64(r[496:504], _z.C0.B0.A1[3])
	binary.BigEndian.PutUint64(r[488:496], _z.C0.B0.A1[4])
	binary.BigEndian.PutUint64(r[480:488], _z.C0.B0.A1[5])

	binary.BigEndian.PutUint64(r[472:480], _z.C0.B1.A0[0])
	binary.BigEndian.PutUint64(r[464:472], _z.C0.B1.A0[1])
	binary.BigEndian.PutUint64(r[456:464], _z.C0.B1.A0[2])
	binary.BigEndian.PutUint64(r[448:456], _z.C0.B1.A0[3])
	binary.BigEndian.PutUint64(r[440:448], _z.C0.B1.A0[4])
	binary.BigEndian.PutUint64(r[432:440], _z.C0.B1.A0[5])

	binary.BigEndian.PutUint64(r[424:432], _z.C0.B1.A1[0])
	binary.BigEndian.PutUint64(r[416:424], _z.C0.B1.A1[1])
	binary.BigEndian.PutUint64(r[408:416], _z.C0.B1.A1[2])
	binary.BigEndian.PutUint64(r[400:408], _z.C0.B1.A1[3])
	binary.BigEndian.PutUint64(r[392:400], _z.C0.B1.A1[4])
	binary.BigEndian.PutUint64(r[384:392], _z.C0.B1.A1[5])

	binary.BigEndian.PutUint64(r[376:384], _z.C0.B2.A0[0])
	binary.BigEndian.PutUint64(r[368:376], _z.C0.B2.A0[1])
	binary.BigEndian.PutUint64(r[360:368], _z.C0.B2.A0[2])
	binary.BigEndian.PutUint64(r[352:360], _z.C0.B2.A0[3])
	binary.BigEndian.PutUint64(r[344:352], _z.C0.B2.A0[4])
	binary.BigEndian.PutUint64(r[336:344], _z.C0.B2.A0[5])

	binary.BigEndian.PutUint64(r[328:336], _z.C0.B2.A1[0])
	binary.BigEndian.PutUint64(r[320:328], _z.C0.B2.A1[1])
	binary.BigEndian.PutUint64(r[312:320], _z.C0.B2.A1[2])
	binary.BigEndian.PutUint64(r[304:312], _z.C0.B2.A1[3])
	binary.BigEndian.PutUint64(r[296:304], _z.C0.B2.A1[4])
	binary.BigEndian.PutUint64(r[288:296], _z.C0.B2.A1[5])

	binary.BigEndian.PutUint64(r[280:288], _z.C1.B0.A0[0])
	binary.BigEndian.PutUint64(r[272:280], _z.C1.B0.A0[1])
	binary.BigEndian.PutUint64(r[264:272], _z.C1.B0.A0[2])
	binary.BigEndian.PutUint64(r[256:264], _z.C1.B0.A0[3])
	binary.BigEndian.PutUint64(r[248:256], _z.C1.B0.A0[4])
	binary.BigEndian.PutUint64(r[240:248], _z.C1.B0.A0[5])

	binary.BigEndian.PutUint64(r[232:240], _z.C1.B0.A1[0])
	binary.BigEndian.PutUint64(r[224:232], _z.C1.B0.A1[1])
	binary.BigEndian.PutUint64(r[216:224], _z.C1.B0.A1[2])
	binary.BigEndian.PutUint64(r[208:216], _z.C1.B0.A1[3])
	binary.BigEndian.PutUint64(r[200:208], _z.C1.B0.A1[4])
	binary.BigEndian.PutUint64(r[192:200], _z.C1.B0.A1[5])

	binary.BigEndian.PutUint64(r[184:192], _z.C1.B1.A0[0])
	binary.BigEndian.PutUint64(r[176:184], _z.C1.B1.A0[1])
	binary.BigEndian.PutUint64(r[168:176], _z.C1.B1.A0[2])
	binary.BigEndian.PutUint64(r[160:168], _z.C1.B1.A0[3])
	binary.BigEndian.PutUint64(r[152:160], _z.C1.B1.A0[4])
	binary.BigEndian.PutUint64(r[144:152], _z.C1.B1.A0[5])

	binary.BigEndian.PutUint64(r[136:144], _z.C1.B1.A1[0])
	binary.BigEndian.PutUint64(r[128:136], _z.C1.B1.A1[1])
	binary.BigEndian.PutUint64(r[120:128], _z.C1.B1.A1[2])
	binary.BigEndian.PutUint64(r[112:120], _z.C1.B1.A1[3])
	binary.BigEndian.PutUint64(r[104:112], _z.C1.B1.A1[4])
	binary.BigEndian.PutUint64(r[96:104], _z.C1.B1.A1[5])

	binary.BigEndian.PutUint64(r[88:96], _z.C1.B2.A0[0])
	binary.BigEndian.PutUint64(r[80:88], _z.C1.B2.A0[1])
	binary.BigEndian.PutUint64(r[72:80], _z.C1.B2.A0[2])
	binary.BigEndian.PutUint64(r[64:72], _z.C1.B2.A0[3])
	binary.BigEndian.PutUint64(r[56:64], _z.C1.B2.A0[4])
	binary.BigEndian.PutUint64(r[48:56], _z.C1.B2.A0[5])

	binary.BigEndian.PutUint64(r[40:48], _z.C1.B2.A1[0])
	binary.BigEndian.PutUint64(r[32:40], _z.C1.B2.A1[1])
	binary.BigEndian.PutUint64(r[24:32], _z.C1.B2.A1[2])
	binary.BigEndian.PutUint64(r[16:24], _z.C1.B2.A1[3])
	binary.BigEndian.PutUint64(r[8:16], _z.C1.B2.A1[4])
	binary.BigEndian.PutUint64(r[0:8], _z.C1.B2.A1[5])

	return
}

// SetBytes interprets e as the bytes of a big-endian GT
// sets z to that value (in Montgomery form), and returns z.
// size(e) == 48 * 12
// z.C1.B2.A1 | z.C1.B2.A0 | z.C1.B1.A1 | ...
func (z *E12) SetBytes(e []byte) error {
	if len(e) != SizeOfGT {
		return errors.New("invalid buffer size")
	}
	z.C0.B0.A0.SetBytes(e[528 : 528+fp.Bytes])

	z.C0.B0.A1.SetBytes(e[480 : 480+fp.Bytes])

	z.C0.B1.A0.SetBytes(e[432 : 432+fp.Bytes])

	z.C0.B1.A1.SetBytes(e[384 : 384+fp.Bytes])

	z.C0.B2.A0.SetBytes(e[336 : 336+fp.Bytes])

	z.C0.B2.A1.SetBytes(e[288 : 288+fp.Bytes])

	z.C1.B0.A0.SetBytes(e[240 : 240+fp.Bytes])

	z.C1.B0.A1.SetBytes(e[192 : 192+fp.Bytes])

	z.C1.B1.A0.SetBytes(e[144 : 144+fp.Bytes])

	z.C1.B1.A1.SetBytes(e[96 : 96+fp.Bytes])

	z.C1.B2.A0.SetBytes(e[48 : 48+fp.Bytes])

	z.C1.B2.A1.SetBytes(e[0 : 0+fp.Bytes])

	return nil
}

// IsInSubGroup ensures GT/E12 is in correct sugroup
func (z *E12) IsInSubGroup() bool {
	var a, b E12

	// check z^(Phi_k(p)) == 1
	a.FrobeniusSquare(z)
	b.FrobeniusSquare(&a).Mul(&b, z)

	if !a.Equal(&b) {
		return false
	}

	// check z^(p+1-t) == 1
	a.Frobenius(z)
	b.Expt(z)

	return a.Equal(&b)
}

func (z *E12) Select(cond int, caseZ *E12, caseNz *E12) *E12 {
	//Might be able to save a nanosecond or two by an aggregate implementation

	z.C0.Select(cond, &caseZ.C0, &caseNz.C0)
	z.C1.Select(cond, &caseZ.C1, &caseNz.C1)

	return z
}

func (z *E12) Div(x *E12, y *E12) *E12 {
	var r E12
	r.Inverse(y).Mul(x, &r)
	return z.Set(&r)
}
