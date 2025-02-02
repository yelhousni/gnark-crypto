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

package fptower

import (
	"math/big"

	"github.com/consensys/gnark-crypto/ecc/cp8-632/fp"
)

// E2 is a degree two finite field extension of fp.Element
type E2 struct {
	A0, A1 fp.Element
}

// Equal returns true if z equals x, fasle otherwise
func (z *E2) Equal(x *E2) bool {
	return z.A0.Equal(&x.A0) && z.A1.Equal(&x.A1)
}

// Cmp compares (lexicographic order) z and x and returns:
//
//   -1 if z <  x
//    0 if z == x
//   +1 if z >  x
//
func (z *E2) Cmp(x *E2) int {
	if a1 := z.A1.Cmp(&x.A1); a1 != 0 {
		return a1
	}
	return z.A0.Cmp(&x.A0)
}

// LexicographicallyLargest returns true if this element is strictly lexicographically
// larger than its negation, false otherwise
func (z *E2) LexicographicallyLargest() bool {
	// adapted from github.com/zkcrypto/bls12_381
	if z.A1.IsZero() {
		return z.A0.LexicographicallyLargest()
	}
	return z.A1.LexicographicallyLargest()
}

// SetString sets a E2 element from strings
func (z *E2) SetString(s1, s2 string) *E2 {
	z.A0.SetString(s1)
	z.A1.SetString(s2)
	return z
}

// SetZero sets an E2 elmt to zero
func (z *E2) SetZero() *E2 {
	z.A0.SetZero()
	z.A1.SetZero()
	return z
}

// Set sets an E2 from x
func (z *E2) Set(x *E2) *E2 {
	z.A0 = x.A0
	z.A1 = x.A1
	return z
}

// SetOne sets z to 1 in Montgomery form and returns z
func (z *E2) SetOne() *E2 {
	z.A0.SetOne()
	z.A1.SetZero()
	return z
}

// SetRandom sets a0 and a1 to random values
func (z *E2) SetRandom() (*E2, error) {
	if _, err := z.A0.SetRandom(); err != nil {
		return nil, err
	}
	if _, err := z.A1.SetRandom(); err != nil {
		return nil, err
	}
	return z, nil
}

// IsZero returns true if the two elements are equal, fasle otherwise
func (z *E2) IsZero() bool {
	return z.A0.IsZero() && z.A1.IsZero()
}

// Add adds two elements of E2
func (z *E2) Add(x, y *E2) *E2 {
	z.A0.Add(&x.A0, &y.A0)
	z.A1.Add(&x.A1, &y.A1)
	return z
}

// Sub two elements of E2
func (z *E2) Sub(x, y *E2) *E2 {
	z.A0.Sub(&x.A0, &y.A0)
	z.A1.Sub(&x.A1, &y.A1)
	return z
}

// Double doubles an E2 element
func (z *E2) Double(x *E2) *E2 {
	z.A0.Double(&x.A0)
	z.A1.Double(&x.A1)
	return z
}

// Neg negates an E2 element
func (z *E2) Neg(x *E2) *E2 {
	z.A0.Neg(&x.A0)
	z.A1.Neg(&x.A1)
	return z
}

// String implements Stringer interface for fancy printing
func (z *E2) String() string {
	return (z.A0.String() + "+" + z.A1.String() + "*u")
}

// ToMont converts to mont form
func (z *E2) ToMont() *E2 {
	z.A0.ToMont()
	z.A1.ToMont()
	return z
}

// FromMont converts from mont form
func (z *E2) FromMont() *E2 {
	z.A0.FromMont()
	z.A1.FromMont()
	return z
}

// MulByElement multiplies an element in E2 by an element in fp
func (z *E2) MulByElement(x *E2, y *fp.Element) *E2 {
	var yCopy fp.Element
	yCopy.Set(y)
	z.A0.Mul(&x.A0, &yCopy)
	z.A1.Mul(&x.A1, &yCopy)
	return z
}

// Conjugate conjugates an element in E2
func (z *E2) Conjugate(x *E2) *E2 {
	z.A0 = x.A0
	z.A1.Neg(&x.A1)
	return z
}

// Exp sets z=x**e and returns it
func (z *E2) Exp(x E2, exponent *big.Int) *E2 {
	z.SetOne()
	b := exponent.Bytes()
	for i := 0; i < len(b); i++ {
		w := b[i]
		for j := 0; j < 8; j++ {
			z.Square(z)
			if (w & (0b10000000 >> j)) != 0 {
				z.Mul(z, &x)
			}
		}
	}

	return z
}

// Mul set z=x*y in E2 and return z
func (z *E2) Mul(x, y *E2) *E2 {
	var a, b, c fp.Element
	a.Add(&x.A0, &x.A1)
	b.Add(&y.A0, &y.A1)
	a.Mul(&a, &b)
	b.Mul(&x.A0, &y.A0)
	c.Mul(&x.A1, &y.A1)
	z.A1.Sub(&a, &b).Sub(&z.A1, &c)
	z.A0.Double(&c).Add(&z.A0, &b)
	return z
}

// Square set z=x*x in E2 and return z
func (z *E2) Square(x *E2) *E2 {

	//Algorithm 22 from https://eprint.iacr.org/2010/354.pdf
	var c0, c2, c3 fp.Element
	c0.Sub(&x.A0, &x.A1)
	c3.Double(&x.A1).Sub(&x.A0, &c3)
	c2.Mul(&x.A0, &x.A1)
	c0.Mul(&c0, &c3).Add(&c0, &c2)
	z.A1.Double(&c2)
	c2.Double(&c2)
	z.A0.Add(&c0, &c2)

	return z
}

// Inverse set z to the inverse of x in E2 and return z
func (z *E2) Inverse(x *E2) *E2 {
	// Algorithm 23 from https://eprint.iacr.org/2010/354.pdf

	var t0, t1, tmp fp.Element
	t0.Square(&x.A0)
	t1.Square(&x.A1)
	tmp.Double(&t1)
	t0.Sub(&t0, &tmp)
	t1.Inverse(&t0)
	z.A0.Mul(&x.A0, &t1)
	z.A1.Mul(&x.A1, &t1).Neg(&z.A1)

	return z
}

// MulByNonResidue mul x by (0,1)
func (z *E2) MulByNonResidue(x *E2) *E2 {
	z.A1, z.A0 = x.A0, x.A1
	z.A0.Double(&z.A0)
	return z
}

// norm sets x to the norm of z
func (z *E2) norm(x *fp.Element) {
	var tmp fp.Element
	tmp.Square(&z.A1).Double(&tmp)
	x.Square(&z.A0).Sub(x, &tmp)
}

// Legendre returns the Legendre symbol of z
func (z *E2) Legendre() int {
	var n fp.Element
	z.norm(&n)
	return n.Legendre()
}

// Sqrt sets z to the square root of and returns z
// The function does not test wether the square root
// exists or not, it's up to the caller to call
// Legendre beforehand.
// cf https://eprint.iacr.org/2012/685.pdf (algo 10)
func (z *E2) Sqrt(x *E2) *E2 {

	// precomputation
	var b, c, d, e, f, x0 E2
	var _b, o fp.Element

	// c must be a non square (works for p=1 mod 12 hence 1 mod 4)
	c.A1.SetOne()

	q := fp.Modulus()
	var exp, one big.Int
	one.SetUint64(1)
	exp.Set(q).Sub(&exp, &one).Rsh(&exp, 1)
	d.Exp(c, &exp)
	e.Mul(&d, &c).Inverse(&e)
	f.Mul(&d, &c).Square(&f)

	// computation
	exp.Rsh(&exp, 1)
	b.Exp(*x, &exp)
	b.norm(&_b)
	o.SetOne()
	if _b.Equal(&o) {
		x0.Square(&b).Mul(&x0, x)
		_b.Set(&x0.A0).Sqrt(&_b)
		z.Conjugate(&b).MulByElement(z, &_b)
		return z
	}
	x0.Square(&b).Mul(&x0, x).Mul(&x0, &f)
	_b.Set(&x0.A0).Sqrt(&_b)
	z.Conjugate(&b).MulByElement(z, &_b).Mul(z, &e)

	return z
}
