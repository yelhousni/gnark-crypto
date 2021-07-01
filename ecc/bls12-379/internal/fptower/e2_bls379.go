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

package fptower

import (
	"github.com/consensys/gnark-crypto/ecc/bls12-379/fp"
)

// Mul sets z to the E2-product of x,y, returns z
func (z *E2) Mul(x, y *E2) *E2 {
	var a, b, c fp.Element
	a.Add(&x.A0, &x.A1)
	b.Add(&y.A0, &y.A1)
	a.Mul(&a, &b)
	b.Mul(&x.A0, &y.A0)
	c.Mul(&x.A1, &y.A1)
	z.A1.Sub(&a, &b).Sub(&z.A1, &c)
	z.A0.Double(&c).Double(&z.A0).Add(&z.A0, &c)
	z.A0.Sub(&b, &z.A0)
	return z
}

// Square sets z to the E2-product of x,x returns z
func (z *E2) Square(x *E2) *E2 {
	//algo 22 https://eprint.iacr.org/2010/354.pdf
	var c0, c2 fp.Element
	c0.Add(&x.A0, &x.A1)
	c2.Double(&x.A1).Double(&c2).Add(&c2, &x.A1).Neg(&c2).Add(&c2, &x.A0)

	c0.Mul(&c0, &c2) // (x1+x2)*(x1+(u**2)x2)
	c2.Mul(&x.A0, &x.A1).Double(&c2)
	z.A1 = c2
	c2.Double(&c2)
	z.A0.Add(&c0, &c2)

	return z
}

// MulByNonResidue multiplies a E2 by (5,1)
func (z *E2) MulByNonResidue(x *E2) *E2 {
	var a, b fp.Element
	c := x.A0
	a.Double(&x.A0).Double(&a).Add(&a, &x.A0)
	b.Double(&x.A1).Double(&b).Add(&b, &x.A1)
	z.A0.Sub(&a, &b)
	z.A1.Add(&b, &c)
	return z
}

// MulByNonResidueInv multiplies a E2 by (5,1)^{-1}
func (z *E2) MulByNonResidueInv(x *E2) *E2 {
	var nonResInv E2
	nonResInv.A0 = fp.Element{
		15001115058802174625,
		8950107502879091364,
		4350176883584445097,
		3470226792022368681,
		368586124710534677,
		144318007828352746,
	}
	nonResInv.A1 = fp.Element{
		8169829964025237780,
		15221973684595632812,
		5549768255880001399,
		5053474046131546590,
		15281643058148911860,
		274254260964319711,
	}
	z.Mul(x, &nonResInv)
	return z
}

// Inverse sets z to the E2-inverse of x, returns z
func (z *E2) Inverse(x *E2) *E2 {
	// Algorithm 8 from https://eprint.iacr.org/2010/354.pdf
	//var a, b, t0, t1, tmp fp.Element
	var t0, t1, tmp fp.Element
	a := &x.A0 // creating the buffers a, b is faster than querying &x.A0, &x.A1 in the functions call below
	b := &x.A1
	t0.Square(a)
	t1.Square(b)
	tmp.Double(&t1).Double(&tmp).Add(&tmp, &t1)
	t0.Add(&t0, &tmp)
	t1.Inverse(&t0)
	z.A0.Mul(a, &t1)
	z.A1.Mul(b, &t1).Neg(&z.A1)

	return z
}

// norm sets x to the norm of z
func (z *E2) norm(x *fp.Element) {
	var tmp fp.Element
	x.Square(&z.A1)
	tmp.Double(x).Double(&tmp).Add(&tmp, x)
	x.Square(&z.A0).Add(x, &tmp)
}
