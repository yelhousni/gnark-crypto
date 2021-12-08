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

// MulBy023 multiplication by sparse element (c0,0,c2,c3)
func (z *E8) MulBy023(c0, c2, c3 *E2) *E8 {

	var c0z0, c2z2, c3z3, c0z1, c2z1, c3z1, z0, z1, z2, z3, A, B, C, t E2

	c0z0.Mul(c0, &z.C0.B0)
	c2z2.Mul(c2, &z.C1.B0)
	c3z3.Mul(c3, &z.C1.B1)
	c0z1.Mul(c0, &z.C0.B1)
	c2z1.Mul(c2, &z.C0.B1)
	c3z1.Mul(c3, &z.C0.B1)

	t.Add(c0, c2)
	A.Add(&z.C0.B0, &z.C1.B0).
		Mul(&A, &t).
		Sub(&A, &c0z0).
		Sub(&A, &c2z2) // A = c0f2 + c2f0

	t.Add(c0, c3)
	B.Add(&z.C0.B0, &z.C1.B1).
		Mul(&B, &t).
		Sub(&B, &c0z0).
		Sub(&B, &c3z3) // B = c0f3 + c3f0

	t.Add(c2, c3)
	C.Add(&z.C1.B0, &z.C1.B1).
		Mul(&C, &t).
		Sub(&C, &c2z2).
		Sub(&C, &c3z3) // C = c2f3 + c3f2

	z0.MulByNonResidue(&C).Add(&z0, &c0z0)
	z1.MulByNonResidue(&c3z3).Add(&z1, &c2z2).Add(&z1, &c0z1)
	z2.MulByNonResidue(&c3z1).Add(&z2, &A)
	z3.Add(&c2z1, &B)

	z.C0.B0.Set(&z0)
	z.C0.B1.Set(&z1)
	z.C1.B0.Set(&z2)
	z.C1.B1.Set(&z3)

	return z
}
