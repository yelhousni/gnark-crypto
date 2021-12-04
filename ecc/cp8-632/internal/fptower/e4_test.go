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
	"testing"

	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"
)

// ------------------------------------------------------------
// tests

func TestE4ReceiverIsOperand(t *testing.T) {

	parameters := gopter.DefaultTestParameters()
	parameters.MinSuccessfulTests = 100

	properties := gopter.NewProperties(parameters)

	genA := GenE4()
	genB := GenE4()

	properties.Property("[CP8-632] Having the receiver as operand (addition) should output the same result", prop.ForAll(
		func(a, b *E4) bool {
			var c, d E4
			d.Set(a)
			c.Add(a, b)
			a.Add(a, b)
			b.Add(&d, b)
			return a.Equal(b) && a.Equal(&c) && b.Equal(&c)
		},
		genA,
		genB,
	))

	properties.Property("[CP8-632] Having the receiver as operand (sub) should output the same result", prop.ForAll(
		func(a, b *E4) bool {
			var c, d E4
			d.Set(a)
			c.Sub(a, b)
			a.Sub(a, b)
			b.Sub(&d, b)
			return a.Equal(b) && a.Equal(&c) && b.Equal(&c)
		},
		genA,
		genB,
	))

	properties.Property("[CP8-632] Having the receiver as operand (mul) should output the same result", prop.ForAll(
		func(a, b *E4) bool {
			var c, d E4
			d.Set(a)
			c.Mul(a, b)
			a.Mul(a, b)
			b.Mul(&d, b)
			return a.Equal(b) && a.Equal(&c) && b.Equal(&c)
		},
		genA,
		genB,
	))

	properties.Property("[CP8-632] Having the receiver as operand (square) should output the same result", prop.ForAll(
		func(a *E4) bool {
			var b E4
			b.Square(a)
			a.Square(a)
			return a.Equal(&b)
		},
		genA,
	))

	properties.Property("[CP8-632] Having the receiver as operand (double) should output the same result", prop.ForAll(
		func(a *E4) bool {
			var b E4
			b.Double(a)
			a.Double(a)
			return a.Equal(&b)
		},
		genA,
	))

	properties.Property("[CP8-632] Having the receiver as operand (mul by non residue) should output the same result", prop.ForAll(
		func(a *E4) bool {
			var b E4
			b.MulByNonResidue(a)
			a.MulByNonResidue(a)
			return a.Equal(&b)
		},
		genA,
	))

	properties.Property("[CP8-632] Having the receiver as operand (Inverse) should output the same result", prop.ForAll(
		func(a *E4) bool {
			var b E4
			b.Inverse(a)
			a.Inverse(a)
			return a.Equal(&b)
		},
		genA,
	))

	properties.Property("[CP8-632] Having the receiver as operand (Conjugate) should output the same result", prop.ForAll(
		func(a *E4) bool {
			var b E4
			b.Conjugate(a)
			a.Conjugate(a)
			return a.Equal(&b)
		},
		genA,
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestE4Ops(t *testing.T) {

	parameters := gopter.DefaultTestParameters()
	parameters.MinSuccessfulTests = 100

	properties := gopter.NewProperties(parameters)

	genA := GenE4()
	genB := GenE4()

	properties.Property("[CP8-632] sub & add should leave an element invariant", prop.ForAll(
		func(a, b *E4) bool {
			var c E4
			c.Set(a)
			c.Add(&c, b).Sub(&c, b)
			return c.Equal(a)
		},
		genA,
		genB,
	))

	properties.Property("[CP8-632] mul & inverse should leave an element invariant", prop.ForAll(
		func(a, b *E4) bool {
			var c, d E4
			d.Inverse(b)
			c.Set(a)
			c.Mul(&c, b).Mul(&c, &d)
			return c.Equal(a)
		},
		genA,
		genB,
	))

	properties.Property("[CP8-632] inverse twice should leave an element invariant", prop.ForAll(
		func(a *E4) bool {
			var b E4
			b.Inverse(a).Inverse(&b)
			return a.Equal(&b)
		},
		genA,
	))

	properties.Property("[CP8-632] square and mul should output the same result", prop.ForAll(
		func(a *E4) bool {
			var b, c E4
			b.Mul(a, a)
			c.Square(a)
			return b.Equal(&c)
		},
		genA,
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

// ------------------------------------------------------------
// benches

func BenchmarkE4Add(b *testing.B) {
	var a, c E4
	a.SetRandom()
	c.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.Add(&a, &c)
	}
}

func BenchmarkE4Sub(b *testing.B) {
	var a, c E4
	a.SetRandom()
	c.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.Sub(&a, &c)
	}
}

func BenchmarkE4Mul(b *testing.B) {
	var a, c E4
	a.SetRandom()
	c.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.Mul(&a, &c)
	}
}

func BenchmarkE4Square(b *testing.B) {
	var a E4
	a.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.Square(&a)
	}
}

func BenchmarkE4Inverse(b *testing.B) {
	var a E4
	a.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.Inverse(&a)
	}
}

func BenchmarkE4MulNonRes(b *testing.B) {
	var a E4
	a.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.MulByNonResidue(&a)
	}
}

func BenchmarkE4Conjugate(b *testing.B) {
	var a E4
	a.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.Conjugate(&a)
	}
}
