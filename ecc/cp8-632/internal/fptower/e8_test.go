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

func TestE8ReceiverIsOperand(t *testing.T) {

	parameters := gopter.DefaultTestParameters()
	parameters.MinSuccessfulTests = 100

	properties := gopter.NewProperties(parameters)

	genA := GenE8()
	genB := GenE8()

	properties.Property("[CP8-632] Having the receiver as operand (addition) should output the same result", prop.ForAll(
		func(a, b *E8) bool {
			var c, d E8
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
		func(a, b *E8) bool {
			var c, d E8
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
		func(a, b *E8) bool {
			var c, d E8
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
		func(a *E8) bool {
			var b E8
			b.Square(a)
			a.Square(a)
			return a.Equal(&b)
		},
		genA,
	))

	properties.Property("[CP8-632] Having the receiver as operand (double) should output the same result", prop.ForAll(
		func(a *E8) bool {
			var b E8
			b.Double(a)
			a.Double(a)
			return a.Equal(&b)
		},
		genA,
	))

	properties.Property("[CP8-632] Having the receiver as operand (Inverse) should output the same result", prop.ForAll(
		func(a *E8) bool {
			var b E8
			b.Inverse(a)
			a.Inverse(a)
			return a.Equal(&b)
		},
		genA,
	))

	properties.Property("[CP8-632] Having the receiver as operand (Conjugate) should output the same result", prop.ForAll(
		func(a *E8) bool {
			var b E8
			b.Conjugate(a)
			a.Conjugate(a)
			return a.Equal(&b)
		},
		genA,
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestE8Ops(t *testing.T) {

	parameters := gopter.DefaultTestParameters()
	parameters.MinSuccessfulTests = 100

	properties := gopter.NewProperties(parameters)

	genA := GenE8()
	genB := GenE8()

	properties.Property("[CP8-632] sub & add should leave an element invariant", prop.ForAll(
		func(a, b *E8) bool {
			var c E8
			c.Set(a)
			c.Add(&c, b).Sub(&c, b)
			return c.Equal(a)
		},
		genA,
		genB,
	))

	properties.Property("[CP8-632] mul & inverse should leave an element invariant", prop.ForAll(
		func(a, b *E8) bool {
			var c, d E8
			d.Inverse(b)
			c.Set(a)
			c.Mul(&c, b).Mul(&c, &d)
			return c.Equal(a)
		},
		genA,
		genB,
	))

	properties.Property("[CP8-632] inverse twice should leave an element invariant", prop.ForAll(
		func(a *E8) bool {
			var b E8
			b.Inverse(a).Inverse(&b)
			return a.Equal(&b)
		},
		genA,
	))

	properties.Property("[CP8-632] square and mul should output the same result", prop.ForAll(
		func(a *E8) bool {
			var b, c E8
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

func BenchmarkE8Add(b *testing.B) {
	var a, c E8
	a.SetRandom()
	c.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.Add(&a, &c)
	}
}

func BenchmarkE8Sub(b *testing.B) {
	var a, c E8
	a.SetRandom()
	c.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.Sub(&a, &c)
	}
}

func BenchmarkE8Mul(b *testing.B) {
	var a, c E8
	a.SetRandom()
	c.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.Mul(&a, &c)
	}
}

func BenchmarkE8Square(b *testing.B) {
	var a E8
	a.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.Square(&a)
	}
}

func BenchmarkE8Inverse(b *testing.B) {
	var a E8
	a.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.Inverse(&a)
	}
}

func BenchmarkE8Conjugate(b *testing.B) {
	var a E8
	a.SetRandom()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.Conjugate(&a)
	}
}
