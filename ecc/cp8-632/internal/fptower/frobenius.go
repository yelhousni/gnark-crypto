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
	"github.com/consensys/gnark-crypto/ecc/cp8-632/fp"
)

// Frobenius
func (z *E8) Frobenius(x *E8) *E8 {
	var b E2

	a := fp.Element{
		3916458128034199340,
		6882541448827412113,
		11224859209293672029,
		1335153420971272167,
		8308071011107334982,
		17446808896553813422,
		10352293250830896992,
		11713874050134427703,
		722244504604808771,
		8930466954225936,
	}
	b.A1 = fp.Element{
		4636478738010456435,
		6326228536455595272,
		1648115418671854160,
		14735586135771377191,
		2863753510355726652,
		16627248998165910150,
		14393081852000620781,
		14420564855755218259,
		12729952670188154943,
		26345310115368755,
	}

	z.C0.B0.Conjugate(&x.C0.B0)
	z.C0.B1.Conjugate(&x.C0.B1).MulByElement(&z.C0.B1, &a)
	z.C1.B0.Conjugate(&x.C1.B0)
	z.C1.B1.Conjugate(&x.C1.B1).MulByElement(&z.C1.B1, &a)

	z.C1.B0.Mul(&z.C1.B0, &b)
	z.C1.B1.Mul(&z.C1.B1, &b)

	return z
}

// FrobeniusSquare
func (z *E8) FrobeniusSquare(x *E8) *E8 {

	a := fp.Element{
		5970287765817275737,
		8878302867479338423,
		16611264740231281568,
		10311898289500844112,
		8371877573987207521,
		9959796394303452291,
		9187638892507884374,
		11930173868072381192,
		12303819203974681632,
		59228732751081619,
	}

	z.C0.Conjugate(&x.C0)
	z.C1.Conjugate(&x.C1).MulByElement(&z.C1, &a)

	return z
}

// FrobeniusCube
func (z *E8) FrobeniusCube(x *E8) *E8 {

	var b E2

	a := fp.Element{
		5970287765817275737,
		8878302867479338423,
		16611264740231281568,
		10311898289500844112,
		8371877573987207521,
		9959796394303452291,
		9187638892507884374,
		11930173868072381192,
		12303819203974681632,
		59228732751081619,
	}
	b.A1 = fp.Element{
		720020609976257095,
		17890431161337734775,
		8870000283087733746,
		13400432714800105023,
		13002426572957943286,
		17627184175321648343,
		4040788601169723788,
		2706690805620790556,
		12007708165583346172,
		17414843161142819,
	}

	z.C0.B0.Conjugate(&x.C0.B0)
	z.C0.B1.Conjugate(&x.C0.B1).MulByElement(&z.C0.B1, &a)
	z.C1.B0.Conjugate(&x.C1.B0)
	z.C1.B1.Conjugate(&x.C1.B1).MulByElement(&z.C1.B1, &a)

	z.C1.B0.Mul(&z.C1.B0, &b)
	z.C1.B1.Mul(&z.C1.B1, &b)

	return z
}
