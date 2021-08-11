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
	"github.com/consensys/gnark-crypto/ecc/bw6-767/fp"
)

/*
var _frobA = fp.Element{
	9193734820520314185,
	15390913228415833887,
	5309822015742495676,
	5431732283202763350,
	17252325881282386417,
	298854800984767943,
	15252629665615712253,
	11476276919959978448,
	6617989123466214626,
	293279592164056124,
	3271178847573361778,
	76563709148138387,
}
var _frobB = fp.Element{
	7467050525960156664,
	11327349735975181567,
	4886471689715601876,
	825788856423438757,
	532349992164519008,
	5190235139112556877,
	10134108925459365126,
	2188880696701890397,
	14832254987849135908,
	2933451070611009188,
	11385631952165834796,
	64130670718986244,
}
var _frobC = fp.Element{
	10159193990637832851,
	5286779382647858051,
	15149190582698529379,
	10172307932521123666,
	7672315572788794062,
	4504265454330324035,
	8586997380578354686,
	5916374020980521403,
	9559933215456904989,
	10407926721244239843,
	3712625600415690514,
	17752318063289862,
}
var _frobAC = fp.Element{
	17481284903592032950,
	10104133845767975835,
	8607375506753517913,
	13706168424391191299,
	9580010308493592354,
	14241333420363995524,
	6665632285037357566,
	5559902898979457045,
	15504799981718861253,
	8332096944629367896,
	18005297320867222879,
	58811391084848524,
}
var _frobBC = fp.Element{
	8432509696077675330,
	1223215890207205731,
	14725840256671635579,
	5566364505741799073,
	9399083757380478269,
	9395645792458112968,
	3468476640422007559,
	15075721871431984968,
	17774199079839826270,
	13048098199691192907,
	11827078705008163532,
	5319279634137719,
}
*/

var _p = fp.Modulus()

// Frobenius set z in E6 to x^q, return z
func (z *E6) Frobenius(x *E6) *E6 {
	z.Exp(x, *_p)
	/*
		z.B0.A0 = x.B0.A0
		z.B0.A1.Mul(&x.B0.A1, &_frobA)
		z.B0.A2.Mul(&x.B0.A2, &_frobB)

		z.B1.A0.Mul(&x.B1.A0, &_frobC)
		z.B1.A1.Mul(&x.B1.A1, &_frobAC)
		z.B1.A2.Mul(&x.B1.A2, &_frobBC)
	*/

	return z
}

// FrobeniusQuad set z in E6 to x^(q^4), return z
func (z *E6) FrobeniusQuad(x *E6) *E6 {
	z.Exp(x, *_p).Exp(z, *_p).Exp(z, *_p).Exp(z, *_p)
	/*
		z.B0.A0 = x.B0.A0
		z.B0.A1.Mul(&x.B0.A1, &_frobA)
		z.B0.A2.Mul(&x.B0.A2, &_frobB)

		z.B1.A0.Mul(&x.B1.A0, &_frobC)
		z.B1.A1.Mul(&x.B1.A1, &_frobAC)
		z.B1.A2.Mul(&x.B1.A2, &_frobBC)
	*/

	return z
}
