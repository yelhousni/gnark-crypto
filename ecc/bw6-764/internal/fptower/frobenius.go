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
	"github.com/consensys/gnark-crypto/ecc/bw6-764/fp"
)

var _frobA = fp.Element{
    16439153956121079689,
    17050947565551831201,
    6541430248152292527,
    15987691203423007698,
    4723498017856209744,
    11008012764507170653,
    2420767866770231066,
    17319827534313390214,
    18288161404050887264,
    12969124009552547627,
    7884759346639203638,
    198874968445229914,
}
var _frobB = fp.Element{
    7920816428325936234,
    9354011573501871825,
    11221799747455232832,
    14503687928435433160,
    5614672112651322354,
    14856992403536051690,
    1807318283793255090,
    4542564508540506816,
    4623588615668520750,
    18059820848367132930,
    1897163546336294903,
    27582194759154967,
}
var _frobC = fp.Element{
    9980992090471785499,
    15855000659121617370,
    14413517482136229143,
    18101285897935477877,
    10875460565813832619,
    10328543569478286171,
    1584640952963340240,
    7636427190363302734,
    12413221449284107090,
    16172590102937706061,
    9756590363442366885,
    784296119889277054,
}
var _frobAC = fp.Element{
    5913226310737464307,
    7958215065344151411,
    17763229995607525360,
    12044635058148889242,
    10338170130507532099,
    7418261094333670727,
    4228086150563486157,
    3415647969144345414,
    4465005946009856399,
    12582200784210128942,
    9781922892975498542,
    226457163204384881,
}
var _frobBC = fp.Element{
    1462654562676642044,
    8158064667071657994,
    647142907729617832,
    16617282622947903340,
    11766634660608945229,
    14177523208507167208,
    971191369986364264,
    13305908238299970952,
    17195392734611292191,
    2816542868042739747,
    3768994563139458151,
    613003346203202107,
}

// Frobenius set z in E6 to Frobenius(x), return z
func (z *E6) Frobenius(x *E6) *E6 {

	z.B0.A0 = x.B0.A0
	z.B0.A1.Mul(&x.B0.A1, &_frobA)
	z.B0.A2.Mul(&x.B0.A2, &_frobB)

	z.B1.A0.Mul(&x.B1.A0, &_frobC)
	z.B1.A1.Mul(&x.B1.A1, &_frobAC)
	z.B1.A2.Mul(&x.B1.A2, &_frobBC)

	return z
}