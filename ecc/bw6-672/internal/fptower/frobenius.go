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
	"github.com/consensys/gnark-crypto/ecc/bw6-672/fp"
)

var _p = fp.Modulus()

// naive Frobenius set z in E6 to Frobenius(x), return z
func (z *E6) Frobenius(x *E6) *E6 {
	z.Exp(x, *_p)
	return z
}

/*
// Frobenius set z in E6 to Frobenius(x), return z
func (z *E6) Frobenius(x *E6) *E6 {
	a := fp.Element{
		13379576892826363940,
		7215010732331695319,
		9587968273680355194,
		17498442511873681277,
		17430577930908868575,
		3900545948475841198,
		4911233234059649485,
		5259203007663136357,
		784833376601091838,
		14101604511763105,
	}
	b := fp.Element{
		14915121275341857583,
		17181477577419481383,
		7084037479344459509,
		15839663182812461911,
		1480269266224751652,
		8132863462598737762,
		18142219247247109524,
		16952785282480028983,
		10364476017837376831,
		17548851641601399,
	}
	c := fp.Element{
		597834311555652830,
		5676150712176383509,
		8459519236066800431,
		11690428270517348528,
		11839864809966557220,
		1185830464157066542,
		5950198841798077595,
		13670804634510857615,
		7801381657215673717,
		65313904097694188,
	}
	ac := fp.Element{
		9847954094458669907,
		5949744236041625087,
		16672005753024814704,
		14891361620976591572,
		464103123424068612,
		12033409411074578961,
		4606708407597207393,
		3765244216433613725,
		11149309394438468670,
		31650456153364504,
	}
	bc := fp.Element{
		2133378694071146473,
		15642617557264169573,
		5955588441730904746,
		10031648941456129162,
		14336300218991991913,
		5418147978279963105,
		734440781275986018,
		6917642835618198626,
		17381024298451958711,
		68761151227532482,
	}

	z.B0.A0.Set(&x.B0.A0)
	z.B0.A1.Mul(&x.B0.A1, &a)
	z.B0.A2.Mul(&x.B0.A2, &b)

	z.B1.A0.Mul(&x.B1.A0, &c)
	z.B1.A1.Mul(&x.B1.A1, &ac)
	z.B1.A2.Mul(&x.B1.A2, &bc)

	return z
}
*/
