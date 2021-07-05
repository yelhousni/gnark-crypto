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

package fp

// MulByNonResidue multiplies a fp.Element by 2
func (z *Element) MulByNonResidue(x *Element) *Element {
	z.Double(x)
	return z
}

// MulByNonResidueInv multiplies a fp.Element by 2**-1
func (z *Element) MulByNonResidueInv(x *Element) *Element {

	nrInv := Element{
		5994291104030128713,
		17848770620494444700,
		13159415653846744115,
		10280169384111010897,
		3075981273978811437,
		18107009476195109375,
		18028680616806106202,
		13605043901734507875,
		15509274096326161528,
		10825105083547355024,
		935915508401581623,
		292710575722023570,
	}

	z.Mul(x, &nrInv)

	return z
}
