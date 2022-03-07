/*
Copyright © 2020 ConsenSys

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

package twistededwards

import (
	"math/big"

	"github.com/consensys/gnark-crypto/ecc/bls12-381/fr"
)

// CurveParams curve parameters: ax^2 + y^2 = 1 + d*x^2*y^2
type CurveParams struct {
	A, D     fr.Element // in Montgomery form
	Cofactor fr.Element // not in Montgomery form
	Order    big.Int
	Base     PointAffine
}

var edwards CurveParams

// GetEdwardsCurve returns the twisted Edwards curve on BLS12-381's Fr
func GetEdwardsCurve() CurveParams {

	// copy to keep Order private
	var res CurveParams

	res.A.Set(&edwards.A)
	res.D.Set(&edwards.D)
	res.Cofactor.Set(&edwards.Cofactor)
	res.Order.Set(&edwards.Order)
	res.Base.Set(&edwards.Base)

	return res
}

func init() {

	// Jubjub curve
	edwards.A.SetOne().Neg(&edwards.A)
	edwards.D.SetString("19257038036680949359750312669786877991949435402254120286184196891950884077233") // -(10240/10241)
	edwards.Cofactor.SetUint64(8).ToMont()
	edwards.Order.SetString("6554484396890773809930967563523245729705921265872317281365359162392183254199", 10)

	edwards.Base.X.SetString("23426137002068529236790192115758361610982344002369094106619281483467893291614")
	edwards.Base.Y.SetString("39325435222430376843701388596190331198052476467368316772266670064146548432123")

}

// mulByA multiplies fr.Element by edwards.A
func mulByA(x *fr.Element) {
	x.Neg(x)
}
