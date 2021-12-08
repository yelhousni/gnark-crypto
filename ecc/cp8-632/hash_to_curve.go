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

package cp8632

import (
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/consensys/gnark-crypto/ecc/cp8-632/fp"
	"github.com/consensys/gnark-crypto/ecc/cp8-632/internal/fptower"
)

// hashToFp hashes msg to count prime field elements.
// https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-06#section-5.2
func hashToFp(msg, dst []byte, count int) ([]fp.Element, error) {

	// 128 bits of security
	// L = ceil((ceil(log2(p)) + k) / 8), where k is the security parameter = 128
	L := 64

	lenInBytes := count * L
	pseudoRandomBytes, err := ecc.ExpandMsgXmd(msg, dst, lenInBytes)
	if err != nil {
		return nil, err
	}

	res := make([]fp.Element, count)
	for i := 0; i < count; i++ {
		res[i].SetBytes(pseudoRandomBytes[i*L : (i+1)*L])
	}
	return res, nil
}

// returns false if u>-u when seen as a bigInt
func sign0(u fp.Element) bool {
	var a, b big.Int
	u.ToBigIntRegular(&a)
	u.Neg(&u)
	u.ToBigIntRegular(&b)
	return a.Cmp(&b) <= 0
}

// ----------------------------------------------------------------------------------------
// G1Affine

// https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-06#section-4.1
// Shallue and van de Woestijne method, works for any elliptic curve in Weierstrass curve
func svdwMapG1(u fp.Element) G1Affine {

	var res G1Affine

	// constants
	// sage script to find z: https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-06#appendix-E.1
	var z, c1, c2, c3, c4 fp.Element
	z.SetUint64(2)
	c1.SetUint64(6)
	c2.SetString("16857842227199225999646786835636637546980511455593181936643515867446449681250963223102065881329922732942164707714798074010252040453159496268796306973476256192757039177873952554508334906059908") // -1
	c3.SetString("1161522585421188978375489229854712031308201942829231321254239305602987112264953234819185572023450683617326826082493819016543191490801004093403702322811785940290484778391978355671865609560986")
	c4.SetString("16857842227199225999646786835636637546980511455593181936643515867446449681250963223102065881329922732942164707714798074010252040453159496268796306973476256192757039177873952554508334906059906") // -3

	var tv1, tv2, tv3, tv4, one, x1, gx1, x2, gx2, x3, x, gx, y fp.Element
	one.SetOne()
	tv1.Square(&u).Mul(&tv1, &c1)
	tv2.Add(&one, &tv1)
	tv1.Sub(&one, &tv1)
	tv3.Mul(&tv2, &tv1).Inverse(&tv3)
	tv4.Mul(&u, &tv1)
	tv4.Mul(&tv4, &tv3)
	tv4.Mul(&tv4, &c3)
	x1.Sub(&c2, &tv4)
	gx1.Square(&x1)
	gx1.Add(&gx1, &aCurveCoeff)
	gx1.Mul(&gx1, &x1)
	e1 := gx1.Legendre()
	x2.Add(&c2, &tv4)
	gx2.Square(&x2)
	gx2.Add(&gx2, &aCurveCoeff)
	gx2.Mul(&gx2, &x2)
	E2 := gx2.Legendre() - e1 // 2 if is_square(gx2) AND NOT e1
	x3.Square(&tv2)
	x3.Mul(&x3, &tv3)
	x3.Square(&x3)
	x3.Mul(&x3, &c4)
	x3.Add(&x3, &z)
	if e1 == 1 {
		x.Set(&x1)
	} else {
		x.Set(&x3)
	}
	if E2 == 2 {
		x.Set(&x2)
	}
	gx.Square(&x)
	gx.Add(&gx, &aCurveCoeff)
	gx.Mul(&gx, &x)
	y.Sqrt(&gx)
	e3 := sign0(u) && sign0(y)
	if !e3 {
		y.Neg(&y)
	}
	res.X.Set(&x)
	res.Y.Set(&y)

	return res
}

// MapToCurveG1Svdw maps an fp.Element to a point on the curve using the Shallue and van de Woestijne map
// https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-06#section-2.2.1
func MapToCurveG1Svdw(t fp.Element) G1Affine {
	res := svdwMapG1(t)
	res.ClearCofactor(&res)
	return res
}

// EncodeToCurveG1Svdw maps an fp.Element to a point on the curve using the Shallue and van de Woestijne map
// https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-06#section-2.2.2
func EncodeToCurveG1Svdw(msg, dst []byte) (G1Affine, error) {
	var res G1Affine
	t, err := hashToFp(msg, dst, 1)
	if err != nil {
		return res, err
	}
	res = MapToCurveG1Svdw(t[0])
	return res, nil
}

// HashToCurveG1Svdw maps an fp.Element to a point on the curve using the Shallue and van de Woestijne map
// https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-06#section-3
func HashToCurveG1Svdw(msg, dst []byte) (G1Affine, error) {
	var res G1Affine
	u, err := hashToFp(msg, dst, 2)
	if err != nil {
		return res, err
	}
	Q0 := MapToCurveG1Svdw(u[0])
	Q1 := MapToCurveG1Svdw(u[1])
	var _Q0, _Q1, _res G1Jac
	_Q0.FromAffine(&Q0)
	_Q1.FromAffine(&Q1)
	_res.Set(&_Q1).AddAssign(&_Q0)
	res.FromJacobian(&_res)
	return res, nil
}

// ----------------------------------------------------------------------------------------
// G2Affine

// https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-06#section-4.1
// Shallue and van de Woestijne method, works for any elliptic curve in Weierstrass curve
func svdwMapG2(u fptower.E2) G2Affine {

	var res G2Affine

	// constants
	// sage script to find z: https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-06#appendix-E.1
	var z, c1, c2, c3, c4 fptower.E2
	z.A0.SetUint64(3)
	z.A1.SetOne()
	c1.A0.SetUint64(8)
	c1.A1.SetString("16857842227199225999646786835636637546980511455593181936643515867446449681250963223102065881329922732942164707714798074010252040453159496268796306973476256192757039177873952554508334906059908")
	c2.A0.SetString("16857842227199225999646786835636637546980511455593181936643515867446449681250963223102065881329922732942164707714798074010252040453159496268796306973476256192757039177873952554508334906059908")
	c3.A0.SetString("1558368830738803275349191582824076350469659905985053863787578252089423397226386761007877077017425397430888291929495114607835172204836958213805200837320679070769936597535825877280079363727884")
	c3.A1.SetString("13970901216447696439790762183452979468072948204674220807153539451802418909813609701163289640364991206264590799408075205392342243124614515256445217727464822437358491433858694389471223006170314")
	c4.A0.SetString("1983275556141085411723151392427839711409471935952139051369825396170170550735407438012007750744696792110842906789976244001206122406254058384564271408644265434442004609161641477000980577183516")
	c4.A1.SetString("5949826668423256235169454177283519134228415807856417154109476188510511652206222314036023252234090376332528720369928732003618367218762175153692814225932796303326013827484924431002941731550556")

	var tv1, tv2, tv3, tv4, one, x1, gx1, x2, gx2, x3, x, gx, y fptower.E2
	one.SetOne()
	tv1.Square(&u).Mul(&tv1, &c1)
	tv2.Add(&one, &tv1)
	tv1.Sub(&one, &tv1)
	tv3.Mul(&tv2, &tv1).Inverse(&tv3)
	tv4.Mul(&u, &tv1)
	tv4.Mul(&tv4, &tv3)
	tv4.Mul(&tv4, &c3)
	x1.Sub(&c2, &tv4)
	gx1.Square(&x1)
	gx1.Add(&gx1, &aTwistCurveCoeff)
	gx1.Mul(&gx1, &x1)
	e1 := gx1.Legendre()
	x2.Add(&c2, &tv4)
	gx2.Square(&x2)
	gx2.Add(&gx2, &aTwistCurveCoeff)
	gx2.Mul(&gx2, &x2)
	e2 := gx2.Legendre() - e1 // 2 if is_square(gx2) AND NOT e1
	x3.Square(&tv2)
	x3.Mul(&x3, &tv3)
	x3.Square(&x3)
	x3.Mul(&x3, &c4)
	x3.Add(&x3, &z)
	if e1 == 1 {
		x.Set(&x1)
	} else {
		x.Set(&x3)
	}
	if e2 == 2 {
		x.Set(&x2)
	}
	gx.Square(&x)
	gx.Add(&gx, &aTwistCurveCoeff)
	gx.Mul(&gx, &x)
	y.Sqrt(&gx)
	e3 := sign0(u.A0) && sign0(y.A0)
	if !e3 {
		y.Neg(&y)
	}
	res.X.Set(&x)
	res.Y.Set(&y)

	return res
}

// MapToCurveG2Svdw maps an fp.Element to a point on the curve using the Shallue and van de Woestijne map
// https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-06#section-2.2.1
func MapToCurveG2Svdw(t fptower.E2) G2Affine {
	res := svdwMapG2(t)
	res.ClearCofactor(&res)
	return res
}

// EncodeToCurveG2Svdw maps an fp.Element to a point on the curve using the Shallue and van de Woestijne map
// https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-06#section-2.2.2
func EncodeToCurveG2Svdw(msg, dst []byte) (G2Affine, error) {
	var res G2Affine
	_t, err := hashToFp(msg, dst, 2)
	if err != nil {
		return res, err
	}
	var t fptower.E2
	t.A0.Set(&_t[0])
	t.A1.Set(&_t[1])
	res = MapToCurveG2Svdw(t)
	return res, nil
}

// HashToCurveG2Svdw maps an fp.Element to a point on the curve using the Shallue and van de Woestijne map
// https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-06#section-3
func HashToCurveG2Svdw(msg, dst []byte) (G2Affine, error) {
	var res G2Affine
	u, err := hashToFp(msg, dst, 4)
	if err != nil {
		return res, err
	}
	var u0, u1 fptower.E2
	u0.A0.Set(&u[0])
	u0.A1.Set(&u[1])
	u1.A0.Set(&u[2])
	u1.A1.Set(&u[3])
	Q0 := MapToCurveG2Svdw(u0)
	Q1 := MapToCurveG2Svdw(u1)
	var _Q0, _Q1, _res G2Jac
	_Q0.FromAffine(&Q0)
	_Q1.FromAffine(&Q1)
	_res.Set(&_Q1).AddAssign(&_Q0)
	res.FromJacobian(&_res)
	return res, nil
}
