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

// Code generated by consensys/gnark-crypto DO NOT EDIT

// Package fp contains field arithmetic operations for modulus = 0x126633...70000d.
//
// The API is similar to math/big (big.Int), but the operations are significantly faster (up to 20x for the modular multiplication on amd64, see also https://hackmd.io/@zkteam/modular_multiplication)
//
// The modulus is hardcoded in all the operations.
//
// Field elements are represented as an array, and assumed to be in Montgomery form in all methods:
// 	type Element [10]uint64
//
// Example API signature
// 	// Mul z = x * y mod q
// 	func (z *Element) Mul(x, y *Element) *Element
//
// and can be used like so:
// 	var a, b Element
// 	a.SetUint64(2)
// 	b.SetString("984896738")
// 	a.Mul(a, b)
// 	a.Sub(a, a)
// 	 .Add(a, b)
// 	 .Inv(a)
// 	b.Exp(b, new(big.Int).SetUint64(42))
//
// Modulus
// 	0x126633cc0f35f63fc1a174f01d72ab5a8fcd8c75d79d2c74e59769ad9bbda2f8152a6c0fadea490b8da9f5e83f57c497e0e8850edbda407d7b5ce7ab839c2253d369bd31147f73cd74916ea4570000d // base 16
// 	20494478644167774678813387386538961497669590920908778075528754551012016751717791778743535050360001387419576570244406805463255765034468441182772056330021723098661967429339971741066259394985997 // base 10
package fp
