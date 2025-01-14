// +build gofuzz

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

package kzg

import (
	"bytes"
	"github.com/consensys/gnark-crypto/ecc/bw6-764/fr"
	bw6764_pol "github.com/consensys/gnark-crypto/ecc/bw6-764/fr/polynomial"
	"github.com/consensys/gnark-crypto/polynomial"
)

const (
	fuzzInteresting = 1
	fuzzNormal      = 0
	fuzzDiscard     = -1
)

func Fuzz(data []byte) int {
	if len(data) == 0 {
		return fuzzDiscard
	}
	size := int(uint8(data[0])) + 2 // TODO fix min size in NewScheme
	if size > (1 << 15) {
		size = 1 << 15
	}
	r := bytes.NewReader(data[1:])
	var alpha, point fr.Element
	alpha.SetRawBytes(r)
	point.SetRawBytes(r)
	s := NewScheme(size, alpha)

	// create polynomials
	f := make([]polynomial.Polynomial, size/2)
	for i := 0; i < len(f); i++ {
		_f := make(bw6764_pol.Polynomial, size)
		for j := 0; j < len(_f); j++ {
			_f[j].SetRawBytes(r)
		}
		f[i] = &_f
	}

	// commit the polynomials
	digests := make([]polynomial.Digest, size/2)
	for i := 0; i < len(digests); i++ {
		digests[i], _ = s.Commit(f[i])

	}

	proof, err := s.BatchOpenSinglePoint(&point, digests, f)
	if err != nil {
		panic(err)
	}

	// verify the claimed values
	_proof := proof.(*BatchProofsSinglePoint)
	for i := 0; i < len(f); i++ {
		expectedClaim := f[i].Eval(point).(fr.Element)
		if !expectedClaim.Equal(&_proof.ClaimedValues[i]) {
			panic("inconsistant claimed values")
		}
	}

	// verify correct proof
	err = s.BatchVerifySinglePoint(digests, proof)
	if err != nil {
		panic(err)
	}

	return fuzzNormal
}
