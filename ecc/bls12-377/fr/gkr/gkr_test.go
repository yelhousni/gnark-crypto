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

package gkr

import (
	"encoding/json"
	"fmt"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fr"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fr/mimc"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fr/polynomial"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fr/sumcheck"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fr/test_vector_utils"
	fiatshamir "github.com/consensys/gnark-crypto/fiat-shamir"
	"github.com/consensys/gnark-crypto/utils"
	"github.com/stretchr/testify/assert"
	"os"
	"path/filepath"
	"reflect"
	"strconv"
	"testing"
	"time"
)

func TestNoGateTwoInstances(t *testing.T) {
	// Testing a single instance is not possible because the sumcheck implementation doesn't cover the trivial 0-variate case
	testNoGate(t, []fr.Element{four, three})
}

func TestNoGate(t *testing.T) {
	testManyInstances(t, 1, testNoGate)
}

func TestSingleMulGateTwoInstances(t *testing.T) {
	testSingleMulGate(t, []fr.Element{four, three}, []fr.Element{two, three})
}

func TestSingleMulGate(t *testing.T) {
	testManyInstances(t, 2, testSingleMulGate)
}

func TestSingleInputTwoIdentityGatesTwoInstances(t *testing.T) {

	testSingleInputTwoIdentityGates(t, []fr.Element{two, three})
}

func TestSingleInputTwoIdentityGates(t *testing.T) {

	testManyInstances(t, 2, testSingleInputTwoIdentityGates)
}

func TestSingleInputTwoIdentityGatesComposedTwoInstances(t *testing.T) {
	testSingleInputTwoIdentityGatesComposed(t, []fr.Element{two, one})
}

func TestSingleInputTwoIdentityGatesComposed(t *testing.T) {
	testManyInstances(t, 1, testSingleInputTwoIdentityGatesComposed)
}

func TestSingleMimcCipherGateTwoInstances(t *testing.T) {
	testSingleMimcCipherGate(t, []fr.Element{one, one}, []fr.Element{one, two})
}

func TestSingleMimcCipherGate(t *testing.T) {
	testManyInstances(t, 2, testSingleMimcCipherGate)
}

func TestATimesBSquaredTwoInstances(t *testing.T) {
	testATimesBSquared(t, 2, []fr.Element{one, one}, []fr.Element{one, two})
}

func TestShallowMimcTwoInstances(t *testing.T) {
	testMimc(t, 2, []fr.Element{one, one}, []fr.Element{one, two})
}
func TestMimcTwoInstances(t *testing.T) {
	testMimc(t, 93, []fr.Element{one, one}, []fr.Element{one, two})
}

func TestMimc(t *testing.T) {
	testManyInstances(t, 2, generateTestMimc(93))
}

func generateTestMimc(numRounds int) func(*testing.T, ...[]fr.Element) {
	return func(t *testing.T, inputAssignments ...[]fr.Element) {
		testMimc(t, numRounds, inputAssignments...)
	}
}

func TestSumcheckFromSingleInputTwoIdentityGatesGateTwoInstances(t *testing.T) {
	circuit := Circuit{Wire{
		Gate:            IdentityGate{},
		Inputs:          []*Wire{},
		nbUniqueOutputs: 2,
	}}

	wire := &circuit[0]

	assignment := WireAssignment{&circuit[0]: []fr.Element{two, three}}
	var o settings
	pool := polynomial.NewPool(256, 1<<11)
	workers := utils.NewWorkerPool()
	o.pool = &pool
	o.workers = workers

	claimsManagerGen := func() *claimsManager {
		manager := newClaimsManager(circuit, assignment, o)
		manager.add(wire, []fr.Element{three}, five)
		manager.add(wire, []fr.Element{four}, six)
		return &manager
	}

	transcriptGen := test_vector_utils.NewMessageCounterGenerator(4, 1)

	proof, err := sumcheck.Prove(claimsManagerGen().getClaim(wire), fiatshamir.WithHash(transcriptGen(), nil))
	assert.NoError(t, err)
	err = sumcheck.Verify(claimsManagerGen().getLazyClaim(wire), proof, fiatshamir.WithHash(transcriptGen(), nil))
	assert.NoError(t, err)
}

var one, two, three, four, five, six fr.Element

func init() {
	one.SetOne()
	two.Double(&one)
	three.Add(&two, &one)
	four.Double(&two)
	five.Add(&three, &two)
	six.Double(&three)
}

var testManyInstancesLogMaxInstances = -1

func getLogMaxInstances(t *testing.T) int {
	if testManyInstancesLogMaxInstances == -1 {

		s := os.Getenv("GKR_LOG_INSTANCES")
		if s == "" {
			testManyInstancesLogMaxInstances = 5
		} else {
			var err error
			testManyInstancesLogMaxInstances, err = strconv.Atoi(s)
			if err != nil {
				t.Error(err)
			}
		}

	}
	return testManyInstancesLogMaxInstances
}

func testManyInstances(t *testing.T, numInput int, test func(*testing.T, ...[]fr.Element)) {
	fullAssignments := make([][]fr.Element, numInput)
	maxSize := 1 << getLogMaxInstances(t)

	t.Log("Entered test orchestrator, assigning and randomizing inputs")

	for i := range fullAssignments {
		fullAssignments[i] = make([]fr.Element, maxSize)
		setRandom(fullAssignments[i])
	}

	inputAssignments := make([][]fr.Element, numInput)
	for numEvals := maxSize; numEvals <= maxSize; numEvals *= 2 {
		for i, fullAssignment := range fullAssignments {
			inputAssignments[i] = fullAssignment[:numEvals]
		}

		t.Log("Selected inputs for test")
		test(t, inputAssignments...)
	}
}

func testNoGate(t *testing.T, inputAssignments ...[]fr.Element) {
	c := Circuit{
		{
			Inputs: []*Wire{},
			Gate:   nil,
		},
	}

	assignment := WireAssignment{&c[0]: inputAssignments[0]}

	proof, err := Prove(c, assignment, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(1, 1)))
	assert.NoError(t, err)

	// Even though a hash is called here, the proof is empty

	err = Verify(c, assignment, proof, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(1, 1)))
	assert.NoError(t, err, "proof rejected")
}

func testSingleMulGate(t *testing.T, inputAssignments ...[]fr.Element) {

	c := make(Circuit, 3)
	c[2] = Wire{
		Gate:   mulGate{},
		Inputs: []*Wire{&c[0], &c[1]},
	}

	assignment := WireAssignment{&c[0]: inputAssignments[0], &c[1]: inputAssignments[1]}.Complete(c)

	proof, err := Prove(c, assignment, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(1, 1)))
	assert.NoError(t, err)

	err = Verify(c, assignment, proof, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(1, 1)))
	assert.NoError(t, err, "proof rejected")

	err = Verify(c, assignment, proof, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(0, 1)))
	assert.NotNil(t, err, "bad proof accepted")
}

func testSingleInputTwoIdentityGates(t *testing.T, inputAssignments ...[]fr.Element) {
	c := make(Circuit, 3)

	c[1] = Wire{
		Gate:   IdentityGate{},
		Inputs: []*Wire{&c[0]},
	}

	c[2] = Wire{
		Gate:   IdentityGate{},
		Inputs: []*Wire{&c[0]},
	}

	assignment := WireAssignment{&c[0]: inputAssignments[0]}.Complete(c)

	proof, err := Prove(c, assignment, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(0, 1)))
	assert.NoError(t, err)

	err = Verify(c, assignment, proof, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(0, 1)))
	assert.NoError(t, err, "proof rejected")

	err = Verify(c, assignment, proof, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(1, 1)))
	assert.NotNil(t, err, "bad proof accepted")
}

func testSingleMimcCipherGate(t *testing.T, inputAssignments ...[]fr.Element) {
	c := make(Circuit, 3)

	c[2] = Wire{
		Gate:   mimcCipherGate{},
		Inputs: []*Wire{&c[0], &c[1]},
	}

	t.Log("Evaluating all circuit wires")
	assignment := WireAssignment{&c[0]: inputAssignments[0], &c[1]: inputAssignments[1]}.Complete(c)
	t.Log("Circuit evaluation complete")
	proof, err := Prove(c, assignment, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(0, 1)))
	assert.NoError(t, err)
	t.Log("Proof complete")
	err = Verify(c, assignment, proof, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(0, 1)))
	assert.NoError(t, err, "proof rejected")

	t.Log("Successful verification complete")
	err = Verify(c, assignment, proof, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(1, 1)))
	assert.NotNil(t, err, "bad proof accepted")
	t.Log("Unsuccessful verification complete")
}

func testSingleInputTwoIdentityGatesComposed(t *testing.T, inputAssignments ...[]fr.Element) {
	c := make(Circuit, 3)

	c[1] = Wire{
		Gate:   IdentityGate{},
		Inputs: []*Wire{&c[0]},
	}
	c[2] = Wire{
		Gate:   IdentityGate{},
		Inputs: []*Wire{&c[1]},
	}

	assignment := WireAssignment{&c[0]: inputAssignments[0]}.Complete(c)

	proof, err := Prove(c, assignment, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(0, 1)))
	assert.NoError(t, err)

	err = Verify(c, assignment, proof, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(0, 1)))
	assert.NoError(t, err, "proof rejected")

	err = Verify(c, assignment, proof, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(1, 1)))
	assert.NotNil(t, err, "bad proof accepted")
}

func mimcCircuit(numRounds int) Circuit {
	c := make(Circuit, numRounds+2)

	for i := 2; i < len(c); i++ {
		c[i] = Wire{
			Gate:   mimcCipherGate{},
			Inputs: []*Wire{&c[i-1], &c[0]},
		}
	}
	return c
}

func testMimc(t *testing.T, numRounds int, inputAssignments ...[]fr.Element) {
	//TODO: Implement mimc correctly. Currently, the computation is mimc(a,b) = cipher( cipher( ... cipher(a, b), b) ..., b)
	// @AlexandreBelling: Please explain the extra layers in https://github.com/ConsenSys/gkr-mimc/blob/81eada039ab4ed403b7726b535adb63026e8011f/examples/mimc.go#L10

	c := mimcCircuit(numRounds)

	t.Log("Evaluating all circuit wires")
	assignment := WireAssignment{&c[0]: inputAssignments[0], &c[1]: inputAssignments[1]}.Complete(c)
	t.Log("Circuit evaluation complete")

	proof, err := Prove(c, assignment, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(0, 1)))
	assert.NoError(t, err)

	t.Log("Proof finished")
	err = Verify(c, assignment, proof, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(0, 1)))
	assert.NoError(t, err, "proof rejected")

	t.Log("Successful verification finished")
	err = Verify(c, assignment, proof, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(1, 1)))
	assert.NotNil(t, err, "bad proof accepted")
	t.Log("Unsuccessful verification finished")
}

func testATimesBSquared(t *testing.T, numRounds int, inputAssignments ...[]fr.Element) {
	// This imitates the MiMC circuit

	c := make(Circuit, numRounds+2)

	for i := 2; i < len(c); i++ {
		c[i] = Wire{
			Gate:   mulGate{},
			Inputs: []*Wire{&c[i-1], &c[0]},
		}
	}

	assignment := WireAssignment{&c[0]: inputAssignments[0], &c[1]: inputAssignments[1]}.Complete(c)

	proof, err := Prove(c, assignment, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(0, 1)))
	assert.NoError(t, err)

	err = Verify(c, assignment, proof, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(0, 1)))
	assert.NoError(t, err, "proof rejected")

	err = Verify(c, assignment, proof, fiatshamir.WithHash(test_vector_utils.NewMessageCounter(1, 1)))
	assert.NotNil(t, err, "bad proof accepted")
}

func setRandom(slice []fr.Element) {
	for i := range slice {
		slice[i].SetRandom()
	}
}

func generateTestProver(path string) func(t *testing.T) {
	return func(t *testing.T) {
		testCase, err := newTestCase(path)
		assert.NoError(t, err)
		proof, err := Prove(testCase.Circuit, testCase.FullAssignment, testCase.transcriptSetting())
		assert.NoError(t, err)
		assert.NoError(t, proofEquals(testCase.Proof, proof))
	}
}

func generateTestVerifier(path string) func(t *testing.T) {
	return func(t *testing.T) {
		testCase, err := newTestCase(path)
		assert.NoError(t, err)
		err = Verify(testCase.Circuit, testCase.InOutAssignment, testCase.Proof, testCase.transcriptSetting())
		assert.NoError(t, err, "proof rejected")
		testCase, err = newTestCase(path)
		assert.NoError(t, err)
		err = Verify(testCase.Circuit, testCase.InOutAssignment, testCase.Proof, fiatshamir.WithHash(&test_vector_utils.MapHash{Map: testCase.Hash}, []byte{1}))
		assert.NotNil(t, err, "bad proof accepted")
	}
}

func TestGkrVectors(t *testing.T) {

	testDirPath := "../../../../internal/generator/gkr/test_vectors"
	dirEntries, err := os.ReadDir(testDirPath)
	assert.NoError(t, err)
	for _, dirEntry := range dirEntries {
		if !dirEntry.IsDir() {

			if filepath.Ext(dirEntry.Name()) == ".json" {
				path := filepath.Join(testDirPath, dirEntry.Name())
				noExt := dirEntry.Name()[:len(dirEntry.Name())-len(".json")]

				t.Run(noExt+"_prover", generateTestProver(path))
				t.Run(noExt+"_verifier", generateTestVerifier(path))

			}
		}
	}
}

func proofEquals(expected Proof, seen Proof) error {
	if len(expected) != len(seen) {
		return fmt.Errorf("length mismatch %d ≠ %d", len(expected), len(seen))
	}
	for i, x := range expected {
		xSeen := seen[i]

		if xSeen.FinalEvalProof == nil {
			if seenFinalEval := x.FinalEvalProof.([]fr.Element); len(seenFinalEval) != 0 {
				return fmt.Errorf("length mismatch %d ≠ %d", 0, len(seenFinalEval))
			}
		} else {
			if err := test_vector_utils.SliceEquals(x.FinalEvalProof.([]fr.Element), xSeen.FinalEvalProof.([]fr.Element)); err != nil {
				return fmt.Errorf("final evaluation proof mismatch")
			}
		}
		if err := test_vector_utils.PolynomialSliceEquals(x.PartialSumPolys, xSeen.PartialSumPolys); err != nil {
			return err
		}
	}
	return nil
}

func benchmarkGkrMiMC(b *testing.B, nbInstances, mimcDepth int) {
	fmt.Println("creating circuit structure")
	c := mimcCircuit(mimcDepth)

	in0 := make([]fr.Element, nbInstances)
	in1 := make([]fr.Element, nbInstances)
	setRandom(in0)
	setRandom(in1)

	fmt.Println("evaluating circuit")
	start := time.Now().UnixMicro()
	assignment := WireAssignment{&c[0]: in0, &c[1]: in1}.Complete(c)
	solved := time.Now().UnixMicro() - start
	fmt.Println("solved in", solved, "μs")

	//b.ResetTimer()
	fmt.Println("constructing proof")
	start = time.Now().UnixMicro()
	_, err := Prove(c, assignment, fiatshamir.WithHash(mimc.NewMiMC()))
	proved := time.Now().UnixMicro() - start
	fmt.Println("proved in", proved, "μs")
	assert.NoError(b, err)
}

func BenchmarkGkrMimc19(b *testing.B) {
	benchmarkGkrMiMC(b, 1<<19, 91)
}

func BenchmarkGkrMimc17(b *testing.B) {
	benchmarkGkrMiMC(b, 1<<17, 91)
}

func TestTopSortTrivial(t *testing.T) {
	c := make(Circuit, 2)
	c[0].Inputs = []*Wire{&c[1]}
	sorted := topologicalSort(c)
	assert.Equal(t, []*Wire{&c[1], &c[0]}, sorted)
}

func TestTopSortDeep(t *testing.T) {
	c := make(Circuit, 4)
	c[0].Inputs = []*Wire{&c[2]}
	c[1].Inputs = []*Wire{&c[3]}
	c[2].Inputs = []*Wire{}
	c[3].Inputs = []*Wire{&c[0]}
	sorted := topologicalSort(c)
	assert.Equal(t, []*Wire{&c[2], &c[0], &c[3], &c[1]}, sorted)
}

func TestTopSortWide(t *testing.T) {
	c := make(Circuit, 10)
	c[0].Inputs = []*Wire{&c[3], &c[8]}
	c[1].Inputs = []*Wire{&c[6]}
	c[2].Inputs = []*Wire{&c[4]}
	c[3].Inputs = []*Wire{}
	c[4].Inputs = []*Wire{}
	c[5].Inputs = []*Wire{&c[9]}
	c[6].Inputs = []*Wire{&c[9]}
	c[7].Inputs = []*Wire{&c[9], &c[5], &c[2]}
	c[8].Inputs = []*Wire{&c[4], &c[3]}
	c[9].Inputs = []*Wire{}

	sorted := topologicalSort(c)
	sortedExpected := []*Wire{&c[3], &c[4], &c[2], &c[8], &c[0], &c[9], &c[5], &c[6], &c[1], &c[7]}

	assert.Equal(t, sortedExpected, sorted)
}

type WireInfo struct {
	Gate   string `json:"gate"`
	Inputs []int  `json:"inputs"`
}

type CircuitInfo []WireInfo

var circuitCache = make(map[string]Circuit)

func getCircuit(path string) (Circuit, error) {
	path, err := filepath.Abs(path)
	if err != nil {
		return nil, err
	}
	if circuit, ok := circuitCache[path]; ok {
		return circuit, nil
	}
	var bytes []byte
	if bytes, err = os.ReadFile(path); err == nil {
		var circuitInfo CircuitInfo
		if err = json.Unmarshal(bytes, &circuitInfo); err == nil {
			circuit := circuitInfo.toCircuit()
			circuitCache[path] = circuit
			return circuit, nil
		} else {
			return nil, err
		}
	} else {
		return nil, err
	}
}

func (c CircuitInfo) toCircuit() (circuit Circuit) {
	circuit = make(Circuit, len(c))
	for i := range c {
		circuit[i].Gate = gates[c[i].Gate]
		circuit[i].Inputs = make([]*Wire, len(c[i].Inputs))
		for k, inputCoord := range c[i].Inputs {
			input := &circuit[inputCoord]
			circuit[i].Inputs[k] = input
		}
	}
	return
}

var gates map[string]Gate

func init() {
	gates = make(map[string]Gate)
	gates["identity"] = IdentityGate{}
	gates["mul"] = mulGate{}
	gates["mimc"] = mimcCipherGate{} //TODO: Add ark
	gates["select-input-3"] = _select(2)
}

type mimcCipherGate struct {
	ark fr.Element
}

func (m mimcCipherGate) Evaluate(input ...fr.Element) (res fr.Element) {
	var sum fr.Element

	sum.
		Add(&input[0], &input[1]).
		Add(&sum, &m.ark)

	res.Square(&sum)    // sum^2
	res.Mul(&res, &sum) // sum^3
	res.Square(&res)    //sum^6
	res.Mul(&res, &sum) //sum^7

	return
}

func (m mimcCipherGate) Degree() int {
	return 7
}

type PrintableProof []PrintableSumcheckProof

type PrintableSumcheckProof struct {
	FinalEvalProof  interface{}     `json:"finalEvalProof"`
	PartialSumPolys [][]interface{} `json:"partialSumPolys"`
}

func unmarshalProof(printable PrintableProof) (Proof, error) {
	proof := make(Proof, len(printable))
	for i := range printable {
		finalEvalProof := []fr.Element(nil)

		if printable[i].FinalEvalProof != nil {
			finalEvalSlice := reflect.ValueOf(printable[i].FinalEvalProof)
			finalEvalProof = make([]fr.Element, finalEvalSlice.Len())
			for k := range finalEvalProof {
				if _, err := test_vector_utils.SetElement(&finalEvalProof[k], finalEvalSlice.Index(k).Interface()); err != nil {
					return nil, err
				}
			}
		}

		proof[i] = sumcheck.Proof{
			PartialSumPolys: make([]polynomial.Polynomial, len(printable[i].PartialSumPolys)),
			FinalEvalProof:  finalEvalProof,
		}
		for k := range printable[i].PartialSumPolys {
			var err error
			if proof[i].PartialSumPolys[k], err = test_vector_utils.SliceToElementSlice(printable[i].PartialSumPolys[k]); err != nil {
				return nil, err
			}
		}
	}
	return proof, nil
}

type TestCase struct {
	Circuit         Circuit
	Hash            *test_vector_utils.ElementMap
	Proof           Proof
	FullAssignment  WireAssignment
	InOutAssignment WireAssignment
}

type TestCaseInfo struct {
	Hash    string          `json:"hash"`
	Circuit string          `json:"circuit"`
	Input   [][]interface{} `json:"input"`
	Output  [][]interface{} `json:"output"`
	Proof   PrintableProof  `json:"proof"`
}

var testCases = make(map[string]*TestCase)

func newTestCase(path string) (*TestCase, error) {
	path, err := filepath.Abs(path)
	if err != nil {
		return nil, err
	}
	dir := filepath.Dir(path)

	tCase, ok := testCases[path]
	if !ok {
		var bytes []byte
		if bytes, err = os.ReadFile(path); err == nil {
			var info TestCaseInfo
			err = json.Unmarshal(bytes, &info)
			if err != nil {
				return nil, err
			}

			var circuit Circuit
			if circuit, err = getCircuit(filepath.Join(dir, info.Circuit)); err != nil {
				return nil, err
			}
			var _hash *test_vector_utils.ElementMap
			if _hash, err = test_vector_utils.ElementMapFromFile(filepath.Join(dir, info.Hash)); err != nil {
				return nil, err
			}
			var proof Proof
			if proof, err = unmarshalProof(info.Proof); err != nil {
				return nil, err
			}

			fullAssignment := make(WireAssignment)
			inOutAssignment := make(WireAssignment)

			sorted := topologicalSort(circuit)

			inI, outI := 0, 0
			for _, w := range sorted {
				var assignmentRaw []interface{}
				if w.IsInput() {
					if inI == len(info.Input) {
						return nil, fmt.Errorf("fewer input in vector than in circuit")
					}
					assignmentRaw = info.Input[inI]
					inI++
				} else if w.IsOutput() {
					if outI == len(info.Output) {
						return nil, fmt.Errorf("fewer output in vector than in circuit")
					}
					assignmentRaw = info.Output[outI]
					outI++
				}
				if assignmentRaw != nil {
					var wireAssignment []fr.Element
					if wireAssignment, err = test_vector_utils.SliceToElementSlice(assignmentRaw); err != nil {
						return nil, err
					}

					fullAssignment[w] = wireAssignment
					inOutAssignment[w] = wireAssignment
				}
			}

			fullAssignment.Complete(circuit)

			for _, w := range sorted {
				if w.IsOutput() {

					if err = test_vector_utils.SliceEquals(inOutAssignment[w], fullAssignment[w]); err != nil {
						return nil, fmt.Errorf("assignment mismatch: %v", err)
					}

				}
			}

			tCase = &TestCase{
				FullAssignment:  fullAssignment,
				InOutAssignment: inOutAssignment,
				Proof:           proof,
				Hash:            _hash,
				Circuit:         circuit,
			}

			testCases[path] = tCase
		} else {
			return nil, err
		}
	}

	return tCase, nil
}

func (c *TestCase) transcriptSetting(initialChallenge ...[]byte) fiatshamir.Settings {
	return fiatshamir.WithHash(&test_vector_utils.MapHash{Map: c.Hash}, initialChallenge...)
}

type mulGate struct{}

func (g mulGate) Evaluate(element ...fr.Element) (result fr.Element) {
	result.Mul(&element[0], &element[1])
	return
}

func (g mulGate) Degree() int {
	return 2
}

type _select int

func (g _select) Evaluate(in ...fr.Element) fr.Element {
	return in[g]
}

func (g _select) Degree() int {
	return 1
}
