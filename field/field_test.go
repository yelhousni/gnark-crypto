package field

import (
	"crypto/rand"
	"fmt"
	"github.com/leanovate/gopter/gen"
	"math/big"
	mrand "math/rand"
	"testing"

	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"
)

//TODO: Use genF.Map to generate ints in field instead of using byte slices

func TestIntToMont(t *testing.T) {

	parameters := gopter.DefaultTestParameters()
	parameters.MinSuccessfulTests = 10
	properties := gopter.NewProperties(parameters)
	genF := genField(t)

	properties.Property("must recover initial non-montgomery value by repeated halving", prop.ForAll(
		func(f *Field, ib [][]uint8) (bool, error) {

			var i big.Int
			i.SetBytes(ib[0])
			i.Mod(&i, f.ModulusBig)

			// turn into mont
			var mont big.Int
			f.ToMont(&mont, &i)
			f.FromMont(&mont, &mont)

			return mont.Cmp(&i) == 0, nil
		}, genF, genUint8SliceSlice(1),
	))

	properties.Property("turning R into montgomery form must match the R value from field", prop.ForAll(
		func(f *Field) (bool, error) {
			// test if using the same R
			i := big.NewInt(1)
			i.Lsh(i, 64*uint(f.NbWords))
			f.ToMont(i, i)

			err := BigIntMatchUint64Slice(i, f.RSquare)
			return err == nil, err
		}, genF),
	)

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestBigIntMatchUint64Slice(t *testing.T) {
	parameters := gopter.DefaultTestParameters()
	parameters.MinSuccessfulTests = 10
	properties := gopter.NewProperties(parameters)
	genF := genField(t)

	properties.Property("random big.int must match uint64 slice made out of .Bytes()", prop.ForAll(
		func(f *Field, ib [][]uint8) (bool, error) {

			var i big.Int
			i.SetBytes(ib[0])
			bytes := i.Bytes()
			ints := make([]uint64, (len(bytes)-1)/8+1)

			for j := 0; j < len(bytes); j++ {
				ints[j/8] ^= uint64(bytes[len(bytes)-1-j]) << (8 * (j % 8))
			}

			err := BigIntMatchUint64Slice(&i, ints)
			return err == nil, err
		}, genF, genUint8SliceSlice(1)))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestQuadExtensionMul(t *testing.T) {

	verifyMul := func(base *Field, x8Slice [][]uint8, y8Slice [][]uint8) (bool, error) {
		var nonRes big.Int
		base.FromMont(&nonRes, &base.NonResidue)
		if !nonRes.IsInt64() {
			return false, fmt.Errorf("non-residue too large: %v", nonRes)
		}

		f := NewTower(base, 2, base.NonResidue.Int64())
		x := uint8SliceSliceToBigIntSlice(&f, x8Slice)
		y := uint8SliceSliceToBigIntSlice(&f, y8Slice)

		z := f.Mul(x, y)

		var z0, z1, u big.Int

		base.
			Mul(&z0, &x[0], &y[0]).
			Mul(&u, &x[1], &y[1]).
			Mul(&u, &u, big.NewInt(base.NonResidue.Int64())).
			Add(&z0, &z0, &u)

		base.
			Mul(&z1, &x[0], &y[1]).
			Mul(&u, &x[1], &y[0]).
			Add(&z1, &z1, &u)

		return z0.Cmp(&z[0]) == 0 && z1.Cmp(&z[1]) == 0, nil
	}
	genF := genField(t)
	parameters := gopter.DefaultTestParameters()

	parameters.MinSuccessfulTests = 10
	properties := gopter.NewProperties(parameters)
	properties.Property("multiplication should yield the correct value", prop.ForAll(verifyMul, genF, genUint8SliceSlice(2), genUint8SliceSlice(2)))
	properties.TestingRun(t, gopter.ConsoleReporter(false))

	parameters.MinSuccessfulTests = 4
	properties = gopter.NewProperties(parameters)
	properties.Property("multiplication should yield the correct value (small cases)", prop.ForAll(
		verifyMul,
		genF,
		genSmallUint8SliceSlice(2, 3),
		genSmallUint8SliceSlice(2, 3),
	))
	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestExponentiationBls12381G2(t *testing.T) {
	base, err := NewField("dummyName", "dummyElement", "0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab", false)

	if err != nil {
		t.Fatal(err)
	}

	f := NewTower(base, 2, -1)
	Z := make([]big.Int, 2)
	Z[0].SetInt64(-2)
	Z[1].SetInt64(-1)

	type expTestCase struct {
		pow *big.Int
		res []big.Int
	}

	cases := []expTestCase{
		{big.NewInt(2), f.FromInt64([]int64{3, 4})},
	}
	for _, c := range cases {
		res := f.Exp(Z, c.pow)
		if !f.Equal(c.res, res) {
			t.Error("Failed on power", c.pow.Int64())
		}
	}
}

const minNbWords = 5
const maxNbWords = 37

func genSmallUint8SliceSlice(outerSize int, max uint8) gopter.Gen {
	return gen.SliceOfN(
		outerSize,
		gen.SliceOfN(1, gen.UInt8Range(0, max)),
	)
}

func genUint8SliceSlice(outerSize int) gopter.Gen {
	return gen.SliceOfN(
		outerSize,
		gen.SliceOfN(maxNbWords*8, gen.UInt8()),
	)
}

func uint8SliceSliceToBigIntSlice(f *Extension, in [][]uint8) []big.Int {
	res := make([]big.Int, f.Degree)
	bytes := make([]byte, f.Base.NbWords*8)

	for i := 0; i < len(res); i++ {

		j := 0
		for ; j < len(bytes) && j < len(in[i]); j++ {
			bytes[j] = in[i][len(in[i])-j-1]
		}

		res[i].SetBytes(bytes[:j]).Mod(&res[i], f.Base.ModulusBig)
	}

	return res
}

func genField(t *testing.T) gopter.Gen {
	return func(genParams *gopter.GenParameters) *gopter.GenResult {

		genField := func() *Field {

			nbWords := minNbWords + mrand.Intn(maxNbWords-minNbWords)
			bitLen := nbWords*64 - 1 - mrand.Intn(64)

			modulus, err := rand.Prime(rand.Reader, bitLen)
			if err != nil {
				t.Error(err)
			}

			var field *Field
			field, err = NewField("dummy", "DummyElement", modulus.Text(10), false)

			if err == nil {
				if field.NbBits != bitLen || field.NbWords != nbWords {
					err = fmt.Errorf("mismatch: field.NbBits = %d, bitLen = %d, field.NbWords = %d, nbWords = %d", field.NbBits, bitLen, field.NbWords, nbWords)
				}
			}

			if err != nil {
				t.Error(err)
			}
			return field
		}

		field := genField()
		genResult := gopter.NewGenResult(field, gopter.NoShrinker)
		return genResult
	}
}
