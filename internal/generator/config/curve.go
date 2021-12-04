package config

import (
	"github.com/consensys/gnark-crypto/ecc"
	"github.com/consensys/gnark-crypto/field"
)

// Curve describes parameters of the curve useful for the template
type Curve struct {
	Name         string
	CurvePackage string
	Package      string // current package being generated
	EnumID       string
	FpModulus    string
	FrModulus    string

	Fp           *field.Field
	Fr           *field.Field
	FpUnusedBits int
	G1           Point
	G2           Point
}

func (c *Curve) ID() ecc.ID {
	switch c.Name {
	case "bn254":
		return ecc.BN254
	case "bls12-381":
		return ecc.BLS12_381
	case "bls12-377":
		return ecc.BLS12_377
	case "bw6-761":
		return ecc.BW6_761
	case "bw6-633":
		return ecc.BW6_633
	case "bls24-315":
		return ecc.BLS24_315
	case "bls12-379":
		return ecc.BLS12_379
	// case "bw6-764":
	// 	return ecc.BW6_764
	// case "bw6-672":
	// 	return ecc.BW6_672
	case "cp8-632":
		return ecc.CP8_632
	default:
		panic("not implemented")
	}
}

type Point struct {
	CoordType        string
	PointName        string
	GLV              bool  // scalar mulitplication using GLV
	CofactorCleaning bool  // flag telling if the Cofactor cleaning is available
	CRange           []int // multiexp bucket method: generate inner methods (with const arrays) for each c
}

var Curves []Curve

func defaultCRange() []int {
	// default range for C values in the multiExp
	return []int{4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 20, 21, 22}
}
