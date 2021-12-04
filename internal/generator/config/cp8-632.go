package config

func init() {
	Curves = append(Curves, Curve{
		Name:         "cp8-632",
		CurvePackage: "cp8632",
		EnumID:       "CP8_632",
		FrModulus:    "39705142709513438335025689890408969744933502416914749335064285505637884093126342347073617133569",
		FpModulus:    "16857842227199225999646786835636637546980511455593181936643515867446449681250963223102065881329922732942164707714798074010252040453159496268796306973476256192757039177873952554508334906059909",
		G1: Point{
			CoordType:        "fp.Element",
			PointName:        "g1",
			GLV:              true,
			CofactorCleaning: true,
			CRange:           []int{4, 5, 8, 16},
		},
		G2: Point{
			CoordType:        "fptower.E2",
			PointName:        "g2",
			GLV:              true,
			CofactorCleaning: true,
			CRange:           []int{4, 5, 8, 16},
		},
	})

}
