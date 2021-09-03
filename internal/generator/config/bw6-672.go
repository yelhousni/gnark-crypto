package config

func init() {
	Curves = append(Curves, Curve{
		Name:         "bw6-672",
		CurvePackage: "bw6672",
		EnumID:       "BW6_672",
		FrModulus:    "39705142709513438335025689890408969744933502416914749335064285505637884093126342347073617133569",
		FpModulus:    "10298672283669831003294162666372200112722413308864335408449041451781560030635870917520584453632890853794509234527426523862901218215709508411167619988692340017052754441213777163725504913151409615827435527",
		G1: Point{
			CoordType:        "fp.Element",
			PointName:        "g1",
			GLV:              true,
			CofactorCleaning: true,
			CRange:           []int{4, 5, 8, 16},
		},
		G2: Point{
			CoordType:        "fp.Element",
			PointName:        "g2",
			GLV:              true,
			CofactorCleaning: true,
			CRange:           []int{4, 5, 8, 16},
		},
	})

}
