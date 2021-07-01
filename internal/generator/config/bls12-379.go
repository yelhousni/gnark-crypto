package config

func init() {
	Curves = append(Curves, Curve{
		Name:         "bls12-379",
		CurvePackage: "bls12379",
		EnumID:       "BLS12_379",
		FrModulus:    "15567573732069904898445906795858855143537221806202697669674256244108007833601",
		FpModulus:    "647455824720115791999401948377863964948475498868568187799869012112522291430466516542182303552703940519176779595777",
		G1: Point{
			CoordType:        "fp.Element",
			PointName:        "g1",
			GLV:              true,
			CofactorCleaning: true,
			CRange:           defaultCRange(),
		},
		G2: Point{
			CoordType:        "fptower.E2",
			PointName:        "g2",
			GLV:              true,
			CofactorCleaning: true,
			CRange:           defaultCRange(),
		},
	})

}
