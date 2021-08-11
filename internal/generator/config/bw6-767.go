package config

func init() {
	Curves = append(Curves, Curve{
		Name:         "bw6-767",
		CurvePackage: "bw6767",
		EnumID:       "BW6_767",
		FrModulus:    "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787",
		FpModulus:    "496597749679620867773432037469214230242402307330180853437434581099336634619713640485778675608223760166307530047354464605410050411581079376994803852937842168733702867087556948851016246640584660942486895230518034810309227309966899431",
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
