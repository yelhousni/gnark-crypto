package ecc

import (
	"path/filepath"
	"strings"

	"github.com/consensys/bavard"
	"github.com/consensys/gnark-crypto/internal/generator/config"
)

func Generate(conf config.Curve, baseDir string, bgen *bavard.BatchGenerator) error {
	packageName := strings.ReplaceAll(conf.Name, "-", "")

	g1 := pconf{conf, conf.G1}
	g2 := pconf{conf, conf.G2}

	entries := []bavard.Entry{
		{File: filepath.Join(baseDir, "doc.go"), Templates: []string{"doc.go.tmpl"}},
		{File: filepath.Join(baseDir, "multiexp.go"), Templates: []string{"multiexp.go.tmpl"}},
		{File: filepath.Join(baseDir, "multiexp_test.go"), Templates: []string{"tests/multiexp.go.tmpl"}},
		{File: filepath.Join(baseDir, "marshal.go"), Templates: []string{"marshal.go.tmpl"}},
		{File: filepath.Join(baseDir, "marshal_test.go"), Templates: []string{"tests/marshal.go.tmpl"}},
	}
	conf.Package = packageName
	if err := bgen.Generate(conf, packageName, "./ecc/template", entries...); err != nil {
		return err
	}

	// fuzz testing
	entries = []bavard.Entry{
		{File: filepath.Join(baseDir, "fuzz.go"), Templates: []string{"fuzz.go.tmpl"}, BuildTag: "gofuzz"},
		{File: filepath.Join(baseDir, "fuzz_test.go"), Templates: []string{"tests/fuzz.go.tmpl"}, BuildTag: "gofuzz"},
	}
	if err := bgen.Generate(conf, packageName, "./ecc/template", entries...); err != nil {
		return err
	}

	// hash To curve
	if conf.HashE1 != nil {
		entries = []bavard.Entry{
			{File: filepath.Join(baseDir, "sswu-fp.go"), Templates: []string{"sswu-fp.go.tmpl"}},
			{File: filepath.Join(baseDir, "sswu-fp_test.go"), Templates: []string{"tests/sswu-fp.go.tmpl"}},
		}

		hashConf := config.NewHashSuiteInfo(conf.Fp, &conf.G1, conf.Name, conf.HashE1)

		if err := bgen.Generate(hashConf, packageName, "./ecc/template", entries...); err != nil {
			return err
		}
	}

	// G1
	entries = []bavard.Entry{
		{File: filepath.Join(baseDir, "g1.go"), Templates: []string{"point.go.tmpl"}},
		{File: filepath.Join(baseDir, "g1_test.go"), Templates: []string{"tests/point.go.tmpl"}},
	}
	if err := bgen.Generate(g1, packageName, "./ecc/template", entries...); err != nil {
		return err
	}

	// G2
	entries = []bavard.Entry{
		{File: filepath.Join(baseDir, "g2.go"), Templates: []string{"point.go.tmpl"}},
		{File: filepath.Join(baseDir, "g2_test.go"), Templates: []string{"tests/point.go.tmpl"}},
	}
	return bgen.Generate(g2, packageName, "./ecc/template", entries...)
}

type pconf struct {
	config.Curve
	config.Point
}
