package element

const Inverse = `

{{/* We use big.Int for Inverse for these type of moduli */}}
{{if eq .NoCarry false}}

// Inverse z = x^-1 mod q 
// note: allocates a big.Int (math/big)
func (z *{{.ElementName}}) Inverse( x *{{.ElementName}}) *{{.ElementName}} {
	inverse(z, x)
	return z
}

func _inverseGeneric(z, x *{{.ElementName}})  {
	var _xNonMont big.Int
	x.ToBigIntRegular( &_xNonMont)
	_xNonMont.ModInverse(&_xNonMont, Modulus())
	z.SetBigInt(&_xNonMont)
}

{{ else }}

// Inverse z = x^-1 mod q 
// Algorithm 16 in "Efficient Software-Implementation of Finite Fields with Applications to Cryptography"
// if x == 0, sets and returns z = x 
func (z *{{.ElementName}}) Inverse(x *{{.ElementName}}) *{{.ElementName}} {
	inverse(z, x)
	return z
}



// _inverseGeneric z = x^-1 mod q 
// Algorithm 16 in "Efficient Software-Implementation of Finite Fields with Applications to Cryptography"
// if x == 0, sets and returns z = x 
func  _inverseGeneric(z, x *{{.ElementName}})  {
	if x.IsZero() {
		z.SetZero()
		return
	}

	// initialize u = q
	var u = {{.ElementName}}{
		{{- range $i := .NbWordsIndexesFull}}
		{{index $.Q $i}},{{end}}
	}

	// initialize s = r^2
	var s = {{.ElementName}}{
		{{- range $i := .RSquare}}
		{{$i}},{{end}}
	}

	// r = 0
	r := {{.ElementName}}{}

	v := *x

	var carry, borrow uint64
	var bigger bool

	for  {
		for v[0]&1 == 0 {
			{{ rsh "v" .NbWords}}
			if s[0]&1 == 1 {
				{{ template "add_q" dict "all" . "V1" "s" }}
			}
			{{ rsh "s" .NbWords}}
		} 
		for u[0]&1 == 0 {
			{{ rsh "u" .NbWords}}
			if r[0]&1 == 1 {
				{{ template "add_q" dict "all" . "V1" "r" }}
			}
			{{ rsh "r" .NbWords}}
		} 
		{{ template "bigger" dict "all" . "V1" "v" "V2" "u"}}
		if bigger  {
			{{ template "sub_noborrow" dict "all" . "V1" "v" "V2" "u" "OmitLast" "true"}}
			{{ template "sub_noborrow" dict "all" . "V1" "s" "V2" "r" "OmitLast" "false"}}
			if borrow == 1 {
				{{ template "add_q" dict "all" . "V1" "s" }}
			}
		} else {
			{{ template "sub_noborrow" dict "all" . "V1" "u" "V2" "v" "OmitLast" "true"}}
			{{ template "sub_noborrow" dict "all" . "V1" "r" "V2" "s" "OmitLast" "false"}}
			if borrow == 1 {
				{{ template "add_q" dict "all" . "V1" "r" }}
			}
		}
		if (u[0] == 1) && ({{- range $i := reverse .NbWordsIndexesNoZero}}u[{{$i}}] {{if eq $i 1}}{{else}} | {{end}}{{end}} ) == 0 {
			z.Set(&r)
			return
		}
		if (v[0] == 1) && ({{- range $i := reverse .NbWordsIndexesNoZero}}v[{{$i}}] {{if eq $i 1}}{{else}} | {{end}}{{end}} ) == 0 {
			z.Set(&s)
			return
		}
	}

}

{{ end }}




{{ define "bigger" }}
	// {{$.V1}} >= {{$.V2}}
	bigger = !({{- range $i := reverse $.all.NbWordsIndexesNoZero}} {{$.V1}}[{{$i}}] < {{$.V2}}[{{$i}}] || ( {{$.V1}}[{{$i}}] == {{$.V2}}[{{$i}}] && (
		{{- end}}{{$.V1}}[0] < {{$.V2}}[0] {{- range $i :=  $.all.NbWordsIndexesNoZero}} )) {{- end}} )
{{ end }}

{{ define "add_q" }}
	// {{$.V1}} = {{$.V1}} + q 
	{{$.V1}}[0], carry = bits.Add64({{$.V1}}[0], {{index $.all.Q 0}}, 0)
	{{- range $i := .all.NbWordsIndexesNoZero}}
		{{- if eq $i $.all.NbWordsLastIndex}}
			{{$.V1}}[{{$i}}], _ = bits.Add64({{$.V1}}[{{$i}}], {{index $.all.Q $i}}, carry)
		{{- else}}
			{{$.V1}}[{{$i}}], carry = bits.Add64({{$.V1}}[{{$i}}], {{index $.all.Q $i}}, carry)
		{{- end}}
	{{- end}}
{{ end }}

{{ define "sub_noborrow" }}
	// {{$.V1}} = {{$.V1}} - {{$.V2}}
	{{$.V1}}[0], borrow = bits.Sub64({{$.V1}}[0], {{$.V2}}[0], 0)
	{{- range $i := .all.NbWordsIndexesNoZero}}
		{{- if and (eq $i $.all.NbWordsLastIndex) (eq "true" $.OmitLast)}}
		{{$.V1}}[{{$i}}], _ = bits.Sub64({{$.V1}}[{{$i}}], {{$.V2}}[{{$i}}], borrow)
		{{- else}}
		{{$.V1}}[{{$i}}], borrow = bits.Sub64({{$.V1}}[{{$i}}], {{$.V2}}[{{$i}}], borrow)
		{{- end}}
	{{- end}}
{{ end }}


{{ define "rsh V nbWords" }}
	// {{$.V}} = {{$.V}} >> 1
	{{$lastIndex := sub .nbWords 1}}
	{{- range $i :=  iterate .nbWords}}
		{{- if ne $i $lastIndex}}
			{{$.V}}[{{$i}}] = {{$.V}}[{{$i}}] >> 1 | {{$.V}}[{{(add $i 1)}}] << 63
		{{- end}}
	{{- end}}
	{{$.V}}[{{$lastIndex}}] >>= 1
{{ end }}

`
