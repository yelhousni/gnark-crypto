package gpoint

const Base = `
// Most algos for points operations are taken from http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html

import (
	{{if eq .CoordType "fp.Element"}}
		"github.com/consensys/gurvy/{{.Fpackage}}/fp"
	{{end}}
	"github.com/consensys/gurvy/{{.Fpackage}}/fr"
	"github.com/consensys/gurvy/utils/parallel"
	"github.com/consensys/gurvy/utils/debug"
	"sync"
)

type {{.PName}}CoordType = {{.CoordType}}

// {{.PName}}Jac is a point with {{.PName}}CoordType coordinates
type {{.PName}}Jac struct {
	X, Y, Z {{.PName}}CoordType
}

// {{.PName}}Proj point in projective coordinates
type {{.PName}}Proj struct {
	X, Y, Z {{.PName}}CoordType
}

// {{.PName}}Affine point in affine coordinates
type {{.PName}}Affine struct {
	X, Y {{.PName}}CoordType
}

// {{toLower .PName}}JacExtended parameterized jacobian coordinates (x=X/ZZ, y=Y/ZZZ, ZZ**3=ZZZ**2)
type {{toLower .PName}}JacExtended struct {
	X, Y, ZZ, ZZZ {{.PName}}CoordType
}

// SetInfinity sets p to O
func (p *{{toLower .PName}}JacExtended) SetInfinity() *{{toLower .PName}}JacExtended {
	p.X.SetOne()
	p.Y.SetOne()
	p.ZZ.SetZero()
	p.ZZZ.SetZero()
	return p
}

// ToAffine sets p in affine coords
func (p *{{toLower .PName}}JacExtended) ToAffine(Q *{{.PName}}Affine) *{{.PName}}Affine {
	Q.X.Inverse(&p.ZZ).MulAssign(&p.X)
	Q.Y.Inverse(&p.ZZZ).MulAssign(&p.Y)
	return Q
}

// ToJac sets p in affine coords
func (p *{{toLower .PName}}JacExtended) ToJac(Q *{{.PName}}Jac) *{{.PName}}Jac {
	Q.X.Mul(&p.ZZ, &p.X).MulAssign(&p.ZZ)
	Q.Y.Mul(&p.ZZZ, &p.Y).MulAssign(&p.ZZZ)
	Q.Z.Set(&p.ZZZ)
	return Q
}

// mAdd
// http://www.hyperelliptic.org/EFD/{{toLower .PName}}p/auto-shortw-xyzz.html#addition-madd-2008-s
func (p *{{toLower .PName}}JacExtended) mAdd(a *{{.PName}}Affine) *{{toLower .PName}}JacExtended {

	//if a is infinity return p
	if a.X.IsZero() && a.Y.IsZero() {
		return p
	}
	// p is infinity, return a
	if p.ZZ.IsZero() {
		p.X = a.X
		p.Y = a.Y
		p.ZZ.SetOne()
		p.ZZZ.SetOne()
		return p
	}

	var U2, S2, P, R, PP, PPP, Q, Q2, RR, X3, Y3 {{.PName}}CoordType

	// p2: a, p1: p
	U2.Mul(&a.X, &p.ZZ)
	S2.Mul(&a.Y, &p.ZZZ)
	if U2.Equal(&p.X) && S2.Equal(&p.Y) {
		return p.double(a)
	}
	P.Sub(&U2, &p.X)
	R.Sub(&S2, &p.Y)
	PP.Square(&P)
	PPP.Mul(&P, &PP)
	Q.Mul(&p.X, &PP)
	RR.Square(&R)
	X3.Sub(&RR, &PPP)
	Q2.AddAssign(&Q).AddAssign(&Q)
	p.X.Sub(&X3, &Q2)
	Y3.Sub(&Q, &p.X).MulAssign(&R)
	R.Mul(&p.Y, &PPP)
	p.Y.Sub(&Y3, &R)
	p.ZZ.MulAssign(&PP)
	p.ZZZ.MulAssign(&PPP)

	return p
}

// double point in ZZ coords
// http://www.hyperelliptic.org/EFD/{{toLower .PName}}p/auto-shortw-xyzz.html#doubling-dbl-2008-s-1
func (p *{{toLower .PName}}JacExtended) double(q *{{.PName}}Affine) *{{toLower .PName}}JacExtended {

	var U, S, M, _M, Y3 {{.PName}}CoordType

	U.Double(&q.Y)
	p.ZZ.Square(&U)
	p.ZZZ.Mul(&U, &p.ZZ)
	S.Mul(&q.X, &p.ZZ)
	_M.Square(&q.X)
	M.Double(&_M).
		AddAssign(&_M) // -> + a, but a=0 here
	p.X.Square(&M).
		SubAssign(&S).
		SubAssign(&S)
	Y3.Sub(&S, &p.X).MulAssign(&M)
	U.Mul(&p.ZZZ, &q.Y)
	p.Y.Sub(&Y3, &U)

	return p
}

// Set set p to the provided point
func (p *{{.PName}}Jac) Set(a *{{.PName}}Jac) *{{.PName}}Jac {
	p.X.Set(&a.X)
	p.Y.Set(&a.Y)
	p.Z.Set(&a.Z)
	return p
}

// Equal tests if two points (in Jacobian coordinates) are equal
func (p *{{.PName}}Jac) Equal(a *{{.PName}}Jac) bool {

	if p.Z.IsZero() && a.Z.IsZero() {
		return true
	}
	_p := {{.PName}}Affine{}
	_p.FromJacobian(p)

	_a := {{.PName}}Affine{}
	_a.FromJacobian(a)

	return _p.X.Equal(&_a.X) && _p.Y.Equal(&_a.Y)
}


// Equal tests if two points (in Affine coordinates) are equal
func (p *{{ .PName}}Affine) Equal(a *{{ .PName}}Affine) bool {
	return p.X.Equal(&a.X) && p.Y.Equal(&a.Y)
}

// Clone returns a copy of self
func (p *{{.PName}}Jac) Clone() *{{.PName}}Jac {
	return &{{.PName}}Jac{
		p.X, p.Y, p.Z,
	}
}


// Neg computes -G
func (p *{{.PName}}Jac) Neg(a *{{.PName}}Jac) *{{.PName}}Jac {
	p.Set(a)
	p.Y.Neg(&a.Y)
	return p
}

// Neg computes -G
func (p *{{.PName}}Affine) Neg(a *{{.PName}}Affine) *{{.PName}}Affine {
	p.X.Set(&a.X)
	p.Y.Neg(&a.Y)
	return p
}

// SubAssign substracts two points on the curve
func (p *{{.PName}}Jac) SubAssign(curve *Curve, a {{.PName}}Jac) *{{.PName}}Jac {
	a.Y.Neg(&a.Y)
	p.AddAssign(curve, &a)
	return p
}

// FromJacobian rescale a point in Jacobian coord in z=1 plane
// WARNING super slow function (due to the division)
func (p *{{.PName}}Affine) FromJacobian(p1 *{{.PName}}Jac) *{{.PName}}Affine {

	var bufs [3]{{.PName}}CoordType

	if p1.Z.IsZero() {
		p.X.SetZero()
		p.Y.SetZero()
		return p
	}

	bufs[0].Inverse(&p1.Z)
	bufs[2].Square(&bufs[0])
	bufs[1].Mul(&bufs[2], &bufs[0])

	p.Y.Mul(&p1.Y, &bufs[1])
	p.X.Mul(&p1.X, &bufs[2])

	return p
}

// FromJacobian converts a point from Jacobian to projective coordinates
func (p *{{ .PName}}Proj) FromJacobian(Q *{{ .PName}}Jac) *{{ .PName}}Proj {
	// memalloc
	var buf {{.PName}}CoordType
	buf.Square(&Q.Z)

	p.X.Mul(&Q.X, &Q.Z)
	p.Y.Set(&Q.Y)
	p.Z.Mul(&Q.Z, &buf)

	return p
}

func (p *{{.PName}}Jac) String(curve *Curve) string {
	if p.Z.IsZero() {
		return "O"
	}
	_p := {{.PName}}Affine{}
	_p.FromJacobian(p)
	_p.X.FromMont()
	_p.Y.FromMont()
	return "E([" + _p.X.String() + "," + _p.Y.String() + "]),"
}

// FromAffine sets p = Q, p in Jacboian, Q in affine
func (p *{{ .PName}}Jac) FromAffine(Q *{{ .PName}}Affine) *{{ .PName}}Jac {
	if Q.X.IsZero() && Q.Y.IsZero() {
		p.Z.SetZero()
		p.X.SetOne()
		p.Y.SetOne()
		return p
	}
	p.Z.SetOne()
	p.X.Set(&Q.X)
	p.Y.Set(&Q.Y)
	return p
}


func (p *{{.PName}}Affine) String(curve *Curve) string {
	var x, y {{.PName}}CoordType
	x.Set(&p.X)
	y.Set(&p.Y)
	return "E([" + x.String() + "," + y.String() + "]),"
}

// IsInfinity checks if the point is infinity (in affine, it's encoded as (0,0))
func (p *{{.PName}}Affine) IsInfinity() bool {
	return p.X.IsZero() && p.Y.IsZero()
}
`