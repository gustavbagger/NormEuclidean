package main

import (
	"errors"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"
)

type VariableState struct {
	q1, q2, h, f, u, sigma, phi uint64
	H, X, E                     float64
	Case, r                     int
}

func VariableStateNew(q1, q2, f uint64, lambda float64, r int) VariableState {
	H := float64(f) / float64(q1*q2)
	var Case int
	var u, sigma, phi uint64
	if q2 < (q1 << 1) {
		H /= 3
	} else {
		H /= 2
	}

	h := hSet(f, lambda, r)

	if q2 < h {
		Case = 1
		u = q1 * q2
		sigma = (q1 + 1) * (q2 + 1)
		phi = (q1 - 1) * (q2 - 1)
	} else if q1 < h {
		Case = 2
		u = q1
		sigma = (q1 + 1)
		phi = (q1 - 1)
	} else {
		Case = 3
		u = 1
		sigma = 1
		phi = 1
	}

	X := H / float64(h)
	E := 1 - (float64(sigma)/4.0+(float64(phi)/float64(u))+(float64(phi)/X))*
		(float64(sigma)/X)*
		math.Pi*
		math.Pi/6.0
	return VariableState{q1, q2, h, f, u, sigma, phi, H, X, E, Case, r}
}

func (vars VariableState) validateE() bool {
	return vars.E > 0
}

func (vars VariableState) validateX() bool {
	return vars.X > float64(vars.u)
}

func (vars VariableState) validateQ1Q2() bool {
	return vars.q1 < vars.q2
}

func (vars VariableState) dSet() (float64, error) {
	switch vars.r {
	case 1:
		fallthrough
	case 2:
		return 1, nil
	case 3:
		return 1 + float64(1)/(float64(6*vars.h)), nil
	case 4:
		return 1 + float64(2)/(float64(3*vars.h)), nil
	case 5:
		return 1 + float64(5)/(float64(3*vars.h)), nil
	case 6:
		return 1 + float64(10)/(float64(3*vars.h)) + float64(5)/(float64(36*vars.h*vars.h)), nil
	default:
		return 0, errors.New("Invalid choice of r")
	}
}

func hSet(f uint64, lambda float64, r int) uint64 {
	return uint64(math.Ceil(lambda * math.Pow(float64(f), 1.0/float64(2*r))))
}

func (vars VariableState) WSet() (float64, error) {
	d, err := vars.dSet()
	if err != nil {
		return 0, err
	}

	factorial := 1.0 / float64(vars.h)
	for i := 2; i <= vars.r; i++ {
		factorial *= float64(i) / float64(vars.h)
	}
	return float64(2*vars.r-1) + math.Sqrt(float64(vars.f))*d*factorial, nil
}

func (vars VariableState) validateH() bool {
	if !vars.validateE() || !vars.validateX() || !vars.validateQ1Q2() {
		return false
	} else {
		W, err := vars.WSet()
		if err != nil {
			log.Printf("error: ", err)
			return false
		}
		return vars.H*vars.H >
			(1.0/vars.E)*
				(math.Pi*math.Pi/6.0)*
				float64(vars.sigma)/float64(vars.phi)*
				float64(vars.u)*
				float64(vars.h)*
				math.Sqrt(float64(vars.f))*
				math.Pow(2.0*float64(vars.h)/float64(vars.h-uint64(3*(3-vars.Case))), float64(2*vars.r))*
				W
	}
}

func largestPossibleQ2(f uint64) uint64 {
	return uint64(math.Floor(1.821 * math.Pow(float64(f), 0.25) * math.Pow(math.Log(float64(f)), 1.5)))
}

func smallestFreeWin(f uint64, lambda float64, r, maxTest int) uint64 {
	var varsQ2Large, varsQ2Small VariableState
	q2Max := largestPossibleQ2(f)
	h := hSet(f, lambda, r) //This is to validateH when q2 ==h

	for q1 := uint64(3); q1 < q2Max; q1 += 2 {
		varsQ2Large = VariableStateNew(q1, q2Max, f, lambda, r)
		varsQ2Small = VariableStateNew(q1, h, f, lambda, r)
		if varsQ2Large.validateH() && (q1 >= h || varsQ2Small.validateH()) {
			continue
		} else {
			return q1 - 2
		}
	}
	return uint64(maxTest)
}

func smallestFreeWinVaryingLambda(f uint64, r, maxTest, steps int) uint64 {
	min := 0.5
	max := 1.5
	length := max - min
	courseness := length / float64(steps)
	q1Boundary := uint64(1)
	var q1Try uint64
	for lambda := min; lambda <= max; lambda += courseness {
		q1Try = smallestFreeWin(f, lambda, r, maxTest)
		if q1Try > q1Boundary {
			q1Boundary = q1Try
		}
	}
	return q1Boundary
}

// Need to fix that f>max(uint64)
func main() {
	if len(os.Args) != 2 {
		log.Fatalln("Wrong arguments. Usage: go run . <steps (int)>")
	}
	steps, err := strconv.Atoi(os.Args[1])
	if err != nil {
		log.Fatalln("Something went wrong with strconv")
	}
	oldq1 := uint64(7)
	for b := 14; b <= 20; b++ {
		for i := 0; i < 9*steps; i++ {
			if b == 20 && i >= 2*steps {
				os.Exit(0)
			}
			a := 1.0 + float64(i)/float64(steps)
			f := uint64(math.Ceil(a * math.Pow10(b)))
			currentq1 := smallestFreeWinVaryingLambda(f, 3, 1000, steps)
			if currentq1 > oldq1 {
				oldq1 = currentq1
				fmt.Printf("%v, %2.f*10^%v\n", currentq1, a, b)
			}
		}
	}
}
