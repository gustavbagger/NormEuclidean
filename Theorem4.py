import math as ma
import numpy as np
import sympy as sy
from sympy import primerange 
import matplotlib.pyplot as plt

#Case 1, q_1,q_2 < h
def VariablesC1(q_1:int,q_2:int,f:int,h:int) -> list:
    u = q_1 * q_2
    sigma = (q_1 + 1) * (q_2 + 1)
    phi = (q_1 - 1) * (q_2 - 1)
    if q_2 < 2 * q_1:
        H = f / (3 * u)
    else:
        H = f/(2 * u)
    l = 2
    return list([u,sigma,phi,H,l])

#Case 2, q_1 < h <= q_2
def VariablesC2(q_1:int,q_2:int,f:int,h:int) -> list:
    u = q_1
    sigma = (q_1 + 1)
    phi = (q_1 - 1)
    if q_2 < 2 * q_1:
        H = f / (3 * u * q_2)
    else:
        H = f/(2 * u * q_2)
    l = 1
    return list([u,sigma,phi,H,l])

#Case 3, h<= q_1,q_2
def VariablesC3(q_1:int,q_2:int,f:int,h:int) -> list:
    u = 1
    sigma = 1
    phi = 1
    if q_2 < 2 * q_1:
        H = f / (3 * q_1 *q_2)
    else:
        H = f / (2 * q_1 * q_2)
    l = 0
    return list([u,sigma,phi,H,l])

#Set variables according to what Case we are in
def Variablechoice(q_1:int,q_2:int,f:int,h:int,Case:int) -> list:
    if Case == 1:
        return VariablesC1(q_1,q_2,f,h)
    elif Case == 2:
        return VariablesC2(q_1,q_2,f,h)
    elif Case == 3:
        return VariablesC3(q_1,q_2,f,h)
    else:
        raise ValueError("Case does not exist")

#Set value of E_u(X) depending on case
def Efunc(q_1:int,q_2:int,f:int,h:int,Case:int) -> float:
    u,sigma,phi,H,l = Variablechoice(q_1,q_2,f,h,Case)
    X = H/h
    return 1 - (ma.pi**2/6) * (sigma/X) * ((sigma/4) + (phi/u) + (phi/X))

#Check if E_u(X)>0 and X>= u
def E_X_validQ(q_1:int,q_2:int,f:int,h:int,Case:int) -> bool:
    u,sigma,phi,H,l = Variablechoice(q_1,q_2,f,h,Case)
    X = H/h
    return (Efunc(q_1,q_2,f,h,Case)>0) and (X >= u)

#Set value of d_r(h), Trevino, Rem 2.1
def dfunc(h:int,r:int) -> float:
    if (r==1) or (r==2):
        return 1
    elif r==3:
        return 1 + (1/(6*h))
    elif r==4:
        return 1 + (2/(3*h))
    elif r==5:
        return 1 + (5/(3*h))
    elif r==6:
        return 1 + (10/(3*h))+(5/(36*h**2))
    else:
        raise ValueError("r is not in the correct range")

#Set value of W(f,h,r)
def Wfunc(f:int,h:int,r:int) -> float:
    return 2 * r - 1 + f**(1/2) * ma.factorial(r) * dfunc(h,r) / (h**r)

#Set value of h = ceiling(lambda * f^(1/2r))
def hfunc(f:int,lamb:float,r:int) -> int:
    return ma.ceil(lamb * f**(1/(2*r)))

#Decide case, dependant on size of q1,q2 vs h
def whatcase(q_1:int,q_2:int,f:int,lamb:float,r:int) -> int:
    h = hfunc(f,lamb,r)
    if q_2 < h:
        return 1
    elif q_1 < h:
        return 2
    else:
        return 3

#Check the main inequality in Th4
def Hcondition(q_1:int,q_2:int,f:int,lamb:float,r:int):
    Case = whatcase(q_1, q_2, f, lamb, r)
    h = hfunc(f,lamb,r)
    if not E_X_validQ(q_1, q_2, f, h, Case):
        return False
    E = Efunc(q_1,q_2,f,h,Case)
    W = Wfunc(f,h,r)
    u,sigma,phi,H,l = Variablechoice(q_1, q_2, f, h, Case)
    return (1 / E) * (ma.pi**2 / 6) * (sigma / phi) * u * h * f**(1/2) * ((2*h)/(h - 3*l))**(2*r) * W < H**2 and q_1<q_2

#Unconditional upper bound on q2
def q_2max(f:int) -> int:
    return sy.prevprime(1.821 * ( f**(1 / 4) ) * ma.log(f)**(3 / 2))


#test when q_2 is close enough to q_1 to give worse H value
def findq1boundarylambdafree(f:int,maxtest:int = 10**3,r:int = 3) -> int:
    q1 = 3
    q2max = q_2max(f)
    while q1<maxtest:
        lamb = 0.1
        test = 0
        while lamb<2:
            h = sy.prevprime(hfunc(f,lamb,r))
            if Hcondition(q1,q2max,f,lamb,r) and Hcondition(q1,h,f,lamb,r):
                test = 1
                break 
            elif Hcondition(q1,q2max,f,lamb,r) and q1>= h:
                test = 1
                break 
            lamb += 0.05
        if test == 0:
            return q1
        q1 = sy.nextprime(q1)

#Find prime p st if q1>=p, Th4 is satisfied. Fixed choice of lambda
def findq1boundarylambdafixed(f:int,maxtest:int = 10**3,r:int = 3,lamb:float = 3/4) -> int:
    q1 = 3
    q2max = q_2max(f)
    h = sy.prevprime(hfunc(f,lamb,r))
    while q1<maxtest:
        test = 0
        if Hcondition(q1,q2max,f,lamb,r) and Hcondition(q1,h,f,lamb,r):
            test = 1
        elif Hcondition(q1,q2max,f,lamb,r) and q1>= h:
            test = 1
        if test == 0:
            return sy.prevprime(q1)
        q1 = sy.nextprime(q1)

#Find the best possible lambda value in order to satisfy Th4 immediately for largest possible prime p
def findq1bound(f:int):
    lamb = 0.5
    q1bound = (1,lamb)
    while lamb< 1.5:
        lamb+=0.01
        q1try = findq1boundarylambdafixed(f,10**3,3,lamb)
        if q1bound[0] < q1try:
            q1bound = (q1try,lamb)
    q1bound = (sy.nextprime(q1bound[0]),q1bound[1])
    return q1bound

#Traverse range 10^14 - 10^20 via a*10^b 
#For each 14<= b(int) <= 20, increments a by 1/finiteness
#Prints "p (a,b)" st for all f>a*10^b, q1>=p implies f not Norm-Euclidean 
def q1bounds(fineness:int):
    c = fineness
    oldq1 = 7
    for b in range(14,20):
        for i in range(0,9*c):
            a = 1 + i/c
            f = ma.ceil(a*10**b)
            currentq1 = findq1bound(f)[0]
            if currentq1 > oldq1:
                oldq1 = currentq1
                print(currentq1,(round(a,2),b))
    for i in range(0,2*c):
        a = 1 + i/c
        f = ma.ceil(a*10**20)
        currentq1 = findq1bound(f)[0]
        if currentq1 > oldq1:
            oldq1 = currentq1
            print(currentq1,(round(a,2),b))
        
        
#------------depricated------------
def findq1boundarylambdafixedq2large(f:int,maxtest:int = 10**3,r:int = 3,lamb:float = 3/4) -> int:
    q1 = 3
    q2max = q_2max(f)

    while q1<maxtest:
        test = 0
        if Hcondition(q1,q2max,f,lamb,r):
            test = 1
        if test == 0:
            return sy.prevprime(q1)
        q1 = sy.nextprime(q1)
       
        
def plot_q1bounds():
    ypoints=np.array([])
    xpoints=np.array([])
    for n in range(69):
        a = ((n) % 10)+1
        b = (n)//10
        f = a * 10**(14+b)
        xpoints = np.append(xpoints,[n])
        ypoints = np.append(ypoints,[findq1bound(f)[0]])
    return plt.plot(xpoints,ypoints),plt.ylabel('log10(q)>y => win'), plt.xlabel('n')