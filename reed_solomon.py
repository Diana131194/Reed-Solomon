
import sys
from sympy import abc
import sympy as sp
import random
import math
import numpy as np


##############################ENCODER##############################


def get_msg(k):  # get m=(m0,m1,...,mk) from user
    msg = []
    print("Enter", k, " values:")
    for i in range(0, k):
        element = int(input())
        msg.append(element)
    return msg

def getEvaluationPoints(msg,n,q): # returns the list of evaluation points -> [(1,fm(1)),...,(n,fm(n))]
    sum = 0
    points = []
    for i in range(0, n+1):  
        for j in range(0, k):  # sum = m0 + m1 * x + ... + m(k-1) * x ^ (k-1)
            sum = sum + msg[j]*(i**j)  # sum = sum + mj * x^j, sum gets another element
        sum = sum % q
        points.append((i,sum))
        sum = 0
    return points

def get_next_prime(n):  # returns next prime greater than n
    x = n+1
    for num in range(x,2*x):
        isPrime = True
        for j in range(2,int(sqrt(num))): # enough to check factors until sqrt(num)
            if num % j == 0:
                isPrime = False
                break
        if isPrime:
            return num

def getNoisyMsg(encoded_msg,numOfErrors): #insert numOfErrors errors to encoded msg, first numOfErrors evaluation points are mutated
    noisy = []
    for i in range(0,int(numOfErrors)): #copy with noise
        tup = encoded_msg[i]
        second = tup[1]
        noisy.append((i,second-1))
    for i in range(int(numOfErrors),len(encoded_msg)): #copy without noise
        noisy.append(encoded_msg[i])
    return noisy


k = int(input("Enter k: \n"))

msg = get_msg(k) #get msg from user

n = random.randint(4*(k**2),5*(k**2)) # in order to be able to fix at least one error, n-2k*sqrt(n) > 0 <===> n > 4*k^2

q = get_next_prime(n)  # q is the field , first prime greater than n, q>n>k

points = getEvaluationPoints(msg,n,q) #points = [(0,fm(0)),...,(n,fm(n))], the encoded message

n_sqrt = int(math.ceil(math.sqrt(n)))-1

maxNumOfErrors = math.floor(n-math.floor(2*k*n_sqrt)) #upper bound for # of errors

noisy_msg = getNoisyMsg(points,maxNumOfErrors)  #insert maxNumOfErrors errors to message

print "k:",k,"q:",q,"n:",n
print "Number of errors:",int(maxNumOfErrors)
print "Message coefficients:",msg




##############################ENCODER-END##############################


##############################DECODER##############################


def buildEquation(point,n): #build a single row i in the matrix
     n_sqrt = int(math.ceil(math.sqrt(n)))
     eq = []
     for i in range(n_sqrt): #deg(Qx),deg(Qy)<=sqrt(n)
            for j in range(n_sqrt):
                eq.append((point[0]**i)*(point[1]**j)) # point[0] = ai, point[1] = yi
     return eq

def buildEquationsSystem(points,n): #build the matrix s.t to all i's, Q(ai,yi)=0
    eq=[]
    for point in points:
        eq.append(buildEquation(point,n))
    return np.array(eq)


def buildPoly(A,n): #build a polynomial corresponds to the null space matrix A, s.t Aij := coefficient of x^i * y^j
    R = GF(q)['x,y']
    x,y = R.gens()
    polynomial = 0
    n_sqrt = int(math.ceil(math.sqrt(n)))

    output = []
    for i in range(n_sqrt):
        for j in range(n_sqrt):
            output.append(x^i * y^j)
    for k in range(len(A)):
        polynomial = polynomial + A[k] * output[k]
    return polynomial

def getFinalList(factorizedPoly): # given Q(x,y) = Q1(x,y)Q2(x,y)...Qt(x,y) find all
    finalList = []                # Qi's of the form y-pi(x) s.t deg(pi(x))=k-1
    for factor in factorizedPoly:
        R = GF(q)['x,y']
        x,y = R.gens()
        poly = 0
        poly = poly + factor[0]
        if poly.degree(y)==1 and poly.degree(x)==k-1:
            poly = (poly*(-1)) + y
            finalList.append(poly)
    return finalList



equations = buildEquationsSystem(noisy_msg,n)

matrix = Matrix(GF(q),equations)

nullSpace = matrix.right_kernel()[1]         #find the the matrix B which zeroes the equation system(matrix X B = 0),
                                             #folded to a vector v=B[:]. Finding the null space means to find the polynomial
                                             # Q(x,y) that gives Q(ai,yi) = 0 for all i's.

polynomial = buildPoly(nullSpace,n) # find the polynomial corresponds to the null space vector

factorizedPoly = list(polynomial.factor()) #decompose Q(x,y) to factors, i.e Q(x,y) = Q1(x,y)Q2(x,y)...Qt(x,y)

finalList = getFinalList(factorizedPoly) #get all polynomials p(x) with deg(p(x)) = k-1

print "Q(x,y) factorized:",factorizedPoly
print "Number of factors: ",len(factorizedPoly)
print "Final list:",finalList



##############################DECODER-END##############################











