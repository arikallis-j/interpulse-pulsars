import numpy as np
from math import ceil
from scipy.special import jv
from scipy.optimize import root_scalar
from scipy.integrate import quad

from constants import *

def rad(x):
	return x*PI/180

def deg(x):
	return x*180/PI

def sin(x):
	return np.sin(x)
	
def cos(x):
	return np.cos(x)

def ln(x):
	return np.log(x)

def log(x):
	return np.log(x)/np.log(10)

def exp(x):
	return np.exp(x)

def exp10(x):
	return 10**x

def Jv(v,x):
	return jv(v,x)

def spherical(r, phi, R_0):
	r_m = r * R_0
	phi_m = rad(phi)
	theta_m = r_m / R_st
	return r_m, phi_m, theta_m

def polar(x,y):
	r = np.sqrt(x*x + y*y)
	phi = arc(x,y)
	return r, phi

def decart(r, phi):
	x = r*cos(phi)
	y = r*sin(phi)
	return (x,y)

def arc(x,y):
	try:
		phi = np.arctan(abs(y/x))
	except ZeroDivisionError:
		if y>=0:
			return PI/2
		else:
			return PI + PI/2

	if ((x>0) and (y>=0)):
		return phi
	if ((x<0) and (y>=0)):
		return PI - phi
	if ((x<0) and (y<=0)):
		return PI + phi
	if ((x>0) and (y<=0)):
		return 2*PI - phi

def hs(x):
	return np.heaviside(x,1)

def zero_div(x,y):
	if y!= 0:
		return x/y
	else:
		return inf

def fon(r):
	if r>1:
		return 1
	else:
		return None

def get_null(f, a, b, h):
	x = a
	while x<=b:
		if f(x+h)*f(x)<=0:
			return x+h/2
		else:
			x+=h
	return b

def get_lambda_iv(v,n):
	root = 0
	k = 0
	lamda_i = []
	J = lambda x: Jv(v, x)
	while (len(lamda_i) < n):
		try:
			t_s = root_scalar(J, bracket = [k, k+1]).root
			if root!=t_s:
				root = t_s
				lamda_i.append(root)
		except:
			pass
		k += 1
	return lamda_i

def get_c_iv(v, n):
	lambda_i = get_lambda_iv(v, n)
	c_i = []
	Jp = lambda x: Jv(v+1, x)
	for k in range(n):
		Ii = lambda x: x**(v+1) * (1-x**2) * Jv(v,lambda_i[k]*x)
		vi = quad(Ii, 0, 1)
		ci = 2*vi[0]/((Jv(v+1,lambda_i[k]))**2)
		c_i.append(ci)
	return c_i

def linsolve(M, v):
	X = np.linalg.solve(M, v)
	return X

def eq_gamma(g, h, r_0, r_c, dpdh):
	dgdh = r_cl/e * dpdh(h) - 2/3 * r_cl/r_c**2 * g**4
	return dgdh

def eq_delta(g, h, r_0, r_c, dpdh):
	dgdh = r_cl/e * dpdh(h)
	return dgdh