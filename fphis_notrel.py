import numpy as np
from scipy.integrate import odeint

from constants import *
from fmath import *


def Omega(P):
	return 2*PI/P

def W50_angle(P, W50_time):
	return 2*PI*(W50_time/P)

def sigma(beta, W50_angle):
	return np.arctan(beta/W50_angle)

def B_bgi(P, P_dot, b_args):
	epsilon = b_args['epsilon']
	const = exp10(12) * epsilon**(-1/2)
	p = P/exp10(0)
	p_dot = P_dot/exp10(-15)
	return const * p**(3/4) * p_dot**(1/2)

def B_mhd(P, P_dot, b_args):
	const = exp10(12) * 1.04/2**0.5  
	p = P/exp10(0)
	p_dot = P_dot/exp10(-15)
	return const * p**(1/2) * p_dot**(1/2)

def R_0(P, f_s):
	omega = Omega(P)
	return f_s**(1/2) * (omega*R_st/c_lt)**(1/2) * R_st

def R_c(r_m):
	try:
		return 4/3 * R_st**2 / r_m
	except ZeroDivisionError:
		return inf

def theta_b(chi, r_m, phi_m):
	return chi - 3/2 * (r_m/R_st) * sin(phi_m)

def H_RS(P, B, R_c, theta_b):
	const = 1.1 * exp10(4)
	cos_t = abs(cos(theta_b))
	p = P/exp10(0)
	r_c = R_c/exp10(7)
	b = B/exp10(12) 
	return const * cos_t**(-3/7) * p**(3/7) * r_c**(2/7) * b**(-4/7)

def pho_GJ(Omega, B, theta_b):
	return Omega*B*cos(theta_b)/(2*PI*c_lt)

def psi_RS(h, pho_GJ):
	psi_rs = 2*PI*pho_GJ*h**2
	return psi_rs * hs(psi_rs)

def psi_VC(Omega, B, R_0, chi, r_m, phi_m, theta_m):
	theta_0 = R_0/R_st
	const = 1/2 * (Omega*B*R_0**2)/c_lt * (1 - (theta_m/theta_0)**2)
	return const * (cos(chi) + 3/4*theta_m*sin(phi_m)*sin(chi))

def field_RS(h, H_RS, pho_GJ, R_0):
	delta_psi_rs = 2*psi_RS(H_RS, pho_GJ)/H_RS*(1 - h/H_RS)
	return delta_psi_rs * hs(delta_psi_rs)

def field_VC(h, lambda_00, Omega, B, R_0, chi, r_m, phi_m, theta_m):
	delta_psi_VC = psi_VC(Omega, B, R_0, chi, r_m, phi_m, theta_m) * lambda_00/R_0 * exp(-lambda_00*h/R_0)
	return delta_psi_VC

def epsilon_e(dif_args,R_0):
	H = np.linspace(0, R_st, exp10(h_int))
	gamma = odeint(eq_gamma, gamma0, H, args=(dif_args))
	delta = odeint(eq_delta, gamma0, H, args=(dif_args))
	gamma1, delta1 = gamma[-1][0], delta[-1][0]
	return (delta1 - gamma1) * mc2

def L_r(L0):
	return L0 - 3*ln(L0)

def L_0(eph, B, R_c):
	return ln(e**2/(h_pl) * (omega_b0*R_c/c_lt) * (B_cr/B)**2  * (mc2/eph)**2)

def e_ph(L0, B, R_c):
	return 8/(3*L_r(L0)) * (mc2) * (B_cr/B) * (R_c/R_st)

def epsilon_ph(B, R_c):
	L0 = L0_mid
	for k in range(h_int):
		eph = e_ph(L0, B, R_c)
		L0 = L_0(eph, B, R_c)
	return eph

def l_cr():
	pass