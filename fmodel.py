from fmath import *
from fphis import *
from scipy.optimize import root_scalar

def lambda_m(system_const, pulsar_const, point_const, point_vars):
	R_0, B = pulsar_const['R_0'], pulsar_const['B']
	dif_args = point_vars['dif_args']
	R_c, r = point_const['R_c'], point_const['r']
	eps_e = epsilon_e(dif_args, R_0)
	eps_ph = epsilon_ph(B, R_c)
	l_m = (eps_e/eps_ph)
	return l_m*hs(l_m-1)*hs(1 - r)


################################################################
#OLD MODEL PSI#
################################################################

def psi_old(system_const, pulsar_const, point_const):
	P, Omega, B = pulsar_const['P'], pulsar_const['Omega'], pulsar_const['B'], 
	R_0, chi = pulsar_const['R_0'], pulsar_const['chi']
	R_c, theta_b = point_const['R_c'], point_const['theta_b']
	r_m, phi_m, theta_m = point_const['r_m'], point_const['phi_m'], point_const['theta_m']

	pho_gj = pho_GJ(Omega, B, theta_b)
	psi_rs = lambda h: psi_RS(h, pho_gj)
	psi_vc = psi_VC(Omega, B, R_0, chi, r_m, phi_m, theta_m)
	psi_h = lambda h: hs(psi_vc - psi_rs(h))*psi_rs(h)

	return psi_h

def field_old(system_const, pulsar_const, point_const):
	lambda_00 = system_const['other_args']['lambda_i'][0][0]
	P, Omega, B = pulsar_const['P'], pulsar_const['Omega'], pulsar_const['B'], 
	R_0, chi = pulsar_const['R_0'], pulsar_const['chi']
	R_c, theta_b = point_const['R_c'], point_const['theta_b']
	r_m, phi_m, theta_m = point_const['r_m'], point_const['phi_m'], point_const['theta_m']

	h_rs = H_RS(P, B, R_c, theta_b)
	pho_gj = pho_GJ(Omega, B, theta_b)
	psi_rs = lambda h: psi_RS(h, pho_gj)
	psi_vc = psi_VC(Omega, B, R_0, chi, r_m, phi_m, theta_m)

	field_rs = lambda h: field_RS(h, h_rs, pho_gj, R_0)
	field_vc = lambda h: 0
	if (psi_rs(h_rs) < psi_vc):
		field_h = field_rs
	else:
		field_h = field_vc

	return field_h


def h_gap_old(system_const, pulsar_const, point_const):
	P, B = pulsar_const['P'], pulsar_const['B']
	R_c, theta_b = point_const['R_c'], point_const['theta_b']

	return H_RS(P, B, R_c, theta_b)

################################################################
#NEW MODEL PSI#
################################################################

def get_psi_args(system_const, pulsar_const):
	n, m = system_const['psi_args']['num_sym'], system_const['psi_args']['num_asym']
	r_base, phi_base = np.array(system_const['other_args']['r_base']), np.array(system_const['other_args']['phi_base'])
	lambda_i, c_i = np.array(system_const['other_args']['lambda_i']), np.array(system_const['other_args']['c_i'])
	P, B, R_0, chi = pulsar_const['P'], pulsar_const['B'], pulsar_const['R_0'], pulsar_const['chi']
	H, R, Phi = [], [], []
	for k in range(n+m):
		r_m, phi_m, theta_m = spherical(r_base[k], phi_base[k], R_0)
		R.append(r_m)
		Phi.append(phi_m)
		r_c = R_c(r_m)
		theta = theta_b(chi, r_m, phi_m)
		h_rs = H_RS(P, B, r_c, theta)
		H.append(h_rs)
	H, R, Phi = np.array(H), np.array(R), np.array(Phi)
	h, r, phi = H/R_0, R/R_0, Phi
	A = np.full((n+m,n+m), 0.0)
	b = np.full((n+m,1), 0.0)

	D = lambda_i
	C = np.array([cos(chi) * (phi*0+1), 3/4 * R_0/R_st * sin(chi) * sin(phi)])

	for j in range(n+m):
		for i in range(n):
			A[i,j] = C[0,j] * D[0,i] * Jv(0, lambda_i[0,i]*r[j]) * (exp(-lambda_i[0,i]*h[j]) + exp(+lambda_i[0,i]*h[j])) 
		for i in range(m):
			A[i+n,j] = C[1,j] * D[1,i] * Jv(1, lambda_i[1,i]*r[j]) * (exp(-lambda_i[1,i]*h[j]) + exp(+lambda_i[1,i]*h[j]))

		b0 = np.sum([c_i[0, i] * C[0,j] * D[0,i] * Jv(0, lambda_i[0,i]*r[j]) * exp(+lambda_i[0,i]*h[j]) for i in range(n)]) 
		b1 = np.sum([c_i[1, i] * C[1,j] * D[1,i] * Jv(1, lambda_i[1,i]*r[j]) * exp(+lambda_i[1,i]*h[j]) for i in range(m)])
		b[j] = b0 + b1

	a_i0, b_i0, a_i1, b_i1 = [], [], [], []
	X = linsolve(A, b)

	for k in range(n):
		a_i0.append(float(X[k]))
		b_i0.append(float(c_i[0,k]) - a_i0[k])
	for k in range(m):
		a_i1.append(float(X[k+n]))
		b_i1.append(float(c_i[1,k]) - a_i1[k])
	a_i = [a_i0, a_i1]
	b_i = [b_i0, b_i1]

	return {'a_i': a_i, 'b_i': b_i}

def psi_iv(h, v, i, system_const, pulsar_const, point_const):
	lambda_i = system_const['other_args']['lambda_i']
	a_i, b_i = pulsar_const['model_args']['a_i'], pulsar_const['model_args']['b_i']
	R_0, r_m = pulsar_const['R_0'], point_const['r_m']
	return (a_i[v][i]*exp(-lambda_i[v][i]*h/R_0) + b_i[v][i]*exp(+lambda_i[v][i]*h/R_0)) * Jv(v,lambda_i[v][i]*r_m/R_0)

def psi_new(system_const, pulsar_const, point_const):
	n, m = system_const['psi_args']['num_sym'], system_const['psi_args']['num_asym']
	Omega, B, R_0 = pulsar_const['Omega'], pulsar_const['B'], pulsar_const['R_0']
	chi, r_m, phi_m = pulsar_const['chi'], point_const['r_m'], point_const['phi_m']

	const0 = 1/2 * (Omega*B*R_0**2)/c_lt * cos(chi)
	const1 = 3/8 * (Omega*B*R_0**2)/c_lt * R_0/R_st * sin(chi)*sin(phi_m)
	f0 = lambda h: np.sum([psi_iv(h, 0, k, system_const, pulsar_const, point_const) for k in range(0,n)])
	f1 = lambda h: np.sum([psi_iv(h, 1, k, system_const, pulsar_const, point_const) for k in range(0,m)])

	psi = lambda h: const0 * (1 - (r_m/R_0)**2 - f0(h)) + const1 * ((r_m/R_0) - (r_m/R_0)**3 - f1(h))
	psi_h = lambda h: psi(h)*hs(psi(h))
	return psi_h

def field_iv(h, v, i, system_const, pulsar_const, point_const):
	lambda_i = system_const['other_args']['lambda_i']
	a_i, b_i = pulsar_const['model_args']['a_i'], pulsar_const['model_args']['a_i']
	R_0, r_m = pulsar_const['R_0'], point_const['r_m']
	delta = delta_h*R_0
	dif_coef = (exp(lambda_i[v][i]*delta/R_0) - exp(-lambda_i[v][i]*delta/R_0))/(2*delta)
	return dif_coef * (a_i[v][i]*exp(-lambda_i[v][i]*h/R_0) - b_i[v][i]*exp(+lambda_i[v][i]*h/R_0)) * Jv(v,lambda_i[v][i]*r_m/R_0)

def field_test(system_const, pulsar_const, point_const):
	n, m = system_const['psi_args']['num_sym'], system_const['psi_args']['num_asym']
	Omega, B, R_0 = pulsar_const['Omega'], pulsar_const['B'], pulsar_const['R_0']
	chi, r_m, phi_m = pulsar_const['chi'], point_const['r_m'], point_const['phi_m']

	const0 = 1/2 * (Omega*B*R_0**2)/c_lt * cos(chi)
	const1 = 3/8 * (Omega*B*R_0**2)/c_lt * R_0/R_st * sin(chi)*sin(phi_m)
	f0 = lambda h: np.sum([field_iv(h, 0, k, system_const, pulsar_const, point_const) for k in range(0,n)])
	f1 = lambda h: np.sum([field_iv(h, 1, k, system_const, pulsar_const, point_const) for k in range(0,m)])

	field = lambda h: const0 * (1 - (r_m/R_0)**2 - f0(h)) + const1 * ((r_m/R_0) - (r_m/R_0)**3 - f1(h))
	field_test = lambda h: field(h)*hs(field(h))
	return field_test

def field_new(system_const, pulsar_const, point_const):
	h_gap = point_const['h_gap']
	field = field_test(system_const, pulsar_const, point_const)
	field_h = lambda h: field(h)*hs(h_gap - h)
	return field_h

def h_gap_new(system_const, pulsar_const, point_const):
	R_0 = pulsar_const['R_0']
	log_field = lambda h: log(field_test(system_const, pulsar_const, point_const)(h) + 1)
	h_gap = get_null(log_field, 0, R_st, delta_f*R_0)
	return h_gap

