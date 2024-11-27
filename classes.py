import yaml
from scipy.integrate import odeint
from rich.progress import track
from constants import *
from fmodel import *
from fphis import *
from fmath import *


class System:
	"""docstring for System"""
	def __init__(self, model_b, model_psi, num_points, args, database_file = 'pulsars.json'):
		self.model_b = model_b
		self.model_psi = model_psi
		self.num_points = num_points

		with open("pulsar.yml", "r") as f:
			database_yaml = yaml.safe_load(f)
			database = database_yaml['parameters']
			units = database_yaml['units']

		if model_psi == 'old':
			psi_args = {}
			lambda_i = [get_lambda_iv(0,1)]
			c_i = [get_c_iv(0,1)]
			other_args = {'lambda_i': lambda_i, 'c_i': c_i}
		elif model_psi == 'new':
			psi_args = {'num_sym': int(args[0]), 'num_asym': int(args[0])}
			lambda_i = [get_lambda_iv(0,psi_args['num_sym']), get_lambda_iv(1,psi_args['num_asym'])]
			c_i = [get_c_iv(0,psi_args['num_sym']), get_c_iv(1,psi_args['num_asym'])]
			r_base = np.hstack([np.linspace(0.1, 0.9, psi_args['num_sym']), np.linspace(0.1, 0.9, psi_args['num_asym'])])
			phi_base = np.hstack([np.linspace(90, 90, psi_args['num_sym']), np.linspace(-90, -90, psi_args['num_sym'])])
			other_args = {'lambda_i': lambda_i, 'c_i': c_i, 'r_base': r_base, 'phi_base': phi_base}

		if model_b == 'bgi':
			B_model = B_bgi
			b_args = {'f_s': 1.96, 'epsilon': float(args[2]) if len(args)>2 else 0.04}
		elif model_b == 'mhd':
			B_model = B_mhd
			b_args = {'f_s': 1.75}

		self.units = units
		self.num_pulsars = len(database)
		self.psi_args = psi_args
		self.b_args = b_args
		self.other_args = other_args
		self.B_model = B_model
		self.const = {	
						'model_b': self.model_b,
						'model_psi': self.model_psi,
						'num_pulsars': self.num_pulsars,
						'num_points': self.num_points,
						'psi_args': self.psi_args,
						'b_args': self.b_args,
						'other_args': self.other_args
					 }

		self.pulsars = []
		for k in database:
			self.pulsars.append(Pulsar(self, k['name'], k['chi'], k['beta'], k['P'], k['P_dot'], k['W50']))

	def __str__(self):
		text = "Model with:\n"
		text += "-- " + self.model_b + " magnetic field model,\n"
		text += "-- " + self.model_psi + " psi model,\n" 
		text += "-- " + str(self.num_points) + " number of points,\n" 
		text += "-- " + str(self.const['psi_args']) + ", " + str(self.const['b_args']) + " coefficients, \n"
		text += "please wait..."
		return text

	def get_gamma(self, num_pulsar, r, phi):
		return self.pulsars[num_pulsar].get_gamma(r, phi)

	def get_field(self, num_pulsar, r, phi):
		return self.pulsars[num_pulsar].get_field(r, phi)

	def get_psi(self, num_pulsar):
		return self.pulsars[num_pulsar].get_psi()
			 
	def get_hgap(self, num_pulsar):
		return self.pulsars[num_pulsar].get_hgap()

	def get_lambda(self, num_pulsar):
		return self.pulsars[num_pulsar].get_lambda()

	def get_profile(self, num_pulsar):
		return self.pulsars[num_pulsar].get_profile()


class Pulsar:
	"""docstring for Pulsar"""
	def __init__(self, system, name, chi, beta, P, P_dot, W50):
		self.system = system
		self.name = name
		self.chi_d = min(chi, 180-chi)
		self.chi = system.units['chi']*self.chi_d
		self.beta_d = abs(beta)
		self.beta = system.units['beta']*self.beta_d
		self.P = system.units['P']*P
		self.P_dot = system.units['P_dot']*P_dot
		self.W50 = W50_angle(self.P, system.units['W50']*W50)

		self.Omega = Omega(self.P)
		self.sigma = sigma(self.beta, self.W50)
		self.R_0 = R_0(self.P, system.b_args['f_s'])
		self.B = system.B_model(self.P, self.P_dot, self.chi, system.b_args)
		self.psi_m = psi_VC(self.Omega, self.B, self.R_0, self.chi, 0, 0, 0)
		self.const = {	
						'name': self.name,
						'chi': self.chi,
						'chi_d': self.chi_d,
						'beta': self.beta,
						'beta_d': self.beta_d,
						'P': self.P,
						'P_dot': self.P_dot,
						'Omega': self.Omega,
						'W50': self.W50,
						'sigma': self.sigma,
						'R_0': self.R_0,
						'B': self.B,
						'psi_m': self.psi_m,
						'model_args': {}
					 }

		if system.model_psi == 'old':
			self.model_args = {}
			self.psi = psi_old
			self.h_gap = h_gap_old
			self.field = field_old
		elif system.model_psi == 'new':
			self.model_args = get_psi_args(system.const, self.const)
			self.const['model_args'] = self.model_args
			self.psi = psi_new
			self.h_gap = h_gap_new
			self.field = field_new

	def __str__(self):
		text = "Pulsar with:\n"
		text += "-- " + self.name + " name,\n"
		text += "-- " + str(self.chi) + " chi,\n"
		text += "-- " + str(self.beta) + " beta,\n"
		text += "-- " + str(self.P) + " P,\n"
		text += "-- " + str(self.P_dot) + " P_dot,\n"
		text += "-- " + str(self.W50) + " W50,\n"
		text += "please wait..."
		return text

	def get_gamma(self, r, phi):
		x, y = decart(r, rad(phi))
		point = Point(self, x, y)
		return point.get_gamma()

	def get_field(self, r, phi):
		x, y = decart(r, rad(phi))
		point = Point(self, x, y)
		return point.get_field()

	def get_psi(self, k = 0, b = 0):
		R = np.linspace(-1,1, self.system.num_points)
		Psi, points = [], []
		for k in track(range(self.system.num_points), description="Calculating psi..."):
			points.append(Point(self, R[k], 0))
		for p in points:
			Psi.append(p.psi/self.psi_m)
		Psi = np.array(Psi)
		return (R, Psi)

	def get_hgap(self):
		R = np.linspace(-1,1, self.system.num_points)
		Hgap, points = [], []
		for k in track(range(self.system.num_points), description="Calculating hgap..."):
			points.append(Point(self, R[k], 0))
		for p in points:
			Hgap.append(p.h_gap)
		Hgap = np.array(Hgap)/self.R_0
		return (R, Hgap)

	def get_lambda(self):
		X = np.linspace(-1,1, self.system.num_points)
		Y = np.linspace(-1,1, self.system.num_points)
		Lambda, points = [], []
		for x in X:
			for y in Y:
				points.append(Point(self, x, y))
		for k in track(range(self.system.num_points**2), description="Calculating lambda..."):
			Lambda.append(lambda_m(points[k].system.const, points[k].pulsar.const, points[k].const, points[k].vars))
		Lambda = np.array(Lambda)
		return (X, Y, Lambda)

	def get_profile(self):
		R = np.linspace(-1,1, self.system.num_points)
		shift = sin(self.sigma)
		LambdaMP, points = [], []
		for k in range(self.system.num_points):
			points.append(Point(self, R[k], +shift))
		for k in track(range(self.system.num_points), description="Calculating MainPulse..."):
			LambdaMP.append(lambda_m(points[k].system.const, points[k].pulsar.const, points[k].const, points[k].vars))
		LambdaMP = np.array(LambdaMP)

		Lambda_m = max(LambdaMP)

		LambdaIP, points = [], []
		for k in range(self.system.num_points):
			points.append(Point(self, R[k], -shift))
		for k in track(range(self.system.num_points), description="Calculating InterPulse..."):
			LambdaIP.append(lambda_m(points[k].system.const, points[k].pulsar.const, points[k].const, points[k].vars))
		LambdaIP = np.array(LambdaIP)

		return (R, LambdaMP/Lambda_m, LambdaIP/Lambda_m)

class Point:
	"""docstring for Point"""
	def __init__(self, pulsar, x, y):
		self.system = pulsar.system
		self.pulsar = pulsar
		self.x = x
		self.y = y

		r, phi = polar(x,y)
		self.r, self.phi = r, deg(phi)
		self.r_m, self.phi_m, self.theta_m = spherical(self.r, self.phi, pulsar.R_0)
		self.R_c = R_c(self.r_m)
		self.theta_b = theta_b(pulsar.chi, self.r_m, self.phi_m)
		self.const = {	
						'x': self.x,
						'y': self.y,
						'r': self.r,
						'phi': self.phi,
						'r_m': self.r_m,
						'phi_m': self.phi_m,
						'theta_m': self.theta_m,
						'R_c': self.R_c,
						'theta_b': self.theta_b,
						'psi': 0,
						'h_gap': 0
					 }

		self.psi_h = pulsar.psi(self.system.const, self.pulsar.const, self.const)
		self.h_gap = pulsar.h_gap(self.system.const, self.pulsar.const, self.const)
		self.psi = self.psi_h(self.h_gap) - self.psi_h(0)
		self.const['h_gap'], self.const['psi'] = self.h_gap, self.psi 
		self.field_h = pulsar.field(self.system.const, self.pulsar.const, self.const)

		self.dif_args = (self.pulsar.R_0, self.R_c, self.field_h)
		self.vars = {
						'h_gap': self.h_gap,
						'psi': self.psi,
						'dif_args':self.dif_args
					}

	def __str__(self):
		text = "Point with:\n"
		text += "-- " + str(self.r) + " r,\n"
		text += "-- " + str(self.phi) + " phi,\n"
		text += "please wait..."
		return text

	def get_field(self):
		H = np.linspace(0, 1, self.system.num_points)
		Field, points = [], []
		for k in track(range(self.system.num_points), description="Calculating field..."):
			Field.append(self.field_h(H[k]*self.pulsar.R_0))
		Field = np.array(Field)
		return (H, Field)
	
	def get_gamma(self):
		h = np.linspace(0, 20, self.system.num_points)
		H = h*self.pulsar.R_0
		Gamma = odeint(eq_gamma, gamma0, H, args=(self.dif_args))
		Delta = odeint(eq_delta, gamma0, H, args=(self.dif_args))
		return (h, Gamma, Delta)
