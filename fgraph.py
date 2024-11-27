import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
from rich.console import Console
from rich.table import Table

from constants import *
from fmath import *

def graph_gamma(H, Gamma, Delta, form, name):
	if form=='table':
		table = Table(title = 'Pulsar ' + name['name'] + ' ' + name['model_b'].upper() + ' (' + name['r'] + ', ' + name['phi'] + '°' + ')', title_style="bold cyan", header_style="cyan")
		table.add_column("N", justify = "right", style="cyan")
		table.add_column("h/R_0", justify = "center", style="red")
		table.add_column("gamma", justify = "center", style="blue")
		table.add_column("delta", justify = "center", style="green")

		for k in range(len(H)):
			table.add_row(str(k+1), str(H[k]), str(Gamma[k][0]), str(Delta[k][0]))
		console = Console()
		console.print('')
		console.print(table)
		return 0
	plt.style.use('seaborn-whitegrid')
	fig, ax = plt.subplots()

	ax.set_xlim(0,20)
	ax.plot(H, Delta, 'b-', linewidth=2, label='$\gamma_{0}$')
	ax.plot(H, Gamma, 'g-', linewidth=2, label='$\gamma_{e}$')

	ax.set_title("Pulsar " + name['name'] + ' ' + name['model_b'].upper() + ' (' + name['r'] + ', ' + name['phi'] + '°' + ')', loc='center', fontsize=20)
	ax.set_xlabel('$h/R_{0}$')
	ax.set_ylabel('$\gamma$')
	ax.legend()

	if form=='save':
		if not os.path.isdir('graph'):
			os.mkdir('graph')
		plt.savefig('graph/' + 'gamma_' + name['name'] + '_' + name['model_b'] + '_' + name['r'] + '_' + name['phi'] + '.png')
	elif form=='print':
		plt.show()
	else: 
		print("FormError")

def graph_field(H, Field, form, name):
	if form=='table':
		table = Table(title = 'Pulsar ' + name['name'] + ' ' + name['model_b'].upper() + ' (' + name['r'] + ', ' + name['phi'] + '°' + ')', title_style="bold cyan", header_style="cyan")
		table.add_column("N", justify = "right", style="cyan")
		table.add_column("h/R_0", justify = "center", style="red")
		table.add_column("E_field", justify = "center", style="blue")

		for k in range(len(H)):
			table.add_row(str(k+1), str(H[k]), str(Field[k]))
		console = Console()
		console.print('')
		console.print(table)
		return 0
	plt.style.use('seaborn-whitegrid')
	fig, ax = plt.subplots()
	ax.set_xlim(0,20)
	#ax.set_ylim(0,2*10**13)
	ax.plot(H, Field, 'b-', linewidth=2, label='$field$')

	ax.set_title("Pulsar " + name['name'] + ' ' + name['model_b'].upper() + ' (' + name['r'] + ', ' + name['phi'] + '°' + ')', loc='center', fontsize=20)
	ax.set_xlabel('$r/R_{0}$')
	ax.set_ylabel('$Field E$')
	ax.legend()

	if form=='save':
		if not os.path.isdir('field'):
			os.mkdir('field')
		plt.savefig('field/'+ 'field_' + name['name'] + '_' + name['model_b'] + '_' + name['r'] + '_' + name['phi'] +'.png')
	elif form=='print':
		plt.show()
	else: 
		print("FormError")

def graph_psi(R, Psi, form, name):
	if form=='table':
		table = Table(title = 'Pulsar ' + name['name'] + ' ' + name['model_b'].upper(), title_style="bold cyan", header_style="cyan")
		table.add_column("N", justify = "right", style="cyan")
		table.add_column("r_m/R_0", justify = "center", style="red")
		table.add_column("psi/psi_max", justify = "center", style="blue")

		for k in range(len(R)):
			table.add_row(str(k+1), str(R[k]), str(Psi[k]))
		console = Console()
		console.print('')
		console.print(table)
		return 0

	plt.style.use('seaborn-whitegrid')
	fig, ax = plt.subplots()

	ax.set_xlim(-1,1)
		
	ax.plot(R, Psi, 'b-', linewidth=2, label='$\psi_{NW}$')

	ax.set_title("Pulsar " + name['name'] + ' ' + name['model_b'].upper(), fontsize=20)
	ax.set_xlabel('$r/R_{0}$')
	ax.set_ylabel('$\psi/\psi_{max}$')
	ax.legend()

	if form=='save':
		if not os.path.isdir('psi'):
			os.mkdir('psi')
		plt.savefig('psi/'+ 'psi_' + name['name'] + '_' + name['model_b'] + '.png')
	elif form=='print':
		plt.show()
	else: 
		print("FormError")

def graph_hgap(R, Hgap, form, name):
	if form=='table':
		table = Table(title = 'Pulsar ' + name['name'] + ' ' + name['model_b'].upper(), title_style="bold cyan", header_style="cyan")
		table.add_column("N", justify = "right", style="cyan")
		table.add_column("r_m/R_0", justify = "center", style="red")
		table.add_column("h_gap/R_0", justify = "center", style="blue")

		for k in range(len(R)):
			table.add_row(str(k+1), str(R[k]), str(Hgap[k]))
		console = Console()
		console.print(table)
		return 0

	plt.style.use('seaborn-whitegrid')
	fig, ax = plt.subplots()

	ax.set_xlim(-1,1)
	#ax.set_ylim(0, 12*10**4)
		
	ax.plot(R, Hgap, 'r-', linewidth=2, label='$h_{gap}$')

	ax.set_title("Pulsar " + name['name'] + ' ' + name['model_b'].upper(), fontsize=20)
	ax.set_xlabel('$r/R_{0}$')
	ax.set_ylabel('$h_{gap}/R_{0}$')
	ax.legend()

	if form=='save':
		if not os.path.isdir('h_gap'):
			os.mkdir('h_gap')
		plt.savefig('h_gap/'+ 'h_gap_' + name['name'] + '_' + name['model_b'] + '.png')
	elif form=='print':
		plt.show()
	else: 
		print("FormError")

def graph_lambda(X, Y, Lambda, form, name):
	if form=='table':
		table = Table(title = 'Pulsar ' + name['name'] + ' ' + name['model_b'].upper(), title_style="bold cyan", header_style="cyan")
		table.add_column("N", justify = "right", style="cyan")
		table.add_column("x/R_0", justify = "center", style="red")
		table.add_column("y/R_0", justify = "center", style="blue")
		table.add_column("lambda", justify = "center", style="green")

		for i in range(len(X)):
			for j in range(len(Y)):
				k = i*len(X)+j
				table.add_row(str(k+1), str(X[i]), str(Y[j]), str(Lambda[k]))
		console = Console()
		console.print('')
		console.print(table)
		return 0

	Xf = np.linspace(-1, 1, N_fon)
	Yf = np.linspace(-1, 1, N_fon)

	x_arr, y_arr = np.meshgrid(X, Y)
	xf_arr, yf_arr = np.meshgrid(Xf, Yf)
	z_arr, zf_arr, N = x_arr*0, xf_arr*0, len(X)

	for i in range(N):
		for j in range(N):
			k = i*N + j
			z_arr[j][i] = Lambda[k]

	for i in range(N_fon):
		for j in range(N_fon):
			k = i*N_fon + j
			r, phi = polar(xf_arr[j][i], yf_arr[j][i])
			zf_arr[j][i] = fon(r)

	plt.style.use('_classic_test_patch')
	fig, ax = plt.subplots()

	norm = mpl.colors.Normalize(round(np.min(z_arr)),round(np.max(z_arr)))

	fig.colorbar(mpl.cm.ScalarMappable(norm = norm))
	ax.contourf(x_arr, y_arr, z_arr, N_lev, norm = norm)
	ax.contourf(xf_arr, yf_arr, zf_arr, N_lev, norm = norm, colors = ['w'])

	"""
	if line:
		x = np.linspace(-cos(pulsar.sigma),cos(pulsar.sigma), N)
		y = +sin(pulsar.sigma)+(x*0)
		ax.plot(x, y, 'b-', linewidth=2, label='$\gamma_{0}$')
	"""

	fig.set_figwidth(9)  #ширина и
	fig.set_figheight(8) #высота "Figure"

	ax.set_title("Pulsar " + name['name'] + ' ' + name['model_b'].upper(), fontsize=20)
	ax.set_xlabel('$X$')
	ax.set_ylabel('$Y$')

	if form=='save':
		if not os.path.isdir('pic'):
			os.mkdir('pic')
		plt.savefig('pic/'+ 'lambda_' + name['name'] + '_' + name['model_b'] + '.png')
	elif form=='print':
		plt.show()
	else:
		print("FormError")

def graph_profile(R, LambdaMP, LambdaIP, form, name):
	if form=='table':
		table = Table(title = 'Pulsar ' + name['name'] + ' ' + name['model_b'].upper(), title_style="bold cyan", header_style="cyan")
		table.add_column("N", justify = "right", style="cyan")
		table.add_column("r_m/R_0", justify = "center", style="red")
		table.add_column("J_mp/J_max", justify = "center", style="blue")
		table.add_column("J_ip/J_max", justify = "center", style="green")

		for k in range(len(R)):
			table.add_row(str(k+1), str(R[k]), str(LambdaMP[k]), str(LambdaIP[k]))
		console = Console()
		console.print('')
		console.print(table)
		return 0
	plt.style.use('seaborn-whitegrid')
	fig, ax = plt.subplots()

	ax.set_xlim(-1,1)
	ax.set_title("Pulsar " + name['name'] + ' ' + name['model_b'].upper(), fontsize=20)
	ax.plot(R, LambdaMP, 'b-', linewidth=2, label='$profile-Main$')
	ax.plot(R, LambdaIP, 'g-', linewidth=2, label='$profile-Inter$')

	ax.set_xlabel('$r/R_{0}$')
	ax.set_ylabel('$J/J_{max}$')
	ax.legend()

	if form=='save':
		if not os.path.isdir('prof'):
			os.mkdir('prof')
		plt.savefig('prof/'+ 'prof_' + name['name'] + '_' + name['model_b'] + '.png')
	elif form=='print':
		plt.show()
	else:
		print("FormError")