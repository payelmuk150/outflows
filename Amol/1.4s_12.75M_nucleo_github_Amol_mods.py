import numpy as np
from math import *
import scipy as sc
from scipy import optimize
import mpmath as gp
import matplotlib.pyplot as plt
from time import process_time 
import bisect

t1_start = process_time()

# Amol: various physical constants and parameters
pi = np.pi
m_pl = 1.221e+22   		# Planck mass in MeV
hbarc = 1.973e-11		# hbarc in MeV cm
cl = 2.99798458e+10		# Speed of light in cm/s
M_sun = 1.988e+33		# Solar mass in g
MeVinerg = 1.602e-6	 	# 1 MeV in erg

g_ps = 2.0				# Relativistic degrees of freedom in the postshock region
Sbr_ps = 6.0 			# Entropy per baryon in the radiation field, in the postshock region

#strprog = 's18.0'
strprog = 's12.75'

# target pressure from Sukhbold progenitor
f_target = np.genfromtxt(strprog+'_presn', skip_header = 2)
#f_target = np.genfromtxt('s18.0_presn', skip_header = 2)
r_target = 0
M_target = 0
Rs = 1.2e+9				# Amol: Shock radius? This should move across the time sequence of snapshots, correct?
z_s = 0

filelen = np.shape(f_target)[0]
print('# filelen = ', filelen)

# Front shock parameters (to calculate the jump velocity)
vshock0 =  7.0e8 #1.0e9
tpb = 2.0


#tps = 1.4
#rshock0 = vshock0*tps
 
def vexp(r): 
	return vshock0*(r/Rs)


# Amol: Extract the radius and enclosed mass at the shock
for r in f_target[:, 2]:
	if r >= Rs:
		r_target = r 
		M_target = f_target[z_s, 1]
		break
	z_s += 1

pns_mass = 1.8			# PNS mass in solar masses

# Question from Amol: "the cell outer total mass" on the last line of the progenitor file is about 11.45 Msun. 
# Is it possible that they may have already subtracted the Chandrasekhar mass? In that case, do we need to subtract 1.8 Msun or ~0.4 Msun? 

rho_target = (M_target / M_sun - pns_mass) * M_sun / (4.0 * np.pi * r_target**3 / 3.0)
print('# Boundary rho', rho_target)

# Boundary pressure
P_target = (Sbr_ps * rho_target)**(4.0/3.0)		# P propto (S rho)^4/3, ignoring the constant factors
print('# Boundary pressure (target)', P_target)
# Amol: more accurate expression with constant factors given below:
#P_target = (1.0/4.0)*(45/(2*pi**2*g_ps))**(1.0/3.0) (Sbr_ps * rho_target_in_nat)^4
# Also, shouldn't rho be converted to natural units before using P \approx (S rho)^(4/3)?


#GM=7.56395e+15 * pns_mass
GM = (M_sun * cl**2)/(m_pl**2 * MeVinerg) * pns_mass	# G*M_PNS in MeV^{-1}

#r_front = 3.0e+9										# Is this the front shock radius? If so, then what is Rs?

# This is the radius of the termination shock. I have set it to a huge number here to make it redundant 
# in the context of subsonic outflows
#R_t = 1.8e+12

# PNS radius
R = 1.7 											# PNS radius in units of 10^6 cm
r_in = R * 1.0e+6									# units of cm
r_in_nat = r_in/hbarc								# Initial radius in MeV^{-1}

strmns = str(pns_mass) + "M"
strrns = str(round(r_in/1.0e5)) + "km"
strtpb = str(tpb) + "sGR"

print ("# PNS radius in km: ", r_in/1.0e5)
print ("# PNS mass in solar masses: ", strmns)
print ("# Initial time post bounce: ", strtpb)

# neutrino parameters. Luminosities and energies.

### Mixing radius
rmix = 1.701e6 											# rmix in cm # 1.0e13 # (infinity) # r_in/1.0e5 # 
rmix_in_nat = rmix/hbarc
strmix = str(round(rmix/1.0e5)) + "km"
if rmix > 1.0e9:
	strmix = "infkm"
print ("# Radius of flavor mixing: ", strmix)


### Bollig parameters 3D
strlum = "_Bollig_3D_lum" + "_" + strrns + "_" + strmns + "_" + strtpb + "_rmix=" + strmix
L_nue = 10.76										# units of 10^51 erg/s
L_nuebar = 10.04 
# Pinching parameters 2.1 for nue and 1.5 for nuebar
# Check below values using Yecalc script
e_nuebar = 21.8771 									# 16.23, sqrt(<E^3>/<E>), ,<E>_nuebar = Hudepohl_11.7MeV_1sec
e_nue = 17.464										# 13.08, sqrt(<E^3>/<E>), ,<E>_nue = Hudepohl_9.7MeV_1sec

# Bollig parameters 3D, with flavor mixing (equilibration) 
#strlum = "_Bollig_3D_lum_mix"
L_nue_mix = 11.34										# units of 10^51 erg/s
L_nuebar_mix = 11.10 
# Pinching parameters 2.1 for nue and 1.5 for nuebar
# Check below values using Yecalc script
e_nuebar_mix = 22.3347 									# 16.23, sqrt(<E^3>/<E>), ,<E>_nuebar = Hudepohl_11.7MeV_1sec
e_nue_mix = 21.0594										# 13.08, sqrt(<E^3>/<E>), ,<E>_nue = Hudepohl_9.7MeV_1sec

def L_nue_f(r):
	return L_nue*(r < rmix_in_nat) + L_nue_mix*(r >= rmix_in_nat)

def L_nuebar_f(r):
	return L_nuebar*(r < rmix_in_nat) + L_nuebar_mix*(r >= rmix_in_nat)

def e_nue_f(r):
	return e_nue*(r < rmix_in_nat) + e_nue_mix*(r >= rmix_in_nat)
	
def e_nuebar_f(r):
	return e_nuebar*(r < rmix_in_nat) + e_nuebar_mix*(r >= rmix_in_nat)
	
## Bollig parameters 1D
#strlum = "_Bollig_lum"
#L_nue = 7.69										# units of 10^51 erg/s
#L_nuebar = 7.37 
## Pinching parameters 2.1 for nue and 1.5 for nuebar
## Check below values using Yecalc script
#e_nuebar = 20.4733 									# 16.23, sqrt(<E^3>/<E>), ,<E>_nuebar = Hudepohl_11.7MeV_1sec
#e_nue = 15.9947										# 13.08, sqrt(<E^3>/<E>), ,<E>_nue = Hudepohl_9.7MeV_1sec


# Bollig parameters 1D, with flavor mixing (equilibration) 
#strlum = "_Bollig_lum_mix"
#L_nue = 8.83										# units of 10^51 erg/s
#L_nuebar = 8.72333 
## Pinching parameters 2.1 for nue and 1.5 for nuebar
## Check below values using Yecalc script
#e_nuebar = 20.6094 									# 16.23, sqrt(<E^3>/<E>), ,<E>_nuebar = Hudepohl_11.7MeV_1sec
#e_nue = 19.4089										# 13.08, sqrt(<E^3>/<E>), ,<E>_nue = Hudepohl_9.7MeV_1sec


## Bollig parameters 1D, with nux parameters in place of nue 
#strlum = "_Bollig_lum_nuxparams"
#L_nue = 9.4											# units of 10^51 erg/s
#L_nuebar = 8.72333 
## Pinching parameters 2.1 for nue and 1.5 for nuebar
## Check below values using Yecalc script
#e_nuebar = 20.6094 									# 16.23, sqrt(<E^3>/<E>), ,<E>_nuebar = Hudepohl_11.7MeV_1sec
#e_nue = 20.6531										# 13.08, sqrt(<E^3>/<E>), ,<E>_nue = Hudepohl_9.7MeV_1sec


## Hudepohl parameters 
#strlum = "_4.0xHudepohl_lum"
#L_nue = 4.0*6.12									# units of 10^51 erg/s
#L_nuebar = L_nue / 1.22 
## Pinching parameters 2.1 for nue and 1.5 for nuebar
## Check below values using Yecalc script
#e_nuebar = 16.23 									# 16.23, sqrt(<E^3>/<E>), ,<E>_nuebar = Hudepohl_11.7MeV_1sec
#e_nue = 13.08										# 13.08, sqrt(<E^3>/<E>), ,<E>_nue = Hudepohl_9.7MeV_1sec

# Initial entropy 
S_in= 6.0											# What is this? Is this entropy at the inner boundary of gain region?

# masses in MeV
m_n = 939.57
m_e = 0.511

def g1(x_th, y, z, w, T):
	x = np.arange(x_th, x_th + 20, 0.1)
	integrand_100 = x**2 * np.sqrt(np.abs(x**2 - x_th**2)) / (np.exp(x + y) + 1)
	g100 = np.trapz(integrand_100, x)
	return g100

def g2(x_th, y, z, w, T):
	x = np.arange(x_th, x_th + 20, 0.1)
	integrand_200 = x**3.0 * np.sqrt(np.abs(x**2.0 - x_th**2.0)) / (np.exp(x + y) + 1)
	g200 = np.trapz(integrand_200, x)
	return g200

def initial_T(T):
#	set R_nu ~ R_gain

#	r_in = (R*1.0e+6)
#	r_nu_nat=r_in*(5.0e+10)
	r_nu_nat=r_in/hbarc					# Neutrinosphere radius in Mev^{-1}	
	r_in_nat=r_nu_nat

	GR_factor_angle = ((1.0 - (2.0 * GM / (r_nu_nat))) / (1.0 - (2.0 * GM / r_in_nat)))**0.5
	factor = 1.0 - (1.0 - ((r_nu_nat / r_in_nat)**2.0 / GR_factor_angle**2.0))**0.5

	rho = (T**3.0 / S_in) * 5.21 * 1e+8 					# Density in g/cm^3 using Qian and Woosley Eq. (9)
#	rho_MeV = rho * 5.6e+26 * (1.97e-11)**3.0				# Energy density in MeV^{-4}
	rho_MeV = rho * cl**2/MeVinerg * hbarc**3.0				# Energy density in MeV^{-4}
	eta = 3.0 * 0.5 * (rho_MeV / m_n) / T**3.0				# Y_e = 0.5 assumed?
		
	qdot_cool_ann = (1.57e-25 * T**9.0 * 
					(g1(m_e / T, eta, 0, 0, T) * g2(m_e / T, - eta, 0, 0, T) + g1(m_e / T, -eta, 0, 0, T) * g2(m_e / T,eta, 0, 0, T)) / rho_MeV)

	return ((3.4e-24 * ((L_nuebar_f(r_in_nat) * e_nuebar_f(r_in_nat)**2.0 + L_nue_f(r_in_nat) * e_nue_f(r_in_nat)**2.0) * GR_factor_1**6.0 * (1.0e+6/r_in)**2) * factor * GR_factor_angle**6.0
			- 1.6e-24 * T**6.0) - qdot_cool_ann)
				   # Includes heating by neutrino absorption, cooling by e-/e+ capture, and cooling by e-e+ annihilation
				   # Generalization of Qian and Woosley, Eqs. (10), (11), (13), and (9) with GR corrections and changing g_star

def fermion_func(x):
	return - (x * (2.0 * x * m_e**2.0 * float(gp.fp.polylog(2, - np.exp(-m_e / x))) + 6.0 * x**2.0 * m_e * float(gp.fp.polylog(3, - np.exp(-m_e / x))) 
			+ 6.0 * x**3.0 * float(gp.fp.polylog(4, - np.exp(-m_e / x)) )) * 4.0 * 15.0 / (np.pi**4.0 * x**4.0))

def P_func(x):
	return (2.0 + fermion_func(x))*pi**2 / 90.0 * x**4

print('# gstar at 3 GK', 2 + fermion_func(0.0617))			# why 0.0617?

# GR factor (initial)
GR_factor_1 = np.sqrt(1.0 / (1.0 - (2.0 * GM / r_in_nat)))
print('# GR_factor', GR_factor_1)

# Initial temperature
T_in = sc.optimize.fsolve(initial_T, 4.0)[0]			# Solve for intial temperature by setting initial_T(x) = 0, i.e., heating rate = cooling rate
print('# T_in',T_in)									# T_in in MeV?

# Initial density
rho_start = 5.21 * T_in**3.0 * 1e+8 / S_in 											# g/cm^3, Qian and Woosley Eq. (9)
print ('# Starting density', rho_start)

# calculate the initial derivate of g_star (logarithmic derivative?)
fermion = fermion_func(T_in)
A = 2.0 + fermion
print('# g_star at T_in = ', A, "\n")
diff_T = 0.0001
fermion_dT = fermion_func(T_in + diff_T)
fermion_derivative = (fermion_dT - fermion) / (diff_T * A)
derivative_A_term = fermion_derivative

# Initial sound speed = sqrt(dP/d\rho)?
vs_in = ((T_in * S_in / (4.0 * m_n)) * (4.0 + (T_in * derivative_A_term)) / (3.0 + (T_in * derivative_A_term)))**0.5
#print ('# Initial sound speed (cm/s)', vs_in * cl)

rho_in = (T_in**3.0 / S_in) * A * 1e+8 / 1.055

## Set radial increment (make this adaptive eventually)
#dr = 1.0e+4 								# dr in cm 
##	dr_nat = dr * 5.0e+10						# dr in MeV^{-1}
#dr_nat = dr/hbarc							# dr in MeV^{-1}

#initial_vs = [(2.00e+6 + 2 * 1e-4 * 1e+6 * i) for i in range(1)]
#initial_vs = [0.1*vs_in*cl]

viter = 0

v_in_lo = 0.0*cl
v_in_hi = vs_in*cl

v_in = (v_in_lo + v_in_hi)/2.0

convergence_fraction = 1.0

print ('# v_in rel error = ', (v_in_hi-v_in_lo)/(vs_in*cl))

#while viter < 30:
while (v_in_hi-v_in_lo)/(vs_in*cl) > 1.0e-9 or v_in_lo != v_in:

	v_in = (v_in_lo + v_in_hi)/2.0

#for v_in in initial_vs:
	viter += 1

#	f = open('Nucleo_Amol.txt' ,'w')
	z = z_s

	# initial time
	t=0.0
#	dt=0.0

	# to calculate accumulated mass
	M = 0.0
#	dM = 0.0

	# set r = initial radius (MeV^{-1})
	r = r_in_nat

	#	Reset sound speed = initial sound speed at the start of each iteration
	vs = vs_in

	printed = False
	
	print ('# Iteration # ', viter)

	print ('# Initial sound speed (cm/s)', vs_in * cl)
	print ('# Initial speed (cm/s)', v_in)
	v_nat_in = v_in/cl														# Initial speed in units of c
	v = v_nat_in
	S = S_in
	T = T_in
	
#	dt = (dr_nat * hbarc) / (v * cl)							# dt in s

	v_list = [] # [v_in]
	r_list = [] # [r]
	t_list = [] #[dt]
	v_list = []
	vexp_list = []
	vs_list = []
	S_list = []
	T_list = [] #[T_in]
	rho_list = [] #[rho_in]
	stage = []
	qdot_list = []
	qdot_cool_list = []
#	dT_list = []
#	dT_list2 = []
	Tferr_list = []
	
	Tferr = GM/(r * 3.0 * vs**2)
	
	# mass outflow rate
	M_dot = 4.0 * np.pi * r_in**2.0 * rho_start * v_in / (M_sun * GR_factor_1 ) 		# in solar masses
	print('# Mass outflow rate' , M_dot)
	
	# Initial qdot set to zero at the gain radius.
	qdot = 0.0

	rsteps = 0
	
	while T > 0 and S > 0 and vs > 0 and not np.isnan(T) and rsteps < 100000 and v < vs and Tferr > 0.025: # and r*hbarc < R_t
	
		rsteps += 1

		# Set radial increment (make this adaptive eventually)
		dr = 1.0e+3*(r/r_in_nat) 					# dr in cm (logarithmically spaced intervals) - initially 10^3 cm, and increases with distance
		#	dr_nat = dr * 5.0e+10					# dr in MeV^{-1}
		dr_nat = dr/hbarc							# dr in MeV^{-1}

		fermion = fermion_func(T)
		A = 2.0 + fermion

		GR_factor = 1.0 / (1.0 - (2.0 * GM / r))
		GR_factor_angle = ((1.0 - (2.0 * GM / r_in_nat)) / (1.0 - (2.0 * GM / r)))**0.5
		y_fac = ((1.0 - (2.0 * GM / r)) / (1.0 - v**2))**0.5
	   
		density8 = (T**3.0 / S) * A / 1.055												# Is this a more general expression than Q&W Eq. (9)?
		X_n = 828 * T**(9.0 / 8) * np.exp(- 7.074 / T) / density8**(3.0 / 4.0)			# Neutron fraction?

		fac = 1.0 if X_n >= 1.0 else X_n	

		dt = (dr_nat * hbarc) / (v * cl)							# dt in s
		t = t + dt

		fermion = fermion_func(T)
		A = 2.0 + fermion
		diff_T = 0.0001
		fermion_dT = fermion_func(T + diff_T)
		fermion_derivative = (fermion_dT - fermion) / ( diff_T * A)
		derivative_A_term = fermion_derivative 

		beta = (1.0 / 4.0) * (1.0 + (3.0 + T * derivative_A_term)**-1)
		dv = ((2 * vs**2 / r) - ((GM / r**2) * (1 - vs**2) * GR_factor) - (qdot * beta / (v * y_fac * (1 + 3 * vs**2)))) * dr_nat * (1.0 - v**2.0) / (v - (vs * vs / v))
		dS = (qdot * m_n / (T * v * y_fac)) * dr_nat
		dvs = ((qdot * dr_nat / (v * y_fac * (1 + 3 * vs**2))) - (v * dv / (1 - v**2)) - (GM * GR_factor * dr_nat / r**2)) * (1 + 3 * vs**2) / (6 * vs)

		S = S + dS
		r = r + dr_nat
		v = v + dv
		vs = vs + dvs

		a = derivative_A_term 
		b = 4 - ((4 * m_n * vs**2 / S) * derivative_A_term)
		c = -12 * m_n * vs**2 / S

#		Told = T
		#T = (vs**2 * 4 * m_n / S) / (1.0 + (3.0 + T * derivative_A_term)**-1)
		T = (-b + (b**2 - 4*a*c)**0.5) / (2 * a)

#		dlogTdlogr = (Told - T) / Told * r / dr_nat

		rho = (T**3.0 / S) * A * 1.0e+8 / 1.055								# rho in g/cm^3
#		dM = 4 * pi * rho * (r / 5.0e+10)**2.0 * (dr_nat / 5.0e+10)
		dM = 4 * pi * rho * (r * hbarc)**2.0 * (dr_nat * hbarc)
		M = M + dM

		Tferr = GM/(r * 3.0 * vs**2)
		
#		if (T<= 0.86):
#			r_list.append(r / 5.0e+10)
		r_list.append(r * hbarc)						# r_list in cm
		t_list.append(t)								# t_list in s
		T_list.append(T)								# T_list in MeV
		rho_list.append(rho)							# rho_list in g/cm**3
		v_list.append(v*cl)								# v_list in cm/s
		vs_list.append(vs*cl)							# vs_list in cm/s
		S_list.append(S)								
		Tferr_list.append(Tferr)
		
#		dT_list.append(dlogTdlogr)
#		dT_list2.append(GM / (r * 3.0 * vs**2))
		
		stage.append("outflow")

		rho_8 = rho / 1e+8
#		rho_MeV = rho * 5.6e+26 * (1.97e-11)**3.0
		rho_MeV = rho * 5.6e+26 * hbarc**3.0
		eta = 3.0 * 0.5 * (rho_MeV / m_n) / T**3.0
		
		qdot_cool_ann = (1.57e-25 * T**9.0 * 
						(g1(m_e / T, eta, 0, 0, T) * g2(m_e / T, - eta, 0, 0, T) + g1(m_e / T, -eta, 0, 0, T) * g2(m_e / T,eta, 0, 0, T)) / rho_MeV)

		factor = 1.0 - (1.0 - ((r_in_nat / r)**2.0 / GR_factor_angle**2.0))**0.5 # gravitational bending
		qdot = ((3.4e-24 * ((L_nuebar_f(r) * e_nuebar_f(r)**2.0 + L_nue_f(r) * e_nue_f(r)**2.0) * GR_factor_1**6.0 * (1.0e+6 / r_in)**2 ) * factor * GR_factor_angle**6.0
			- 1.6e-24 * T**6.0) * fac - qdot_cool_ann)

		qdot_list.append(qdot)
		qdot_cool_list.append(qdot_cool_ann + 1.6e-24 * T**6.0 * fac)

#	Pressure-matching at the far boundary
	P_calc = (S * (T**3.0 / S) * A * 1e+8 / 1.055)**(4.0/3.0) # P propto (S rho)^4/3, ignoring the constant factors
	convergence_fraction = abs(P_calc - P_target) / P_target

	print('# Calculated far boundary S: ', S)
	print('# Calculated far boundary rho: ', (T**3.0 / S) * A * 1e+8)
	print('# Calculated far pressure: ', P_calc)
	print('# Target far pressure: ', P_target)
	print('# fractional convergence: ', convergence_fraction)

#	f.close()

	if (vs < 0):
		print ("# Negative sound speed encountered", vs, "\n")
		v_in_lo = v_in
		continue

	elif (v > vs):
		print("# Outflow speed exceeds sound speed", v, vs, "\n")
		v_in_hi = v_in
		continue
		
	elif (T < 0):
		print ('# Negative temperature encountered', T, "\n")
		v_in_lo = v_in
		continue
		
	elif (S < 0):
		print ("# Negative entropy encountered", S, "\n")
		v_in_lo = v_in
		continue
		
	elif (np.isnan(T)):
		print('# T is nan', "\n")
		break

#			break
	if (convergence_fraction < 0.1):
		print("# Far pressure match achieved (at 10% level)", P_calc, P_target, "\n")
		break
		
	elif (P_calc < P_target):
		print("# Calculated pressure less than target pressure", P_calc, P_target, "\n")
		v_in_hi = v_in
		continue

	elif (P_calc > P_target):
		print("# Calculated pressure exceeds target pressure", P_calc, P_target, "\n")
		v_in_lo = v_in
		continue	

#############################################################################################################################################
################################### 															#############################################
###################################  Extra calculations in case subsonic solution is not found  #############################################
################################### 															#############################################
#############################################################################################################################################

###################################### Jump across the sonic point ##############################################
# if viter == 30:
if (v_in_hi-v_in_lo)/(vs_in*cl) <= 1.0e-9 and v_in_lo == v_in:
	icrit = np.argmax(v_list)
	vcrit = np.max(v_list)
	rcrit = r_list[icrit]
	Scrit = S_list[icrit]
	vscrit = vs_list[icrit]
	
	print("# vcrit = ", vcrit)
	print("# vscrit = ", vscrit)
	print("# rcrit = ", rcrit, "\n")
	
###################################### Run another iteration, but with a long jump over the sonic point #######################

#	f = open('Nucleo_Amol.txt' ,'w')
	z = z_s

	# initial time
	t=0.0
#	dt=0.0

	# to calculate accumulated mass
	M = 0.0
#	dM = 0.0

	# set r = initial radius (MeV^{-1})
	r = r_in_nat

	#	Reset sound speed = initial sound speed at the start of each iteration
	vs = vs_in

	printed = False
	
	print ('# Iteration # ', viter+1)

	print ('# Initial sound speed (cm/s)', vs_in * cl)
	print ('# Initial speed (cm/s)', v_in)
	v_nat_in = v_in/cl														# Initial speed in units of c
	v = v_nat_in
	S = S_in
	T = T_in
	
#	dt = (dr_nat * hbarc) / (v * cl)							# dt in s

	v_list = [] # [v_in]
	r_list = [] # [r]
	t_list = [] #[dt]
	v_list = []
	vexp_list = []
	vs_list = []
	S_list = []
	T_list = [] #[T_in]
	rho_list = [] #[rho_in]
	stage = []
	qdot_list = []
#	dT_list = []
#	dT_list2 = []
	Tferr_list = []
	
	Tferr = GM/(r * 3.0 * vs**2)
	
	# mass outflow rate
	M_dot = 4.0 * np.pi * r_in**2.0 * rho_start * v_in / (M_sun * GR_factor_1 ) 		# in solar masses
	print('# Mass outflow rate' , M_dot)
	
	# Initial qdot set to zero at the gain radius.
	qdot = 0.0

	rsteps = 0
	
	while T > 0 and S > 0 and vs > 0 and not np.isnan(T) and rsteps < 100000 and Tferr > 0.025:
	
		rsteps += 1

		# Set radial increment (make this adaptive eventually)
		dr = 1.0e+3*(r/r_in_nat) 					# dr in cm (logarithmically spaced intervals) - initially 10^3 cm, and increases with distance
	
		#### Take longer step if near the critical point #####
		
		if (v*cl > 0.99*vcrit and v*cl < 1.01*vcrit):
			dr = 2.0*(rcrit - r*hbarc)

		#	dr_nat = dr * 5.0e+10					# dr in MeV^{-1}
		dr_nat = dr/hbarc							# dr in MeV^{-1}

		fermion = fermion_func(T)
		A = 2.0 + fermion

		GR_factor = 1.0 / (1.0 - (2.0 * GM / r))
		GR_factor_angle = ((1.0 - (2.0 * GM / r_in_nat)) / (1.0 - (2.0 * GM / r)))**0.5
		y_fac = ((1.0 - (2.0 * GM / r)) / (1.0 - v**2))**0.5
	   
		density8 = (T**3.0 / S) * A / 1.055												# Is this a more general expression than Q&W Eq. (9)?
		X_n = 828 * T**(9.0 / 8) * np.exp(- 7.074 / T) / density8**(3.0 / 4.0)			# Neutron fraction?

		fac = 1.0 if X_n >= 1.0 else X_n

		dt = (dr_nat * hbarc) / (v * cl)							# dt in s
		t = t + dt

		fermion = fermion_func(T)
		A = 2.0 + fermion
		diff_T = 0.0001
		fermion_dT = fermion_func(T + diff_T)
		fermion_derivative = (fermion_dT - fermion) / ( diff_T * A)
		derivative_A_term = fermion_derivative 

		beta = (1.0 / 4.0) * (1.0 + (3.0 + T * derivative_A_term)**-1)
		dv = ((2 * vs**2 / r) - ((GM / r**2) * (1 - vs**2) * GR_factor) - (qdot * beta / (v * y_fac * (1 + 3 * vs**2)))) * dr_nat * (1.0 - v**2.0) / (v - (vs * vs / v))
		dS = (qdot * m_n / (T * v * y_fac)) * dr_nat
		dvs = ((qdot * dr_nat / (v * y_fac * (1 + 3 * vs**2))) - (v * dv / (1 - v**2)) - (GM * GR_factor * dr_nat / r**2)) * (1 + 3 * vs**2) / (6 * vs)

		S = S + dS
		r = r + dr_nat
		v = v + dv
		vs = vs + dvs

		a = derivative_A_term 
		b = 4 - ((4 * m_n * vs**2 / S) * derivative_A_term)
		c = -12 * m_n * vs**2 / S

#		Told = T
		#T = (vs**2 * 4 * m_n / S) / (1.0 + (3.0 + T * derivative_A_term)**-1)
		T = (-b + (b**2 - 4*a*c)**0.5) / (2 * a)

		if T <= 0:
			print ('# Non-positive temperature encountered: ', T, "\n")
			break
		
#		dlogTdlogr = (Told - T) / Told * r / dr_nat

		rho = (T**3.0 / S) * A * 1.0e+8 / 1.055								# rho in g/cm^3
#		dM = 4 * pi * rho * (r / 5.0e+10)**2.0 * (dr_nat / 5.0e+10)
		dM = 4 * pi * rho * (r * hbarc)**2.0 * (dr_nat * hbarc)
		M = M + dM

		Tferr = GM/(r * 3.0 * vs**2)
		
#		if (T<= 0.86):
#			r_list.append(r / 5.0e+10)
		r_list.append(r * hbarc)						# r_list in cm
		t_list.append(t)								# t_list in s
		T_list.append(T)								# T_list in MeV
		rho_list.append(rho)							# rho_list in g/cm**3
		v_list.append(v*cl)								# v_list in cm/s
		vs_list.append(vs*cl)							# vs_list in cm/s
		S_list.append(S)								
		Tferr_list.append(Tferr)
		
#		dT_list.append(dlogTdlogr)
#		dT_list2.append(GM / (r * 3.0 * vs**2))
		
		stage.append("outflow")

#		print(t, r * hbarc, T, rho, v*cl, vs*cl, S, Tferr, a)

		rho_8 = rho / 1e+8
#		rho_MeV = rho * 5.6e+26 * (1.97e-11)**3.0
		rho_MeV = rho * 5.6e+26 * hbarc**3.0
		eta = 3.0 * 0.5 * (rho_MeV / m_n) / T**3.0
		
		qdot_cool_ann = (1.57e-25 * T**9.0 * 
						(g1(m_e / T, eta, 0, 0, T) * g2(m_e / T, - eta, 0, 0, T) + g1(m_e / T, -eta, 0, 0, T) * g2(m_e / T,eta, 0, 0, T)) / rho_MeV)

		factor = 1.0 - (1.0 - ((r_in_nat / r)**2.0 / GR_factor_angle**2.0))**0.5 # gravitational bending
		qdot = ((3.4e-24 * ((L_nuebar_f(r) * e_nuebar_f(r)**2.0 + L_nue_f(r) * e_nue_f(r)**2.0) * GR_factor_1**6.0 * (1.0e+6 / r_in)**2 ) * factor * GR_factor_angle**6.0
			- 1.6e-24 * T**6.0) * fac - qdot_cool_ann)

		qdot_list.append(qdot)

	R_t_lo = rcrit
	R_t_hi = r * hbarc
	
	rtiter = 0
	convergence_fraction = 1.0
	
	while rtiter < 30:	
		
		rtiter += 1
		R_t = (R_t_lo + R_t_hi)/2.0		
	#	R_t = 0.8e8																# Wind-termination radius in cm
	
		print ('# Termination shock Iteration # ', rtiter)
		print ('# Termination shock radius [km] = ', R_t/1.0e5)
		
		idxRt = bisect.bisect(r_list, R_t)
		
		r = r_list[idxRt]/hbarc
		t = t_list[idxRt]
		v = v_list[idxRt]/cl
		vs = vs_list[idxRt]/cl
		T = T_list[idxRt]
		S = S_list[idxRt]
		rho = rho_list[idxRt]
		Tferr = Tferr_list[idxRt]
		qdot = qdot_list[idxRt]
		
		P = (2.0 + fermion_func(T))*pi**2 / 90.0 * T**4
		
		print(v,vs,T,S,rho,P)
			
		rsteps = 0
		
		# Shock jump conditions
		rho2 = rho * 7.0*v**2 / (v**2 + 6.0*vs**2)														# Density
		v2 = v * rho / rho2  																			# Velocity
		vs2 = np.sqrt(v**2/6.0 + vs**2 - v2**2/6.0)														# Sound speed
		P2 = (rho*v**2 - rho2*v2**2)*cl**2/MeVinerg*hbarc**3 + P_func(T)								# Pressure
		
		T2 = sc.optimize.fsolve(lambda x: P_func(x) - P2, T)[0]
		
		S2 = 3.0*vs2**2 * m_n/T2
		
		print(v2,vs2,T2,S2,rho2,P2)
		
		rho = rho2
		v = v2
		T = T2
		S = S2
		vs = vs2
		P = P2
		
		v2_list = [] # [v_in]
		r2_list = [] # [r]
		t2_list = [] #[dt]
		v2_list = []
		vexp2_list = []
		vs2_list = []
		S2_list = []
		T2_list = [] #[T_in]
		rho2_list = [] #[rho_in]
		stage2 = []
		qdot2_list = []
		Tferr2_list = []

		rsteps = 0
		
		while T > 0 and S > 0 and vs > 0 and not np.isnan(T) and rsteps < 100000 and Tferr > 0.025:
			
			rsteps += 1

			# Set radial increment (make this adaptive eventually)
			dr = 1.0e+3*(r/r_in_nat) 					# dr in cm (logarithmically spaced intervals) - initially 10^3 cm, and increases with distance
		
			#### Take longer step if near the critical point #####
			
			#	dr_nat = dr * 5.0e+10					# dr in MeV^{-1}
			dr_nat = dr/hbarc							# dr in MeV^{-1}

			fermion = fermion_func(T)
			A = 2.0 + fermion

			GR_factor = 1.0 / (1.0 - (2.0 * GM / r))
			GR_factor_angle = ((1.0 - (2.0 * GM / r_in_nat)) / (1.0 - (2.0 * GM / r)))**0.5
			y_fac = ((1.0 - (2.0 * GM / r)) / (1.0 - v**2))**0.5
		   
			density8 = (T**3.0 / S) * A / 1.055												# Is this a more general expression than Q&W Eq. (9)?
			X_n = 828 * T**(9.0 / 8) * np.exp(- 7.074 / T) / density8**(3.0 / 4.0)			# Neutron fraction?

			fac = 1.0 if X_n >= 1.0 else X_n

			dt = (dr_nat * hbarc) / (v * cl)							# dt in s
			t = t + dt

			fermion = fermion_func(T)
			A = 2.0 + fermion
			diff_T = 0.0001
			fermion_dT = fermion_func(T + diff_T)
			fermion_derivative = (fermion_dT - fermion) / ( diff_T * A)
			derivative_A_term = fermion_derivative 

			beta = (1.0 / 4.0) * (1.0 + (3.0 + T * derivative_A_term)**-1)
			dv = ((2 * vs**2 / r) - ((GM / r**2) * (1 - vs**2) * GR_factor) - (qdot * beta / (v * y_fac * (1 + 3 * vs**2)))) * dr_nat * (1.0 - v**2.0) / (v - (vs * vs / v))
			dS = (qdot * m_n / (T * v * y_fac)) * dr_nat
			dvs = ((qdot * dr_nat / (v * y_fac * (1 + 3 * vs**2))) - (v * dv / (1 - v**2)) - (GM * GR_factor * dr_nat / r**2)) * (1 + 3 * vs**2) / (6 * vs)

			S = S + dS
			r = r + dr_nat
			v = v + dv
			vs = vs + dvs

			a = derivative_A_term 
			b = 4 - ((4 * m_n * vs**2 / S) * derivative_A_term)
			c = -12 * m_n * vs**2 / S

	#		Told = T
			#T = (vs**2 * 4 * m_n / S) / (1.0 + (3.0 + T * derivative_A_term)**-1)
			T = (-b + (b**2 - 4*a*c)**0.5) / (2 * a)

			if T <= 0:
				print ('# Non-positive temperature encountered: ', T, "\n")
				break
			
	#		dlogTdlogr = (Told - T) / Told * r / dr_nat

			rho = (T**3.0 / S) * A * 1.0e+8 / 1.055								# rho in g/cm^3
	#		dM = 4 * pi * rho * (r / 5.0e+10)**2.0 * (dr_nat / 5.0e+10)
			dM = 4 * pi * rho * (r * hbarc)**2.0 * (dr_nat * hbarc)
			M = M + dM

			Tferr = GM/(r * 3.0 * vs**2)
			
	#		if (T<= 0.86):
	#			r_list.append(r / 5.0e+10)
			r2_list.append(r * hbarc)						# r_list in cm
			t2_list.append(t)								# t_list in s
			T2_list.append(T)								# T_list in MeV
			rho2_list.append(rho)							# rho_list in g/cm**3
			v2_list.append(v*cl)								# v_list in cm/s
			vs2_list.append(vs*cl)							# vs_list in cm/s
			S2_list.append(S)								
			Tferr2_list.append(Tferr)
			
	#		dT_list.append(dlogTdlogr)
	#		dT_list2.append(GM / (r * 3.0 * vs**2))
			
			stage.append("outflow")

	#		print(t, r * hbarc, T, rho, v*cl, vs*cl, S, Tferr, a)

			rho_8 = rho / 1e+8
	#		rho_MeV = rho * 5.6e+26 * (1.97e-11)**3.0
			rho_MeV = rho * 5.6e+26 * hbarc**3.0
			eta = 3.0 * 0.5 * (rho_MeV / m_n) / T**3.0
			
			qdot_cool_ann = (1.57e-25 * T**9.0 * 
							(g1(m_e / T, eta, 0, 0, T) * g2(m_e / T, - eta, 0, 0, T) + g1(m_e / T, -eta, 0, 0, T) * g2(m_e / T,eta, 0, 0, T)) / rho_MeV)

			factor = 1.0 - (1.0 - ((r_in_nat / r)**2.0 / GR_factor_angle**2.0))**0.5 # gravitational bending
			qdot = ((3.4e-24 * ((L_nuebar_f(r) * e_nuebar_f(r)**2.0 + L_nue_f(r) * e_nue_f(r)**2.0) * GR_factor_1**6.0 * (1.0e+6 / r_in)**2 ) * factor * GR_factor_angle**6.0
				- 1.6e-24 * T**6.0) * fac - qdot_cool_ann)

			qdot2_list.append(qdot)
			
	##	Pressure-matching at the far boundary
		P_calc = (S * (T**3.0 / S) * A * 1e+8 / 1.055)**(4.0/3.0) # P propto (S rho)^4/3, ignoring the constant factors
		convergence_fraction = abs(P_calc - P_target) / P_target

		print('# Calculated far boundary S: ', S)
		print('# Calculated far boundary rho: ', (T**3.0 / S) * A * 1e+8)
		print('# Calculated far pressure: ', P_calc)
		print('# Target far pressure: ', P_target)
		print('# fractional convergence: ', convergence_fraction)


		if (convergence_fraction < 0.1):
			print("# Far pressure match achieved (at 10% level)", P_calc, P_target, "\n")
			break
			
		elif (P_calc < P_target):
			print("# Calculated pressure less than target pressure", P_calc, P_target, "\n")
			R_t_hi = R_t
			continue

		elif (P_calc > P_target):
			print("# Calculated pressure exceeds target pressure", P_calc, P_target, "\n")
			R_t_lo = R_t
			continue	


#################################################################################################################
####################################### End supersonic calculations #############################################
#################################################################################################################



#################################################################################################################
################################ Blend into homologous expansion (with plowing) #################################
#################################################################################################################

#print(len(t_list), len(r_list), len(T_list), len(rho_list))
if viter < 30:
	f = open('vin_tests_'+strprog+strlum+'.txt' ,'w')
	for i in range(len(t_list)):
		f.write(f'{r_list[i]} {T_list[i]} {S_list[i]} {v_list[i]} {vs_list[i]} {Tferr_list[i]}\n')
	f.close()

	f = open('qdot_test_'+strprog+strlum+'.txt' ,'w')
	for i in range(len(qdot_list)):
		f.write(f'{r_list[i]} {qdot_list[i]} {qdot_cool_list[i]} \n')
	f.close()

	f = open('Nucleo_'+strprog+strlum+'.txt' ,'w')
	for i in range(len(t_list)):
		vexp_list.append(vexp(r_list[i]))
		f.write(f'{t_list[i]} {r_list[i]} {T_list[i]} {rho_list[i]} {v_list[i]}\n')
		if (vexp_list[i] > v_list[i]):
			T_break = T_list[i]
	#		rho_break = (T**3.0 / S) * A *1e+8 / 1.055				# rho_break in g/cm^3
			rho_break = rho_list[i]
			time_break = t_list[i]
			v_break = v_list[i]										# v_break in cm/s
			r_break = r_list[i]										# r_break in cm
	#		A_break = A
			r_hubble = r_break										# r_hubble in cm
			rho_ex = rho_break
			T_ex = T_break	
			print("# Break:", time_break, r_break, T_break, rho_break, v_break, "\n")
			break;
		

	while z < filelen-1: 												# Expansion based on plowed mass? But why stop at 1160 specifically?
		dt_break = (f_target[z+1, 2] - f_target[z, 2]) / vshock0	# dt_break in s
		time_break += dt_break
	#		r_hubble += v_break * 3e+10 * dt_break
		r_hubble += v_break * dt_break
	#		rho_factor = (f_target[z + 1, 1] - pns_mass * 1.98e+33) * f_target[z, 2]**3 / ((f_target[z, 1] - pns_mass * 1.98e+33) * f_target[z + 1, 2]**3)
		rho_factor = (f_target[z + 1, 1] - pns_mass * M_sun) * f_target[z, 2]**3 / ((f_target[z, 1] - pns_mass * M_sun) * f_target[z + 1, 2]**3)
		rho_ex *= rho_factor
	#		T_ex *= rho_factor**0.3333
		T_ex *= rho_factor**(1.0/3.0)
		z+=1
		t_expansion = time_break
		
		f.write(f'{time_break} {r_hubble} {T_ex} {rho_ex} {v_break}\n')
			
	f.close()

##				print(r*hbarc,S,v,vs,T)
else:
	f = open('vin_tests_'+strprog+strlum+'.txt' ,'w')
	for i in range(idxRt):
		f.write(f'{r_list[i]} {T_list[i]} {S_list[i]} {v_list[i]} {vs_list[i]} {Tferr_list[i]}\n')
	for i in range(len(t2_list)):
		f.write(f'{r2_list[i]} {T2_list[i]} {S2_list[i]} {v2_list[i]} {vs2_list[i]} {Tferr2_list[i]}\n')
	f.close()
	
	f = open('Nucleo_'+strprog+strlum+'.txt' ,'w')
	for i in range(idxRt):
		vexp_list.append(vexp(r_list[i]))
		f.write(f'{t_list[i]} {r_list[i]} {T_list[i]} {rho_list[i]} {v_list[i]}\n')
	for i in range(len(t2_list)):
		vexp_list.append(vexp(r2_list[i]))
		f.write(f'{t2_list[i]} {r2_list[i]} {T2_list[i]} {rho2_list[i]} {v2_list[i]}\n')
		if (vexp_list[i+idxRt] > v2_list[i]):
			T_break = T2_list[i]
	#		rho_break = (T**3.0 / S) * A *1e+8 / 1.055				# rho_break in g/cm^3
			rho_break = rho2_list[i]
			time_break = t2_list[i]
			v_break = v2_list[i]										# v_break in cm/s
			r_break = r2_list[i]										# r_break in cm
	#		A_break = A
			r_hubble = r_break										# r_hubble in cm
			rho_ex = rho_break
			T_ex = T_break	
			print("# Break:", time_break, r_break, T_break, rho_break, v_break, "\n")
			break;

	while z < filelen-1: 												# Expansion based on plowed mass? But why stop at 1160 specifically?
		dt_break = (f_target[z+1, 2] - f_target[z, 2]) / vshock0	# dt_break in s
		time_break += dt_break
	#		r_hubble += v_break * 3e+10 * dt_break
		r_hubble += v_break * dt_break
	#		rho_factor = (f_target[z + 1, 1] - pns_mass * 1.98e+33) * f_target[z, 2]**3 / ((f_target[z, 1] - pns_mass * 1.98e+33) * f_target[z + 1, 2]**3)
		rho_factor = (f_target[z + 1, 1] - pns_mass * M_sun) * f_target[z, 2]**3 / ((f_target[z, 1] - pns_mass * M_sun) * f_target[z + 1, 2]**3)
		rho_ex *= rho_factor
	#		T_ex *= rho_factor**0.3333
		T_ex *= rho_factor**(1.0/3.0)
		z+=1
		t_expansion = time_break
		
		f.write(f'{time_break} {r_hubble} {T_ex} {rho_ex} {v_break}\n')


print('# S', S)
print('# rho_f',(T**3.0/S)*A*1e+8 / 1.055)

P_calc = (S * (pow(T,3.0)/S)*A*1e+8 / 1.055)**(4.0/3.0) # P propto (S rho)^4/3, ignoring the constant factors
print('# Calculated boundary pressure', P_calc)

t1_stop = process_time()
print('Elapsed time = ', t1_stop - t1_start)

