import numpy as np
from math import *
import scipy as sc
from scipy import optimize
import mpmath as gp
import matplotlib.pyplot as plt

# target pressure from Sukhbold 12.75 M parameter
f_target = np.genfromtxt('s12.75_presn', skip_header = 2)
r_target = 0
M_target = 0
Rs = 1.2e+9
z_s = 0

for r in f_target[:, 2]:
    if r >= Rs:
        r_target = r 
        M_target = f_target[z_s, 1]
        
        break
    z_s += 1

print('z_s', z_s)
pns_mass = 1.8

rho_target = (M_target / 1.98e+33 - pns_mass) * 1.98e+33 / (4 * np.pi * r_target**3 / 3)
print('Boundary rho', rho_target)

# Boundary pressure
P_target = (6 * rho_target)**(4.0/3) # P propto (S rho)^4/3, ignorong the constant factors
print('Boundary pressure', P_target)

pns_mass = 1.8

GM=7.56395e+15 * pns_mass

# Termination Shock Radius. Set to a large value for the first part

R_t = 0.6e+13

# Mixed neutrino parameters. Luminosities and energies.
L_nue = 11.34
L_nuebar = 11.1 # 10^51 erg/s
R = 1.7 # 19 km
r_in= R * 1.0e+6 

e_nue = 21.0594    # 13.08, sqrt(<E^3>/<E>), ,<E>_nue = Hudepohl_9.7MeV_1sec
e_nuebar = 22.3347 # 16.23, sqrt(<E^3>/<E>), ,<E>_nuebar = Hudepohl_11.7MeV_1sec

eps_nue = 18.9042    # 13.08, sqrt(<E^3>/<E>), ,<E>_nue = Hudepohl_9.7MeV_1sec
eps_nuebar = 20.1877 # 16.23, sqrt(<E^3>/<E>), ,<E>_nuebar = Hudepohl_11.7MeV_1sec

eavg_nue = 14.5213    # 13.08, sqrt(<E^3>/<E>), ,<E>_nue = Hudepohl_9.7MeV_1sec
eavg_nuebar = 15.719 # 16.23, sqrt(<E^3>/<E>), ,<E>_nuebar = Hudepohl_11.7MeV_1sec

# Unmixed neutrino parameters. Luminosities and energies.
#L_nue = 10.76
#L_nuebar = 10.04 # 10^51 erg/s
#R = 1.7 # 19 km
#r_in= R * 1.0e+6 

#e_nue = 17.464    # 13.08, sqrt(<E^3>/<E>), ,<E>_nue = Hudepohl_9.7MeV_1sec
#e_nuebar = 21.8771 # 16.23, sqrt(<E^3>/<E>), ,<E>_nuebar = Hudepohl_11.7MeV_1sec

#eps_nue = 15.844    # 13.08, sqrt(<E^3>/<E>), ,<E>_nue = Hudepohl_9.7MeV_1sec
#eps_nuebar = 19.8813 # 16.23, sqrt(<E^3>/<E>), ,<E>_nuebar = Hudepohl_11.7MeV_1sec

#eavg_nue = 12.48    # 13.08, sqrt(<E^3>/<E>), ,<E>_nue = Hudepohl_9.7MeV_1sec
#eavg_nuebar = 15.74 # 16.23, sqrt(<E^3>/<E>), ,<E>_nuebar = Hudepohl_11.7MeV_1sec



# Initial entropy 
S_in= 6.0

# masses
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

def initial_T(T, L_nue=L_nue, e_nue=e_nue, L_nuebar=L_nuebar, e_nuebar=e_nuebar):

    Ye = initial_Ye(T)

    # set R_nu ~ R_gain
    r_in = (R*1.0e+6)
    r_nu_nat=r_in*(5.0e+10)
    r_in_nat=r_nu_nat

    GR_factor_angle = ((1.0 - (2.0 * GM / (r_nu_nat))) / (1.0 - (2.0 * GM / r_in_nat)))**0.5
    factor = 1.0 - (1.0 - ((r_nu_nat / r_in_nat)**2.0 / GR_factor_angle**2.0))**0.5

    rho = (T**3.0 / S_in) * 5.21 * 1e+8 
    rho_MeV = rho * 5.6e+26 * (1.97e-11)**3.0
    eta = 3.0 * 0.5 * (rho_MeV / m_n) / T**3.0
        
    qdot_cool_ann = (1.57e-25 * T**9.0 * 
                    (g1(m_e / T, eta, 0, 0, T) * g2(m_e / T, - eta, 0, 0, T) + g1(m_e / T, -eta, 0, 0, T) * g2(m_e / T,eta, 0, 0, T)) / rho_MeV)

    return ((6.8e-24 * ((Ye * L_nuebar * e_nuebar**2.0 + (1 - Ye) * L_nue * e_nue**2.0) * GR_factor_1**6.0 * (1.0e+6/r_in)**2) * factor * GR_factor_angle**6.0
            - 1.6e-24 * T**6.0) - qdot_cool_ann)


#initial_vs = [(1.88e+6 + 2 * 1e-4 * 1e+6 * i) for i in range(10)] 

#initial_vs = [ (11.9e+6 + 1e-2 * 1e+6 * i) for i in range(20) ]

# Range for the mixed 3D Bollig with Ye evol
vmin = 12.0e+6
vmax = 13.0e+6

# Range for the unmixed 3D Bollig with Ye evol
#vmin = 8.5e+6
#vmax = 10.5e+6

v_in = vmin

mu = 0.511

def fermion_func(x):
    if x == 0:
        return 0
    try:
        # Prevent overflow in the exponential function
        exp_limit = 700  # Value close to the limit where exp() overflows
        exp_arg = -mu / x
        if exp_arg < -exp_limit:
            polylog2 = 0.0
            polylog3 = 0.0
            polylog4 = 0.0
        else:
            polylog2 = float(gp.fp.polylog(2, -np.exp(exp_arg)))
            polylog3 = float(gp.fp.polylog(3, -np.exp(exp_arg)))
            polylog4 = float(gp.fp.polylog(4, -np.exp(exp_arg)))

        return - (x * (2.0 * x * mu**2.0 * polylog2 + 6.0 * x**2.0 * mu * polylog3 + 6.0 * x**3.0 * polylog4) * 4.0 * 15.0 / (np.pi**4.0 * x**4.0))
    except OverflowError:
        return float('inf')

print('gstar at 3 GK', 2 + fermion_func(0.0617))
max_mach = []

while True:

    factor = 1
    f = open('near_critical_solution.txt' ,'w')
    # initial time
    t=0.0
    dt=0.0
    dv_temp = 0.0

    # to calculate accumulated mass
    M = 0.0
    dM = 0.0

    # initial radius
    r_in_nat = r_in * 5.0e+10
    r = r_in_nat

    v_list = [v_in]
    r_list = [r]
    q_list = [0.0]

    # GR factor 
    GR_factor_1 = np.sqrt(1.0 / (1.0 - (2.0 * GM / r)))
    print('GR_factor', GR_factor_1)

    def lam_nue_n(x, GR_factor_1=GR_factor_1):
        alpha = 1.26
        GF = 1.166e-11 # MeV^-2
        delta = 1.293
        sec_to_MeVinv = 6.58e-22
        erg_toMeV = 6.24e+2 * 1e+3

        L_eff = L_nue * GR_factor_1**4 * 1e+51 * erg_toMeV * sec_to_MeVinv
        e_eff = eps_nue * GR_factor_1
        eavg_eff = eavg_nue * GR_factor_1

        return ((1 + 3*alpha**2) / (2 * np.pi**2)) * GF**2 * (L_eff/r_in_nat**2) * (e_eff + 2*delta + delta**2/eavg_eff) * (1 - x) 

    def lam_nuebar_p(x, GR_factor_1=GR_factor_1):
        alpha = 1.26
        GF = 1.166e-11 # MeV^-2
        delta = 1.293
        sec_to_MeVinv = 6.58e-22
        erg_toMeV = 6.24e+2 * 1e+3

        L_eff = L_nuebar * GR_factor_1**4 * 1e+51 * erg_toMeV * sec_to_MeVinv
        e_eff = eps_nuebar * GR_factor_1
        eavg_eff = eavg_nuebar * GR_factor_1

        return ((1 + 3*alpha**2) / (2 * np.pi**2)) * GF**2 * (L_eff/r_in_nat**2) * (e_eff - 2*delta + delta**2/eavg_eff) * (1 - x) 

    def other_lam(T):
        sec_to_MeVinv = 6.58e-22
        return 0.448 * T**5 * sec_to_MeVinv

    def initial_Ye(T):
        x = np.sqrt(1 - (r_in_nat**2/r_in_nat**2))
        lam1 = lam_nue_n(x) + other_lam(T)
        lam2 = lam1 + lam_nuebar_p(x) + other_lam(T)
        #print(lam1)
        #print(lam2)
        return lam1/lam2

    # Initial temperature
    T_in = sc.optimize.fsolve(initial_T, 4.0)[0]
    print('T_in',T_in)

    Ye = initial_Ye(T=T_in)
    Ye_list = [Ye]


    print('Initial Ye', Ye)

    printed = False

    # calculate the initial derivate of g_star
    fermion = fermion_func(T_in)
    A = 2.0 + fermion
    diff_T = 0.0001
    fermion_dT = fermion_func(T_in + diff_T)
    fermion_derivative = (fermion_dT - fermion) / (diff_T * A)
    derivative_A_term = fermion_derivative 

    # mass outflow rate
    rho_start = 5.21 * T_in**3.0 * 1e+8 / S_in # g/cm^3
    M_dot = 4.0 * np.pi * r_in**2.0 * rho_start * v_in / (2.0e+33 * GR_factor_1 ) # in solar masses
    print ('Starting density', rho_start)
    print('Mass outflow rate' , M_dot)

    print ('Initial speed (cm/s)', v_in)
    v_nat_in = v_in/(3.0e+10)
    v = v_nat_in
    S = S_in
    
    # Initial qdot set to zero at the gain radius.
    qdot = 0.0

    T = T_in

    dr = 1.0e+3 # cm 
    dr_nat = dr * 5.0e+10

    vs = ((T_in * S_in / (4.0 * m_n)) * (4.0 + (T_in * derivative_A_term)) / (3.0 + (T_in * derivative_A_term)))**0.5
    vs_list = [vs * 3e+10]
    print ('Initial sound speed (cm/s)', vs * 3e+10)

    mach = [v_nat_in / vs]

    # loop over for radii
    for p in range (50000):

        fermion = fermion_func(T)
        A = 2.0 + fermion

        GR_factor = 1.0 / (1.0 - (2.0 * GM / r))
        GR_factor_angle = ((1.0 - (2.0 * GM / r_in_nat)) / (1.0 - (2.0 * GM / r)))**0.5
        y_fac = ((1.0 - (2.0 * GM / r)) / (1.0 - v**2))**0.5
       
        density8 = (T**3.0 / S) * A / 1.055
        X_n = 828 * T**(9.0 / 8) * np.exp(- 7.074 / T) / density8**(3.0 / 4.0)

        if X_n >= 1.0:
            fac = 1.0
        else:
            fac = X_n
        
        rtxt = str(r/5.0e+10)

        if(r / 5.0e+10 < R_t):

            # Find dt corresponding to a given dr 
            dt = (dr_nat / 5.0e+10) / (v * 3.0e+10)
            t = t + dt

            if (t <= 0.6 and T<= 2.0):
                f.write(f'{t} {r / 5.0e+10} {T} {(T**3.0 / S) * A * 1e+8 / 1.055}\n')

            if T : 
                fermion = fermion_func(T)
                A = 2.0 + fermion
                diff_T = 0.0001
                fermion_dT = fermion_func(T + diff_T)
                fermion_derivative = (fermion_dT - fermion) / ( diff_T * A)
                derivative_A_term = fermion_derivative 

                beta = (1.0 / 4.0) * (1.0 + (3.0 + T * derivative_A_term)**-1)

                # Mach number
                m = v / vs
                x = 1 - factor
                lam1 = lam_nue_n(x) + other_lam(T)
                lam2 = lam1 + lam_nuebar_p(x) + other_lam(T)

                #print('lam_nu_e', lam_nue_n(x))
                #print('other lam', other_lam(T))
                #print('lam_nubar', lam_nuebar_p(x))

                dv = ((2 * vs**2 / r) - ((GM / r**2) * (1 - vs**2) * GR_factor) - (qdot * beta / (v * y_fac * (1 + 3 * vs**2)))) * dr_nat * (1.0 - v**2.0) / (v - (vs * vs / v))
                dvs = ((qdot * dr_nat / (v * y_fac * (1 + 3 * vs**2))) - (v * dv / (1 - v**2)) - (GM * GR_factor * dr_nat / r**2)) * (1 + 3 * vs**2) / (6 * vs)
                dS = (qdot * m_n / (T * v * y_fac)) * dr_nat
                dYe = (lam1 - Ye * lam2) * dr_nat / v
                #print('Ye', Ye)
                #break

                #print('dvs', dvs * 3e+10)
                S = S + dS
                r = r + dr_nat
                v = v + dv
                vs = vs + dvs
                Ye = Ye + dYe

                a = derivative_A_term 
                b = 4 - ((4 * m_n * vs**2 / S) * derivative_A_term)
                c = -12 * m_n * vs**2 / S

                T = (-b + (b**2 - 4*a*c)**0.5) / (2 * a)

                rho = (T**3.0 / S) * A * 1e+8 / 1.055
                dM = 4 * pi * rho * (r / 5.0e+10)**2.0 * (dr_nat / 5.0e+10)
                M = M + dM
                mach.append(m)

        if(T < 0):
            print ('Negative T alert!!', T)
            break

        if(S < 0):
            print ("Negative entropy alert!!", S)

        if(np.isnan(T)):
            print('T is nan !!')

        vs_list.append(vs * 3e+10)
        v_list.append(v * 3e+10)
        r_list.append(r)
        Ye_list.append(Ye)


        rho_8 = rho / 1e+8
        rho_MeV = rho * 5.6e+26 * (1.97e-11)**3.0
        eta = 3.0 * 0.5 * (rho_MeV / m_n) / T**3.0
        
        qdot_cool_ann = (1.57e-25 * T**9.0 * 
                        (g1(m_e / T, eta, 0, 0, T) * g2(m_e / T, - eta, 0, 0, T) + g1(m_e / T, -eta, 0, 0, T) * g2(m_e / T,eta, 0, 0, T)) / rho_MeV)

        factor = 1.0 - (1.0 - ((r_in_nat / r)**2.0 / GR_factor_angle**2.0))**0.5 # gravitational bending
        qdot = ((6.8e-24 * ((Ye*L_nuebar * e_nuebar**2.0 + (1 - Ye)*L_nue * e_nue**2.0) * GR_factor_1**6.0 * (1.0e+6 / r_in)**2 ) * factor * GR_factor_angle**6.0
            - 1.6e-24 * T**6.0) * fac - qdot_cool_ann)
        q_list.append(qdot)

    print('maximum mach number', max(mach))
    max_mach = mach
    #break

    if max(mach) > 1: 
        vmax = v_in
        v_in = (vmin + vmax) / 2
    elif max(mach) < 1 and max(mach) > 0.99: 
        break 
    else: 
        vmin = v_in
        v_in = (vmin + vmax) / 2

plt.semilogx(np.array(r_list)/5e+10, Ye_list)
plt.title('Ye')
plt.show()

plt.semilogx(np.array(r_list)/5e+10, max_mach)
plt.ylim([0.01, 1])
plt.title('Mach number')
plt.show()

plt.semilogx(np.array(r_list)/5e+10, np.array(v_list) )
plt.semilogx(np.array(r_list)/5e+10, np.array(vs_list) )
plt.title('Near critical plot')
plt.show()

plt.semilogx(np.array(r_list)/5e+10, np.gradient(np.array(v_list)))
plt.title('First derivative of speed')
plt.show()