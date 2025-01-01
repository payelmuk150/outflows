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
pns_mass = 1.4

rho_target = (M_target / 1.98e+33 - pns_mass) * 1.98e+33 / (4 * np.pi * r_target**3 / 3)
print('Boundary rho', rho_target)

# Boundary pressure
P_target = (6 * rho_target)**(4.0/3) # P propto (S rho)^4/3, ignorong the constant factors
print('Boundary pressure', P_target)

GM=7.56395e+15 * pns_mass

# Termination Shock Radius. Set to a large value for the first part

R_t = 0.6e+13

# Unmixed neutrino parameters. Luminosities and energies.
L_nue = 8.00
L_nuebar = 8.00 # 10^51 erg/s
R = 2.0 # 19 km
r_in= R * 1.0e+6 

e_nue = 20.0    # 13.08, sqrt(<E^3>/<E>), ,<E>_nue = Hudepohl_9.7MeV_1sec
e_nuebar = 20.0 # 16.23, sqrt(<E^3>/<E>), ,<E>_nuebar = Hudepohl_11.7MeV_1sec


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

def g020(x_th, y, z):
    x = np.arange(x_th, x_th + 20, 0.1)
    integrand_020 = x * (x + z)**2.0 * np.sqrt(np.abs(x**2.0 - x_th**2.0)) / (np.exp(x + y) + 1)
    g020 = np.trapz(integrand_020, x)
    return g020

print('g020', g020(0, 0, 0))


def initial_T(T, Ye, eta, L_nue=L_nue, e_nue=e_nue, L_nuebar=L_nuebar, e_nuebar=e_nuebar):
    fermion_ = fermion_func(T)
    A = 2.0 + fermion_

    rho = (T**3.0 / S_in) * A * 1e+8 / 1.055
    rho_MeV = rho * 5.6095e+26 * (1.98e-11)**3.0

    # set R_nu ~ R_gain
    r_in = (R*1.0e+6)
    r_nu_nat=r_in*(5.0e+10)
    r_in_nat=r_nu_nat
        
    qdot_cool_ann = (1.57e-25 * T**9.0 * 
                    (g1(m_e / T, eta, 0, 0, T) * g2(m_e / T, - eta, 0, 0, T) + g1(m_e / T, -eta, 0, 0, T) * g2(m_e / T,eta, 0, 0, T)) / rho_MeV)

    return ((6.8e-24 * ((Ye * L_nuebar * e_nuebar**2.0 + (1 - Ye) * L_nue * e_nue**2.0) * (1.0e+6/r_in)**2)
            - 1.6e-24 * T**6.0) - qdot_cool_ann) / 1e-20

#vmin = 11.5e+6
#vmax = 12.5e+6

#vmin = 3.5473e+6
#vmin = 3.550091e+6
#vmin = 7.034545e+6
#vmin = 7.142578e+6

#vmin = 5.391210e+6
#vmax = 5.4e+6

#vmin = 0.445003e+6
#vmin = 6.524609e+6
#vmin = 6.524536e+6
#vmin = 2.5e+6
#vmax = 6.8e+6

#vmin = 4.974121e+6
#vmin = 4.9e+6
#vmax = 5.0e+6

vmin = 5.2e+6
#vmin = 5.407031e+6 
vmax = 5.5e+6

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
    #return 3.5


def solve_for_T(T, rho, S): 
    fermion = fermion_func(T)
    A = 2.0 + fermion
    return T - (rho * 1.055 * S / (A * 1e+8))**0.333333

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

    def initial_eta(T, Ye):
        fermion_ = fermion_func(T)
        A = 2.0 + fermion_

        rho = (T**3.0 / S_in) * A * 1e+8 / 1.055
        rho_MeV = rho * 5.6095e+26 * (1.98e-11)**3.0
        #rho_MeV = rho * 5.6e+26 * (1.97e-11)**3.0

        return 3.0 * Ye * (rho_MeV / m_n) / T**3.0


    def initials(vars):
        T = vars
        eta = initial_eta(T, 0.5)

        eq1 = initial_T(T, 0.5, eta=eta)

        return eq1

    # Initial temperature
    T_in = sc.optimize.fsolve(initials, 3.9, xtol=1e-30)[0]
    Ye = 0.5
    print('T_in',T_in)
    Ye_list = [Ye]
    print('Initial Ye', Ye)
    print('Initial eta', initial_eta(T_in, Ye))

    qdot = 0.0
    print('initial T error calculation', initial_T(T_in, Ye, initial_eta(T_in, Ye))*1e-20)
    q_list = [qdot]

    printed = False

    # calculate the initial derivate of g_star
    fermion = fermion_func(T_in)
    A = 2.0 + fermion
    diff_T = 0.0001
    fermion_dT = fermion_func(T_in + diff_T)
    fermion_derivative = (fermion_dT - fermion) / (diff_T * A)
    derivative_A_term = fermion_derivative 

    # mass outflow rate
    rho_start = A * T_in**3.0 * 1e+8 / (S_in * 1.055) # g/cm^3
    M_dot = 4.0 * np.pi * r_in**2.0 * rho_start * v_in / (2.0e+33) # in solar masses
    M_dot_list = [M_dot]
    print ('Starting density', rho_start)
    print('Mass outflow rate' , M_dot)

    print ('Initial speed (cm/s)', v_in)
    v_nat_in = v_in/(3.0e+10)
    v = v_nat_in
    S = S_in

    T = T_in

    dr = 1.0e+3 # cm 
    dr_nat = dr * 5.0e+10

    vs = ((T_in * S_in / (4.0 * m_n)) * (4.0 + (T_in * derivative_A_term)) / (3.0 + (T_in * derivative_A_term)))**0.5
    vs_list = [vs * 3e+10]
    S_list = [S]
    print ('Initial sound speed (cm/s)', vs * 3e+10)

    mach = [v_nat_in / vs]
    rho_list = [rho_start]
    T_list = [T_in]

    fermion = fermion_func(T)
    A = 2.0 + fermion

    density8 = (T**3.0 / S) * A / 1.055
    # gm_to_MeV = 5.6095e+26, cm^-1 to MeV = 1.98e-11. Taken from Saha conversion table: https://www.saha.ac.in/theory/palashbaran.pal/conv.html
    rho_MeV = density8 * 1e+8 * 5.6095e+26 * (1.98e-11)**3.0
    rho = rho_MeV / (5.6095e+26 * (1.98e-11)**3.0)

    # loop over for radii
    for p in range (60000):
        #print(density8 * 1e+8)
        #print(rho_MeV)
        #break
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
                #T = vs**2 * m_n / (beta * S)
                fermion = fermion_func(T)
                A = 2.0 + fermion
                diff_T = 0.0001
                fermion_dT = fermion_func(T + diff_T)
                fermion_derivative = (fermion_dT - fermion) / ( diff_T * A)
                derivative_A_term = fermion_derivative 

                beta = (1.0 / 4.0) * (1.0 + (3.0 + T * derivative_A_term)**-1)
                vs = ((T * beta * S)/m_n)**0.5

                # Mach number
                m = v / vs
                x = 1 - factor

                dv = ((2 * vs**2 / r) - (GM / r**2) - (qdot * beta / v)) * dr_nat / (v - (vs * vs / v))
                dS = (qdot * m_n / (T * v)) * dr_nat
                drho_MeV = - ((2 * r * rho_MeV * v * dr_nat) + (r**2 * rho_MeV * dv)) / (r**2 * v)

                S = S + dS
                r = r + dr_nat
                v = v + dv
                rho_MeV = rho_MeV + drho_MeV 

                rho = rho_MeV / (5.6095e+26 * (1.98e-11)**3.0)
                density8 = rho / 1e+8

                T_guess = (rho * 1.055 * S / (A * 1e+8))**0.333333
                T = sc.optimize.fsolve(solve_for_T, T_guess, args=(rho, S))[0]
                
                dM = 4 * pi * rho * (r / 5.0e+10)**2.0 * (dr_nat / 5.0e+10)
                M = M + dM
                M_dot_ = 4.0 * np.pi * (r / 5e+10)**2.0 * (rho / 2e+33) * (v * 3e+10) # in solar masses

        if(T < 0):
            print ('Negative T alert!!', T)
            break

        if(S < 0):
            print ("Negative entropy alert!!", S)

        if(np.isnan(T)):
            print('T is nan !!')
            break

        rho8 = rho/1e+8
        vs_list.append(vs * 3e+10)
        v_list.append(v * 3e+10)
        r_list.append(r)
        Ye_list.append(Ye)
        mach.append(m)
        M_dot_list.append(M_dot_)
        rho_list.append(rho)
        S_list.append(S)
        T_list.append(T)

        eta = 3.0 * Ye * (rho_MeV / m_n) / T**3.0
        
        qdot_cool_ann = (1.57e-25 * T**9.0 * 
                        (g1(m_e / T, eta, 0, 0, T) * g2(m_e / T, - eta, 0, 0, T) + g1(m_e / T, -eta, 0, 0, T) * g2(m_e / T,eta, 0, 0, T)) / rho_MeV)

        factor = 1.0 - (1.0 - ((r_in_nat / r)**2.0))**0.5 # gravitational bending
        qdot = ((6.8e-24 * ((Ye*L_nuebar * e_nuebar**2.0 + (1 - Ye)*L_nue * e_nue**2.0) * (1.0e+6 / r_in)**2 ) * factor
            - 1.6e-24 * T**6.0) * fac - qdot_cool_ann)
        q_list.append(qdot)

    print('maximum mach number', max(mach))
    max_mach = mach
    #break

    if max(mach) > 1: 
        vmax = v_in
        v_in = (vmin + vmax) / 2
        print('new vin', v_in)
    elif max(mach) < 1 and max(mach) > 0.97: 
        break 
    else: 
        #print('wtf')
        vmin = v_in
        v_in = (vmin + vmax) / 2
plt.semilogx(np.array(r_list)/5e+10, q_list)
plt.title('Mdot')
plt.show()

plt.semilogx(np.array(r_list)/5e+10, M_dot_list)
plt.title('Mdot')
plt.show()

plt.loglog(np.array(r_list)/5e+10, np.array(rho_list))
plt.title('Density')
plt.show()

plt.loglog(np.array(r_list)/5e+10, np.array(T_list))
plt.title('Temperature')
plt.show()

plt.semilogx(np.array(r_list)/5e+10, S_list)
plt.title('S')
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