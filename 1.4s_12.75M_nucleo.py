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
Rs = 1.4e+9
z_s = 0

for r in f_target[:, 2]:
    if r >= Rs:
        r_target = r 
        M_target = f_target[z_s, 1]
        
        break
    z_s += 1

pns_mass = 1.8

rho_target = (M_target / 1.98e+33 - pns_mass) * 1.98e+33 / (4 * np.pi * r_target**3 / 3)
print('Boundary rho', rho_target)

# Boundary pressure
P_target = (6 * rho_target)**(4.0/3) # P propto (S rho)^4/3, ignorong the constant factors
print('Boundary pressure', P_target)

GM=7.56395e+15 * pns_mass
r_front = 3.0e+9

# This is the radius of the termination shock. I have set it to a huge number here to make it redundant 
# in the context of subsonic outflows
R_t = 1.8e+12

# neutrino parameters. Luminosities and energies.
L_nue = 6.12
L_nuebar = L_nue / 1.22 # 10^51 erg/s
R = 1.9 # 19 km
r_in= R * 1.0e+6 
e_nuebar = 16.23 # 16.23, sqrt(<E^3>/<E>), ,<E>_nuebar = Hudepohl_11.7MeV_1sec
e_nue = 13.08    # 13.08, sqrt(<E^3>/<E>), ,<E>_nue = Hudepohl_9.7MeV_1sec

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

def initial_T(T):
    fermion_ = fermion_func(T)
    A = 2.0 + fermion_

    rho = (T**3.0 / S_in) * A * 1e+8 / 1.055
    rho_MeV = rho * 5.6095e+26 * (1.98e-11)**3.0

    # set R_nu ~ R_gain
    r_in = (R*1.0e+6)
    r_nu_nat=r_in*(5.0e+10)
    r_in_nat=r_nu_nat

    GR_factor_angle = ((1.0 - (2.0 * GM / (r_nu_nat))) / (1.0 - (2.0 * GM / r_in_nat)))**0.5
    factor = 1.0 - (1.0 - ((r_nu_nat / r_in_nat)**2.0 / GR_factor_angle**2.0))**0.5
    eta = 3.0 * 0.5 * (rho_MeV / m_n) / T**3.0
        
    qdot_cool_ann = (1.57e-25 * T**9.0 * 
                    (g1(m_e / T, eta, 0, 0, T) * g2(m_e / T, - eta, 0, 0, T) + g1(m_e / T, -eta, 0, 0, T) * g2(m_e / T,eta, 0, 0, T)) / rho_MeV)

    return ((3.4e-24 * ((L_nuebar * e_nuebar**2.0 + L_nue * e_nue**2.0) * GR_factor_1**6.0 * (1.0e+6/r_in)**2) * factor * GR_factor_angle**6.0
            - 1.6e-24 * T**6.0) - qdot_cool_ann)

#initial_vs = [(1.88e+6 + 2 * 1e-4 * 1e+6 * i) for i in range(10)] 
initial_vs = [(2.15e+6 + 2 * 1e-4 * 1e+6 * i) for i in range(10)] 
mu = 0.511

def fermion_func(x) :
    return - (x * (2.0 * x * mu**2.0 * float(gp.fp.polylog(2, - np.exp(-mu / x))) + 6.0 * x**2.0 * mu * float(gp.fp.polylog(3, - np.exp(-mu / x))) 
            + 6.0 * x**3.0 * float(gp.fp.polylog(4, - np.exp(-mu / x)) )) * 4.0 * 15.0 / (np.pi**4.0 * x**4.0))

print('gstar at 3 GK', 2 + fermion_func(0.0617))

for v_in in initial_vs:

    f = open('Alex_checks.txt' ,'w')
    z = z_s

    # initial time
    t=0.0
    dt=0.0

    # to calculate accumulated mass
    M = 0.0
    dM = 0.0

    # initial radius
    r_in_nat = r_in * 5.0e+10
    r = r_in_nat

    v_list = [v_in]
    r_list = [r]

    # GR factor 
    GR_factor_1 = np.sqrt(1.0 / (1.0 - (2.0 * GM / r)))
    print('GR_factor', GR_factor_1)

    # Initial temperature
    T_in = sc.optimize.fsolve(initial_T, 4.0)[0]
    print('T_in',T_in)

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
    #M_dot = 4.0 * np.pi * r_in**2.0 * rho_start * v_in / (2.0e+33 * GR_factor_1 ) # in solar masses
    print ('Starting density', rho_start)

    print ('Initial speed (cm/s)', v_in)
    v_nat_in = v_in/(3.0e+10)
    v = v_nat_in
    S = S_in
    
    # Initial qdot set to zero at the gain radius.
    qdot = 0.0

    T = T_in
    T1 = T_in

    dr = 5.0e+2 # cm 
    dr_nat = dr * 5.0e+10

    vs = ((T_in * S_in / (4.0 * m_n)) * (4.0 + (T_in * derivative_A_term)) / (3.0 + (T_in * derivative_A_term)))**0.5
    print ('Initial sound speed (cm/s)', vs * 3e+10)

    A = 2.0 + fermion
    density8 = (T**3.0 / S) * A / 1.055
    rho_MeV = density8 * 1e+8 * 5.6095e+26 * (1.98e-11)**3.0
    rho = rho_MeV / (5.6095e+26 * (1.98e-11)**3.0)
    y_fac_in = ((1.0 - (2.0 * GM / r)) / (1.0 - v**2 ))**0.5

    M_dot = 4.0 * np.pi * (r / 5e+10)**2.0 * (rho / 2e+33) * (v * 3e+10) * y_fac_in
    M_dot_list = [M_dot]
    S_list = [S + 11 + np.log(T**1.5 / density8)]
    T_list = [T_in * 11.604]
    print('Mass outflow rate' , M_dot)

    # loop over for radii
    for p in range (1000000):

        fermion = fermion_func(T)
        A = 2.0 + fermion

        GR_factor = 1.0 / (1.0 - (2.0 * GM / r))
        GR_factor_angle = ((1.0 - (2.0 * GM / r_in_nat)) / (1.0 - (2.0 * GM / r)))**0.5
        y_fac = ((1.0 - (2.0 * GM / r)) / (1.0 - v**2))**0.5
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

            if(t > 0.6 and not printed):
                printed = True
                T_break = T
                rho_break = (T**3.0 / S) * A *1e+8 / 1.055
                time_break = t
                v_break = v
                r_break = r
                A_break = A
                r_hubble = r_break / 5.0e+10
                T1 = T
                rho_ex = rho_break
                T_ex = T_break

            if (t <= 0.6 and T<= 0.86):
                f.write(f'{t} {r / 5.0e+10} {T} {(T**3.0 / S) * A * 1e+8 / 1.055}\n')

            elif (t > 0.6 and T <= 0.86):

                if z < 1160: 
                    dt_break = (f_target[z+1, 2] - f_target[z, 2]) / 1e+9
                    time_break += dt_break
                    r_hubble += v_break * 3e+10 * dt_break
                    rho_factor = (f_target[z + 1, 1] - pns_mass * 1.98e+33) * f_target[z, 2]**3 / ((f_target[z, 1] - pns_mass * 1.98e+33) * f_target[z + 1, 2]**3)
                    rho_ex *= rho_factor
                    T_ex *= rho_factor**0.3333
                    z+=1
                    t_expansion = time_break
                    f.write(f'{time_break} {r_hubble} {T_ex} {rho_ex}\n')

                else: 
                    t_expansion = t_expansion * 1.0001
                    radius_expansion = r_hubble * (t_expansion / time_break)
                    rho_expansion = rho_ex * (time_break / t_expansion)**3
                    T_expansion = T_ex * (time_break / t_expansion)
                    f.write(f'{t_expansion} {radius_expansion} {T_expansion} {rho_expansion}\n')

            if T : 
                fermion = fermion_func(T)
                A = 2.0 + fermion
                diff_T = 0.0001
                fermion_dT = fermion_func(T + diff_T)
                fermion_derivative = (fermion_dT - fermion) / ( diff_T * A)
                derivative_A_term = fermion_derivative 

                beta = (1.0 / 4.0) * (1.0 + (3.0 + T * derivative_A_term)**-1)
                vs = ((T * beta * S)/m_n)**0.5

                dv = ((2 * vs**2 / r) - ((GM / r**2) * (1 - vs**2) * GR_factor) - (qdot * beta / (v * y_fac * (1 + ((T * S)/m_n) )))) * dr_nat * (1.0 - v**2.0) / (v - (vs * vs / v))
                dS = (qdot * m_n / (T * v * y_fac)) * dr_nat
                dy_fac = (1 / (2.0*y_fac)) * ((1-v**2)*(2*GM/r**2)*dr_nat + (1 - (2*GM/r)) * 2*v*dv) / (1 - v**2)**2
                drho_MeV = - ((2 * r * rho_MeV * v * dr_nat * y_fac) + (r**2 * rho_MeV * y_fac * dv) + (r**2 * rho_MeV * v * dy_fac)) / (r**2 * v * y_fac)

                S = S + dS
                r = r + dr_nat
                v = v + dv
                rho_MeV = rho_MeV + drho_MeV
                rho = rho_MeV / (5.6095e+26 * (1.98e-11)**3.0)
                density8 = rho / 1e+8
                T = (rho * 1.055 * S / (A * 1e+8))**0.333333
                dM = 4 * pi * rho * (r / 5.0e+10)**2.0 * (dr_nat / 5.0e+10)
                M = M + dM
                M_dot_ = 4.0 * np.pi * (r / 5e+10)**2.0 * (rho / 2e+33) * (v * 3e+10) * y_fac # in solar masses

                eta = 3.0 * 0.5 * (rho_MeV / m_n) / T**3.0
        
                qdot_cool_ann = (1.57e-25 * T**9.0 * 
                                (g1(m_e / T, eta, 0, 0, T) * g2(m_e / T, - eta, 0, 0, T) + g1(m_e / T, -eta, 0, 0, T) * g2(m_e / T,eta, 0, 0, T)) / rho_MeV)

                factor = 1.0 - (1.0 - ((r_in_nat / r)**2.0 / GR_factor_angle**2.0))**0.5 # gravitational bending
                qdot = ((3.4e-24 * ((L_nuebar * e_nuebar**2.0 + L_nue * e_nue**2.0) * GR_factor_1**6.0 * (1.0e+6 / r_in)**2 ) * factor * GR_factor_angle**6.0
                    - 1.6e-24 * T**6.0) * fac - qdot_cool_ann)

                M_dot_list.append(M_dot_)


        if(T < 0):
            print ('Negative T alert!!', T)
            break

        if(S < 0):
            print ("Negative entropy alert!!", S)

        if(np.isnan(T)):
            print('T is nan !!')

        v_list.append(v * 3e+10)
        r_list.append(r)
        S_list.append(S + 11 + np.log(T**1.5 / density8))
        T_list.append(T * 11.604)
        #rho_MeV = rho * 5.6e+26 * (1.97e-11)**3.0

    P_calculation = (S * (T**3.0 / S) * A * 1e+8)**(4.0/3) # P propto (S rho)^4/3, ignorong the constant factors
    convergence_fraction = abs(P_calculation - P_target) / P_target

    print('Calculated far boundary S: ', S)
    print('Calculated far boundary rho: ', (T**3.0 / S) * A * 1e+8)
    print('Calculated far pressure: ', P_calculation)
    print('Target far pressure: ', P_target)
    print('fractional convergence: ', convergence_fraction)

    f.close()

    plt.loglog(np.array(r_list)[1:]/5e+10, M_dot_list[1:])
    plt.show()

    plt.loglog(np.array(r_list)/5e+10, S_list)
    plt.show()

    plt.loglog(np.array(r_list)/5e+10, T_list)
    plt.show()

    plt.loglog(np.array(r_list)/5e+10, v_list)
    plt.show()

    if convergence_fraction < 0.1:
        break

f = np.loadtxt('Alex_checks.txt')
plt.loglog(f[:, 0], f[:, 1])
plt.show()
plt.loglog(f[:, 1], f[:, 3])
plt.show()
plt.loglog(f[:, 1], f[:, 2])
plt.show()

plt.loglog(f[:, 1], np.gradient(f[:, 1], f[:, 0]) / 1e+5, label = 'speed')
plt.loglog(f[:, 1], (10000) * f[:, 1] / 1e+9, label = 'Hubble' )
plt.show()

print('S', S)
print('rho_f',(T**3.0/S)*A*1e+8 / 1.055)

P_calculation = (S * (pow(T,3.0)/S)*A*1e+8 / 1.055)**(4.0/3) # P propto (S rho)^4/3, ignorong the constant factors
print('Calculated boundary pressure', P_calculation)




