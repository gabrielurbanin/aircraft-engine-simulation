import math

tipo = 2

if tipo == 1:
    # Inputs do Rolls-Royce Nene
    M0 = 2
    T0 = 216.7 #K
    gamma_c = 1.4
    cp_c = 1004 #J/(kg.K)
    gamma_t = 1.3
    cp_t = 1239 #J/(kg.K)
    hpr = 42.8*10**6 #J/kg
    pi_d_max = 0.95
    pi_b = 0.94
    pi_n = 0.96
    e_c = 0.9
    e_t = 0.9
    eta_b = 0.98
    eta_m = 0.99
    P0_P9 = 0.5
    Tt4 = 1800 #K
    pi_c = 10
    # Inputs extras do Rolls-Royce Nene
    P0 = 101325 #Pa
    tau_d = 1

if tipo == 2:
    # Inputs do VT80
    M0 = 0.5
    T0 = 288.15 #K
    gamma_c = 1.4
    cp_c = 1004 #J/(kg.K)
    gamma_t = 1.3
    cp_t = 1239 #J/(kg.K)
    hpr = 42.8*10**6 #J/kg
    pi_d_max = 0.95
    pi_b = 0.95
    pi_n = 0.987
    e_c = 0.7
    e_t = 0.81
    eta_b = 0.99
    eta_m = 0.98
    P0_P9 = 1
    Tt4 = 920 #K
    pi_c = 3.1
    # Inputs extras do VT80
    P0 = 101325 #Pa
    tau_d = 1

# Equations
R_c = (gamma_c - 1)/gamma_c*cp_c #J/(kg.K)
R_t = (gamma_t - 1)/gamma_t*cp_t #J/(kg.K)
a0 = (gamma_c*R_c*T0)**(1/2) #m/s
V0 = a0*M0
tau_r = 1 + (gamma_c - 1)/2*M0**2
pi_r = tau_r**(gamma_c/(gamma_c - 1))
if M0 <= 1:
    eta_r = 1
else:
    eta_r = 1 - 0.075*(M0 - 1)**1.35

pi_d = pi_d_max*eta_r
tau_lambda = cp_t*Tt4/(cp_c*T0)
tau_c = pi_c**((gamma_c - 1)/(gamma_c*e_c))
eta_c = (pi_c**((gamma_c - 1)/gamma_c) - 1)/(tau_c - 1)
f = (tau_lambda - tau_r*tau_c)/(hpr*eta_b/(cp_c*T0) - tau_lambda) #kgFuel/kgAir
tau_t = 1 - 1/(eta_m*(1 + f))*tau_r/tau_lambda*(tau_c - 1)
pi_t = tau_t**(gamma_t/((gamma_t - 1)*e_t))
eta_t = (1 - tau_t)/(1 - tau_t**(1/e_t))
Pt9_P9 = P0_P9*pi_r*pi_d*pi_c*pi_b*pi_t*pi_n
M9 = (2/(gamma_t - 1)*(Pt9_P9**((gamma_t - 1)/gamma_t) - 1))**(1/2)
T9_T0 = tau_lambda*tau_t/(Pt9_P9**((gamma_t - 1)/gamma_t))*cp_c/cp_t
V9_a0 = M9*(gamma_t*R_t/(gamma_c*R_c)*T9_T0)**(1/2)
F_m0 = a0*((1 + f)*V9_a0 - M0 + (1 + f)*R_t*T9_T0/(R_c*V9_a0)*(1 - P0_P9)/gamma_c) #N/(kg/s)
S = f/F_m0 #(kgFuel/s)/N
eta_T = a0**2*((1 + f)*V9_a0**2 - M0**2)/(2*f*hpr)
eta_P = 2*V0*F_m0/(a0**2*((1 + f)*V9_a0**2 - M0**2))
eta_Total = eta_P*eta_T

# Outputs extras
V9 = V9_a0*a0 #m/s
AF = 1/f  #kgAir/kgFuel
if tipo == 1:
    m0 = 22240/V9 #kg/s (?) Rolls-Royce Nene

if tipo == 2:
    m0 = 0.28917 #kg/s (?) VT80

F = F_m0*m0  #N
FC = S*F*60 #(kgFuel/min) seria mf (?)
MFP_c = math.sqrt(gamma_c/R_t)*0.7*(1+(gamma_c-1)/2*0.7**2)**((gamma_c+1)/2/(gamma_c-1)) #considerando M_c = 0.7  (?) e gamma_c e R_t (?)
T5 = Tt4*tau_t #K
P5 = P0*pi_c*pi_t #K
A8 = m0*(T5)**(0.5)/P5/MFP_c #m2
D8 = ((A8/pi)*4)**(0.5) #m
Pt4 = P0*pi_r*pi_d*pi_c*pi_b #Pa
Pt9 = P0*pi_r*pi_d*pi_c*pi_b*pi_t #Pa
tau_b = Tt4/(T0*tau_r*tau_d*tau_c)
Tt9 = T0*tau_r*tau_d*tau_c*tau_b*tau_t #K
T9 = T0*T9_T0 #K

# Todas as Outputs
print("F_m0:",F_m0) #N/(kg/s)
print("f:",f) #kgFuel/kgAir
print("S:",S) #(kgFuel/s)/N
#S*1000*3600  #(kgFuel/h)/kN
print("eta_T:",eta_T)
print("eta_P:",eta_P)
print("eta_Total:",eta_Total)
print("eta_c:",eta_c)
print("eta_t:",eta_t)
print("V9:",V9) #m/s
print("AF:",AF)  #kgAir/kgFuel
print("m0:",m0) #kg/s
print("F:",F)  #N
print("FC:",FC) #(kgFuel/min) seria mf (?)
print("MFP_c:",MFP_c) #considerando M_c = 0.7  (?) e gamma_c e R_t (?)
print("T5:",T5) #K
print("P5:",P5) #Pa
print("A8:",A8) #m2
print("D8:",D8) #m
print("Pt4:",Pt4) #Pa
print("Pt9:",Pt9) #Pa
print("tau_b:",tau_b)
print("Tt9:",Tt9) #K
print("T9:",T9) #K



