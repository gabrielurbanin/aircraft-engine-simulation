% Versões
% Tipo = 1: tau_t conhecido
% Tipo = 2:potência (W) conhecida
tipo = 2;

% Inputs do PT6A A8
M0 = 0.03;
T0 = 288.15; %K
gamma_c = 1.4;
cp_c = 1004; %J/(kg.K)
gamma_t = 1.3;
cp_t = 1108; %J/(kg.K)
hpr = 42*10^6; %J/kg
pi_d_max = 0.99;
pi_b = 1;
pi_n = 1;
e_c = 0.95;
e_tL = 0.98;
e_tH = 0.98;
eta_b = 0.98;
eta_g = 0.9409;
eta_mL = 0.98;
eta_mH = 0.98;
eta_prop = 0.24;
Tt4 = 1350; %K
pi_c = 8.7;
if tipo == 1
    tau_t = 0.646;
end
if tipo == 2
    W_prop = 1164*10^3; %W
end
% Inputs extra do PT6A A8
P0 = 101325; %Pa
m0 = 4.78; %m/s

%Equations
R_c = (gamma_c - 1)/gamma_c*cp_c; %J/(kg.K)
R_t = (gamma_t - 1)/gamma_t*cp_t; %J/(kg.K)
a0 = (gamma_c*R_c*T0)^(1/2); %m/s
V0 = a0*M0; %m/s
tau_r = 1 + (gamma_c - 1)/2*M0^2;
pi_r = tau_r^(gamma_c/(gamma_c - 1));
if M0 <= 1
    eta_r = 1;
else
    eta_r = 1 - 0.075*(M0 - 1)^1.35;
end
pi_d = pi_d_max*eta_r;
tau_lambda = cp_t*Tt4/(cp_c*T0);
tau_c = pi_c^((gamma_c - 1)/(gamma_c*e_c));
eta_c = (pi_c^((gamma_c - 1)/gamma_c) - 1)/(tau_c - 1);
f = (tau_lambda - tau_r*tau_c)/(hpr*eta_b/(cp_c*T0) - tau_lambda); %kgFuel/kgAir
tau_tH = 1 - tau_r*(tau_c - 1)/(eta_mH*(1 + f)*tau_lambda);
pi_tH = tau_tH^(gamma_t/((gamma_t - 1)*e_tH));
eta_tH = (1 - tau_tH)/(1 - tau_tH^(1/e_tH));
if tipo == 1
    tau_tL = tau_t/tau_tH;
    C_prop = eta_prop*eta_g*eta_mL*(1 + f)*tau_lambda*tau_tH*(1 - tau_tL);
end
if tipo == 2
    C_prop = eta_prop*W_prop/(m0*cp_c*T0);
    tau_tL = 1 - C_prop/(eta_prop*eta_g*eta_mL*(1 + f)*tau_lambda*tau_tH);
    tau_t = tau_tH*tau_tL;
end
pi_tL = tau_tL^(gamma_t/((gamma_t - 1)*e_tL));
eta_tL = (1 - tau_tL)/(1 - tau_tL^(1/e_tL));
Pt9_P0 = pi_r*pi_d*pi_c*pi_b*pi_tH*pi_tL*pi_n;
if Pt9_P0 > ((gamma_t + 1)/2)^(gamma_t/(gamma_t - 1))
    M9 = 1;
    Pt9_P9 = ((gamma_t + 1)/2)^(gamma_t/(gamma_t - 1));
    P0_P9 = Pt9_P9/Pt9_P0;
else
    P0_P9 = 1;
    Pt9_P9 = Pt9_P0;
    M9 = (2/(gamma_t - 1)*(Pt9_P0^((gamma_t - 1)/gamma_t) - 1))^(1/2);
end
V9_a0 = sqrt(2*tau_lambda*tau_tH*tau_tL/(gamma_c - 1)*(1 - (Pt9_P9)^(-1*(gamma_t - 1)/gamma_t)));
Tt9_T0 = tau_lambda*tau_tH*tau_tL;
T9_T0 = Tt9_T0/(Pt9_P9^((gamma_t - 1)/gamma_t));
C_c = (gamma_c - 1)*M0*((1 + f)*V9_a0 - M0 + (1 + f)*R_t/R_c*T9_T0/V9_a0*(1 - P0_P9)/gamma_c);
C_Total = C_prop + C_c;
F_m0 = C_Total*cp_c*T0/V0;
S = f/F_m0; %kg/s/N
S_P = f/(C_Total*cp_c*T0); %kgFuel/s/W vazão mássica por potência
F = m0*F_m0; %N
F_core = C_c*m0*cp_c*T0/V0;
F_prop = C_prop*m0*cp_c*T0/V0;
W_core = F_core*V0;
if tipo == 1
    W_prop = F_prop*V0;
    W_m0 = C_Total*cp_c*T0; %W/(kg/s)
    W = m0*W_m0; %W
end
if tipo == 2
    W = W_prop + W_core;
    W_m0 = W/m0;
end
eta_P = C_Total/(C_prop/eta_prop + ((gamma_c - 1)/2)*((1 + f)*V9_a0^2 - M0^2));
eta_T = C_Total*cp_c*T0/(f*hpr);
eta_Total = eta_P*eta_T;
AF = 1/f;
mf = S*F; %kgFuel/s


