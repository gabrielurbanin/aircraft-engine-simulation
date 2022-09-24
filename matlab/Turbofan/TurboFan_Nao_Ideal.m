% Inputs do GE CF34-10
M0 = 0.78;
H = 38000/3.281; %m %1 m = 3.28084 ft
[T0, a0, P0] = atmosisa(H);
% T0 = 288.15; %K
gamma_c = 1.4;
cp_c = 1004; %J/(kg.K)
gamma_t = 1.3;
cp_t = 1100; %J/(kg.K)
hpr = 42*10^6; %J/kg
pi_d_max = 0.97;
pi_b = 0.98;
pi_n = 0.92;
pi_fn = 0.99;
e_cL = 0.95;
e_cH = 0.95;
e_f = 0.98;
e_tL = 0.95;
e_tH = 0.95;
eta_b = 0.7;
eta_mL = 0.98;
eta_mH = 0.98;
P0_P9 = 0.46;
P0_P19 = 1;
Tt4 = 1850; %K
tau_n = 1;
tau_fn = 1;
pi_cL = 2.02;
pi_cH = 8.29;
pi_f = 1.7;
alfa = 5;
% Inputs extra do GE CF34-10
% P0 = 101325; %Pa
m0 = 242; %m/s

% Equations
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
tau_d = pi_d^((gamma_c - 1)/gamma_c);
tau_f = pi_f^((gamma_c - 1)/(gamma_c*e_f));
eta_f = (pi_f^((gamma_c - 1)/gamma_c) - 1)/(tau_f - 1);
tau_lambda = cp_t*Tt4/(cp_c*T0);
tau_cL = pi_cL^((gamma_c - 1)/(gamma_c*e_cL));
tau_cH = pi_cH^((gamma_c - 1)/(gamma_c*e_cH));
eta_cL = (pi_cL^((gamma_c - 1)/gamma_c) - 1)/(tau_cL - 1);
eta_cH = (pi_cH^((gamma_c - 1)/gamma_c) - 1)/(tau_cH - 1);
f = (tau_lambda - tau_r*tau_d*tau_f*tau_cL*tau_cH)/(hpr*eta_b/(cp_c*T0) - tau_lambda); %kgFuel/kgAir
tau_tH = 1 - (tau_cH - 1)/(1 + f)/tau_lambda*tau_r*tau_d*tau_f*tau_cL*eta_mH;
tau_tL = 1 - ((alfa*(tau_f - 1) + (tau_cL - 1))*eta_mL/(1 + f)*tau_r*tau_d/tau_lambda/tau_tH);
pi_tH = tau_tH^(gamma_t/((gamma_t - 1)*e_tH));
pi_tL = tau_tL^(gamma_t/((gamma_t - 1)*e_tL));
eta_tH = (tau_tH - 1)/(pi_tH^((gamma_t - 1)/gamma_t)-1);
eta_tL = (tau_tL - 1)/(pi_tL^((gamma_t - 1)/gamma_t)-1);
Pt9_P9 = P0_P9*pi_r*pi_d*pi_f*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_n;
M9 = (2/(gamma_t - 1)*(Pt9_P9^((gamma_t - 1)/gamma_t) - 1))^(1/2);
Tt9_T0 = cp_c/cp_t*tau_lambda*tau_tH*tau_tL*tau_n;
T9_T0 = Tt9_T0/Pt9_P9^((gamma_t - 1)/gamma_t);
V9_a0 = M9*(T9_T0)^(1/2);
Pt19_P19 = P0_P19*pi_r*pi_d*pi_f*pi_fn;
M19 = (2/(gamma_c - 1)*(Pt19_P19^((gamma_c - 1)/gamma_c) - 1))^(1/2);
Tt19_T0 = tau_r*tau_d*tau_f*tau_fn;
T19_T0 = Tt19_T0/Pt19_P19^((gamma_c - 1)/gamma_c);
V19_a0 = M19*(T19_T0)^(1/2);
FF_m0 = alfa/(1 + alfa)*a0*(V19_a0 - M0 + 0*T19_T0/V19_a0*(1 - P0_P19)/gamma_c); %N/(kg/s)
FC_m0 = 1/(1 + alfa)*a0*((1 + f)*V9_a0 - M0 + 0*(1 + f)*R_t/R_c*T9_T0/V9_a0*(1 - P0_P9)/gamma_c); %N/(kg/s)
F_m0 = FF_m0  + FC_m0; %N/(kg/s)
S = f/((1 + alfa)*F_m0); %(kgFuel/s)/N
FR = FF_m0/FC_m0;
eta_T = a0^2*((1 + f)*V9_a0^2 + alfa*(V19_a0^2)- (1 + alfa)*M0^2)/(2*f*hpr);
eta_P = 2*M0*((1 + f)*V9_a0 + alfa*V19_a0 - (1 + alfa)*M0)/((1 + f)*(V9_a0^2) + alfa*V19_a0^2 - (1 + alfa)*M0^2);
eta_Total = eta_P*eta_T;
F = F_m0*m0; %N
mf = S*F; %kgFuel/s vazao de combustivel
AF = 1/f; %kgAir/kgFuel
mC = m0*1/(1 + alfa); %kgAir/s vazao de ar pelo Core
mF = m0*alfa/(1 + alfa); %kgAir/s vazao de ar pelo Fan




