tipo = 2;

if tipo == 1
    % Inputs do Rolls-Royce Nene
    M0 = 0.03;
    T0 = 288.15; %K
    gamma = 1.4;
    cp = 1005; %J/(kg.K)
    hpr = 41*10^6; %J/kg
    Tt4 = 1073.15; %K
    pi_c = 65/15;
    % Inputs extras do Rolls-Royce Nene
    P0 = 101325; %Pa
end

if tipo == 2
    % Inputs do VT80
    M0 = 0.03;
    T0 = 288.15; %K
    gamma = 1.4;
    cp = 1004; %J/(kg.K)
    hpr = 42.8*10^6; %J/kg
    Tt4 = 1100; %K
    pi_c = 3;
    % Inputs extras do VT80
    P0 = 101325; %Pa
end

% Equations
R = (gamma - 1)/gamma*cp; %J/(kg.K)
a0 = (gamma*R*T0)^(1/2); %m/s
tau_r = 1 + (gamma - 1)/2*M0^2;
tau_lambda = Tt4/T0;
tau_c = pi_c^((gamma - 1)/gamma);
tau_t = 1 - tau_r/tau_lambda*(tau_c - 1);
V9_a0 = (2/(gamma - 1)*tau_lambda/(tau_r*tau_c)*(tau_r*tau_c*tau_t - 1))^(1/2);

% Outputs
F_m0 = a0*(V9_a0 - M0); %N/(kg/s)
f = cp*T0/hpr*(tau_lambda - tau_r*tau_c); %kgFuel/kgAir
S = f/F_m0; %(kgFuel/s)/N
%S = S*1000*3600;  %(kgFuel/h)/kN
eta_T = 1 - 1/(tau_r*tau_c);
eta_P = 2*M0/(V9_a0 + M0);
eta_Total = eta_P*eta_T;

% Outputs extras
V9 = V9_a0*a0; %m/s
M9 = V9/(gamma*R*(630+273.15))^(1/2);  %630+273.15 K é a temperatura estática na seçao de saída (?)
AF = 1/f;  %kgAir/kgFuel
if tipo == 1
    m0 = 22240/V9; %kg/s (?) Rolls-Royce Nene
end
if tipo == 2
    m0 = 80/sqrt(gamma*R*(650+273.15)); %kg/s (?) VT80
end
F = F_m0*m0;  %N
FC = S*F*60; %(kgFuel/min) seria mf (?)
MFP_c = sqrt(gamma/R)*0.7*(1+(gamma-1)/2*0.7^2)^((gamma+1)/2/(gamma-1)); %considerando M_c = 0.7  (?)
T5 = Tt4*tau_t; %K
pi_t = tau_t^(gamma/(gamma-1));
P5 = P0*pi_c*pi_t; %K
A8 = m0*(T5)^(0.5)/P5/MFP_c; %m2
D8 = ((A8/pi)*4)^(0.5); %m

% Todos os Outputs
F_m0 %N/(kg/s)
f %kgFuel/kgAir
S %(kgFuel/s)/N
%S*1000*3600;  %(kgFuel/h)/kN
eta_T
eta_P
eta_Total
V9 %m/s
M9  %630+273.15 K é a temperatura estática na seçao de saída (?)
AF  %kgAir/kgFuel
m0 %kg/s
F  %N
FC %(kgFuel/min) seria mf (?)
MFP_c %considerando M_c = 0.7  (?)
T5 %K
pi_t
P5 %Pa
A8 %m2
D8 %m


