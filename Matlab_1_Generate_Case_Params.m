clear
clc
PScase='Case9_3IBR';

%-------------------------------------------------------------------
%-------------------------- Network base values --------------------
%-------------------------------------------------------------------

S_b=100*(10^6);                          % Base power (three phase VA rating)
V_bLL=230*(10^3);                        % L-L rms voltage
V_b=sqrt(2/3)*V_bLL;                     % Base voltage (Vph peak value)
I_b=(2/3)*S_b/V_b;                       % Base current
f_b=50;w_b=2*pi*f_b;                     % Base frequency
Z_b=V_b/I_b;L_b=Z_b/w_b;C_b=1/(w_b*Z_b); % Base imp., ind., cap.
w_0=1;

%-------------------------------------------------------------------
%------------------------- Bus parameters  -------------------------
%-------------------------------------------------------------------

file="CSVs\"+PScase+"\Bus.csv"; opts = detectImportOptions(file); opts.VariableNamingRule='preserve';
bus_pu = table2array(readtable(file,opts)); 
Nn = size(bus_pu,1); 

% nominal voltages of each bus
Bus_VbLL=bus_pu(:,3);
Bus_Vb=sqrt(2/3).*Bus_VbLL;
Bus_Ib=(2/3)*S_b./Bus_Vb;
Bus_Zb=Bus_Vb./Bus_Ib; 
Bus_Cb=1 ./(w_b*Bus_Zb);

% shunt capacitors
busC    =   bus_pu(:,2);
busC_SI =   busC.*Bus_Cb;
Nsh=sum(busC~=0);
C_sh=busC(busC~=0);

%-------------------------------------------------------------------
%-------------------------- Branch parameters   --------------------
%-------------------------------------------------------------------

file="CSVs\"+PScase+"\Line.csv"; opts = detectImportOptions(file); opts.VariableNamingRule='preserve';
br_pu =  table2array(readtable(file,opts)); 
Ne=size(br_pu,1);

line_from    =   br_pu(:,1);   
line_to      =   br_pu(:,2);
% RL lines
Re         =     br_pu(:,3);
Le         =     br_pu(:,4);
Re_SI      =     Re*Z_b;
Le_SI      =     Le*L_b;

%-------------------------------------------------------------------
%---------------------------- Load parameters   --------------------
%-------------------------------------------------------------------

file="CSVs\"+PScase+"\Loads.csv"; opts = detectImportOptions(file); opts.VariableNamingRule='preserve';
load_pu =  table2array(readtable(file,opts)); 
Nload=round(load_pu(end,1)); load_bus=(load_pu(:,2));
a_p=load_pu(:,3); b_p=load_pu(:,4); c_p=load_pu(:,5);
a_q=load_pu(:,6); b_q=load_pu(:,7); c_q=load_pu(:,8);
Pload_max=load_pu(:,9);Pload_min=load_pu(:,10);
Qload_max=load_pu(:,11);Qload_min=load_pu(:,12);
Pdif_max=load_pu(:,13); Qdif_max=load_pu(:,14);

%-------------------------------------------------------------------
%--------------------------- IBRs data -----------------------------
%-------------------------------------------------------------------

file="CSVs\"+PScase+"\IBR.csv"; opts = detectImportOptions(file); opts.VariableNamingRule='preserve';
gen_pu =  table2array(readtable(file,opts)); 
Nibr=size(gen_pu,1); ibrbus=gen_pu(:,1);

% Base values AC side with generators
V_b_busgen=Bus_Vb(ibrbus); I_b_busgen=(2/3)*S_b./V_b_busgen;
Z_b_busgen=V_b_busgen./I_b_busgen; L_b_busgen=Z_b_busgen./w_b; C_b_busgen=1./(w_b*Z_b_busgen);

R_f=gen_pu(:,5); L_f=gen_pu(:,6); C_f=gen_pu(:,7);
R_f_SI=R_f.*Z_b_busgen;
L_f_SI=L_f.*L_b_busgen;
C_f_SI=C_f.*C_b_busgen;

% Base values DC side
V_b_dc=   (2)*V_b_busgen;
I_b_dc= (3/4).*I_b_busgen;
R_b_dc= (8/3)*Z_b_busgen;
C_b_dc= (3/8)*C_b_busgen;
L_b_dc= (8/3)*L_b_busgen;

Cdc=gen_pu(:,3);
Gdc=gen_pu(:,4);
Gdc_SI=Gdc./R_b_dc;
Cdc_SI=Cdc.*C_b_dc;

% vdc nominal value
Vdc_n=gen_pu(:,2);
Vdc_n_SI=Vdc_n.*V_b_dc;

% Max. apparent power and current
s_ibr_nom=gen_pu(:,8);

% Max. IBR current
ix_ibr_max=gen_pu(:,9) .* s_ibr_nom ;

%-------------------------------------------------------------------
%--------------------- Fixed PI parameters -------------------------
%-------------------------------------------------------------------

% ************* P,Q Measurement LPF *************
w_c     = 0.1 * w_b .* ones(Nibr,1);
t_r     = log(9)/(0.1*w_b);

% ***************** DC Control ******************
tau_dc  = 5e-3 .* ones(Nibr,1)  ; % DC response time-constant

fBW_DC     = 25                   ; % in [Hz]
tau_DC     = 1/(2*pi*fBW_DC)      ;
w_DC       = 1/tau_DC             ;
gamma_DC   = 0.9                  ;
Pm_0       = 1                    ;
c          = -Pm_0*w_b./(Cdc.*Vdc_n) ;

% Define the equation
omega_DC = zeros(Nibr,1);
for n_ibr = 1:Nibr
    equation = @(omega) (w_DC^2 * (2*gamma_DC*omega - c(n_ibr))^2 + omega^4) ./ ...
                         ((2*gamma_DC*omega*w_DC)^2 + (omega^2 - w_DC^2)^2) - 1/2;
    omega0   = 1;                                 % Initial guess for omega  
    %omega_DC(n_ibr,1) = fsolve(equation, omega0); % Solve for omega using fsolve
    omega_DC(n_ibr,1) = 63.4187;
end

Kdc_P      = (2*gamma_DC.*omega_DC.*(Cdc./w_b).*Vdc_n + Pm_0./Vdc_n);
Kdc_I      = omega_DC.^2 .* (Cdc./w_b) .* Vdc_n;

% *************** Current Control ***************
% (Closed loop current control simplified to a 1st order time lag)
fBW_CC  = 1000                   ; % in [Hz]
tau_CC  = 1/(2*pi*fBW_CC)       ;
Ki_P    = L_f./(w_b .* tau_CC)  ;
Ki_I    = R_f./tau_CC           ;

% ***************** PLL Control *****************
tau_lpf_pll  = ones(Nibr,1)/(2*pi*100)  ; % LPF
fBW_PLL      = 0.1                     ; % in [Hz]
tau_PLL      = 1/(2*pi*fBW_PLL)         ;
gamma_dPLL   = 0.9                      ;
omega_dPLL   = (1/tau_PLL)*sqrt( sqrt(4*gamma_dPLL^4 + 4*gamma_dPLL^2 + 2) - (2*gamma_dPLL^2 + +1));
Kpll_P  = -2*gamma_dPLL*omega_dPLL *ones(Nibr,1) ;
Kpll_I  = -omega_dPLL^2            *ones(Nibr,1) ;

% ****************** PQ Control *****************
fBW_PQ     = 0.1                    ; % in [Hz]
tau_PQ     = 1/(2*pi*fBW_PQ)      ;
Kpq_P      = tau_CC/tau_PQ *ones(Nibr,1) ;
Kpq_I      = 1/tau_PQ      *ones(Nibr,1) ;

% **************** Volt. Control ****************
fBW_VC     = 100                  ; % in [Hz]
tau_VC     = 1/(2*pi*fBW_VC)      ;
Kv_P       = C_f./(w_b*tau_VC)    ;
Kv_I       = zeros(Nibr,1)        ;            

%-------------------------------------------------------------------
%----------------  Incidence Matrices Construction  ----------------
%-------------------------------------------------------------------

% Lines Inc. Matrix. 
Be=zeros(Nn,Ne);
for k=1:Ne
    Be(line_from(k),k) =  1;
    Be(line_to(k),k)   = -1;
end

% Shunt Inc. Matrix. 
Bsh=zeros(Nn,Nsh); j=0;
for k = 1:Nn
    if busC(bus_pu(k,1))~=0
        j=j+1;
        Bsh(bus_pu(k,1),j)=1;
    end
end

% Load Inc. Matrix.
Bload=zeros(Nn,Nload);
for l = 1:Nload
    for k=1:Nn
        if load_pu(l,2)==bus_pu(k,1)
            Bload(k,l)=1;
        end
    end
end

% IBR Inc. Matrix.

Bibr=zeros(Nn,Nibr);
for k=1:Nibr
    Bibr(ibrbus(k),k)=1;
end

save('Generated_files\Case_params.mat')