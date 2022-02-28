
function pm30D_dialysis_FIP2021
% Driver function for the pm3 0D dialysis model including whole body lumped parameter model
% blood flow (Heldt), dialyzer unit (Ursino), Baroreflex (Lin), detailed kidney (PM3 developed)
% Copyright (C) 2021. PM3 Lab (Dr SR Kharche).
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% 15 Feb, 2021
% Author: Timothy J. Hunter (thunte27@uwo.ca)
%
% Lead Developer: Sanjay R. Kharche
% Affiliations: Medical Biphysics, Western University and Lawson Health
% Research Institute.
% Part of PM3 platforms. (url: https://kcru.lawsonresearch.ca/research/srk/index.html )
%
%
% References.
% The model is adapted from the following papers.
% The composite model is desciribed in this paper. In the pdfs/ directory, refs refer to references in this paper:
% Heldt, T., Shim, E. B., Kamm, R. D., Mark, R. G., & Massachusetts. (2002). Computational modeling of cardiovascular response to orthostatic stress. Journal of Applied Physiology, 92(3), 1239–1254. https://doi.org/10.1152/japplphysiol.00241.2001
% Ursino, M., Colí, L., Brighenti, C., Chiari, L., De Pascalis, A., & Avanzolini, G. (2000). Prediction of solute kinetics, acid-base status, and blood volume changes during profiled hemodialysis. Annals of Biomedical Engineering, 28(2), 204–216. https://doi.org/10.1114/1.245
% Lin, J., Ngwompo, R. F., & Tilley, D. G. (2012). Development of a cardiopulmonary mathematical model incorporating a baro–chemoreceptor reflex control system. Proceedings of the Institution of Mechanical Engineers, Part H: Journal of Engineering in Medicine, 226(10), 787-803. doi:10.1177/0954411912451823
% Levy, M. N., & Zieske, H. (1969). Autonomic control of cardiac pacemaker activity and atrioventricular transmission. Journal of Applied Physiology, 27(4), 465-470. doi:10.1152/jappl.1969.27.4.465
%
% Note: This MATLAB code is organised in a manner that you can easily port it to the PM3 high order ODE solver.
%
% Our ODE codes need: a) initial conditions; b) the ODE solver; c) outputs for making figures.
%
% system constants.
% global DELTAT;
NEQ     =44; 			% number of ODEs.
NP      = 160;			% number of paramters.
DELTAT  = 0.01; 		% units: seconds. integration time step.
WRITE_INTERVAL = 0.01;  % The interval between writing output to file. As it stands, its 10 ms.
%
%
% declare the state vector and paramter vector.
y0 		= zeros(1,NEQ)	; % simple initial condition, allocation of memory.
yfinal 	= zeros(1,NEQ)	; % PM3 solver provides the solver at every time step given by the user.
%
%
% define the time independent paramters to pass to RHS, also one of p has to be time. I decided
%

n=1;

p 		= zeros(NP,1)	;

p(1)    = 0.0;          % time, seconds.
p(2)    = 150.0;        % Meqic
p(3)    = 150.0;        % Meqex
p(4)    = 500/(60*60);  % % Qf; Ultrafiltration rate. units:mL/s ; ref. Ursino
p(5)    = 0;            % Qinf: infused fluid
p(6)    = 0;            % cud: concentration of urea in dialysate; units: mmol/L
p(7)    = 142;          % cnad: concentration of Na in dialysate; units: mmol/L; ref. ursino 1997 p666
p(8)    = 2;            %ckd: concentration of K in dialysate; units: mmol
p(9)    = 0;            % ccld: concentration of Cl in dialysate; units: mmol
p(10)   = 35;           % chco3d: concentration of HCO3 in dialysate; units: mmol
p(11)   = 0.800;        % t_n_1: Cardiac cycle interval; units (s), baseline duration of cardiac cycle
p(12)   = 12/60;        % baseline respiration rate; units: breath/s; ref. Heldt 2002
p(13)   = 0;            % cuinf: concentration of U in infused fluid; units: mmol/L
p(14)   = 0;            % cnainf: concentration of Na in infused fluid; units: mmol/L
p(15)   = 0;            % ckinf: concentration of K in infused fluid; units: mmol/L
p(16)   = 0;            % cclinf: concentration of Cl in infused fluid; units: mmol/L

p(17)   = 2.67;         % D_s: dialysance (or clearance) of solute; units: ml/s; ref. ursino 2008
p(18)   = 2.67;         % D_U: dialysance (or clearance) of urea ; units: ml/s; ref. ursino 2008
p(19)   = 0.13/60;      % D_HCO3: dialysance (or clearance) of HCO3 ; units: ml/s; ref. ursino 2000
p(20)   = 0.94;         % F_p: plasma water fraction
p(21)   = 0.72;         % F_R: RBC water fraction
p(22)   = 1;            % gamma_u: fraction of red blood cell water that participates in the transfer through the dialyzer
p(23)   = 1;            % R_DU: Donnan ratio for Urea.
p(24)   = 0;            % gamma_Na: fraction of red blood cell water that participates in the transfer through the dialyzer
p(25)   = 0;            % gamma_K: fraction of red blood cell water that participates in the transfer through the dialyzer
p(26)   = 0;            % gamma_Cl: fraction of red blood cell water that participates in the transfer through the dialyzer
p(27)   = 0;            % gamma_HCO3: fraction of red blood cell water that participates in the transfer through the dialyzer
p(28)   = 25;           % k_Na: mass transfer coefficient for Na; units (ml/s)
p(29)   = 0.0704;       % beta_Na: mass transfer coefficient for Na; units (n/a)
p(30)   = 6.67*10^(-2); % k_K: mass transfer coefficient for K; units (ml/s)
p(31)   = 28.2  ;       % beta_K: mass transfer coefficient for K; units (n/a)
p(32)   = 13    ;       % k_U: mass transfer coefficient for U; units (ml/s)
p(33)   = 1     ;       % beta_U: mass transfer coefficient for U; units (n/a)
p(34)   = 4 * 10^(-3) ; % k_f: water exchange coefficient; units (L^2 s^-1 mmol^-1)
p(35)   = 2.45  ;       % E_is: elastance of the interstitial space; units (mmHg/L)
p(36)   = 11    ;       % V_isn: basal volume of interstitial compartments; units (L)
p(37)   = 3.25  ;       % V_pln: basal volume of blood plasma; units (L)
p(38)   = 7.4   ;       % c_ppln: basal protein concentration in plasma; units (g/dl)
p(39)   = 1.37  ;       % c_pisn: basal protein concentration in interstitial compartment; units (g/dl)
p(40)   = 1.05;         % Gibbs Donnan ratio for anions. Ursino 2000, table 1.
p(41)   = 0.95;         % Gibbs Donnan ratio for cations. Ursino 2000, table 1.
p(42)   = 0.2 /60   ;   % eta_hco3: bicarbonate mass transfer coefficient. table 1, Ursino 2000. units: L/s.
p(43)   = 0.03 /60	;   % eta_h: hydrogen ion mass transfer coefficient. table 1, Ursino 2000. units: L/s.
p(44)   = 0.4	;         % g_hco3: bicarbonate equilibrium ratio. table 1, Ursino 2000. units: dimensionless.
p(45)   = 3.5	;         % g_h: hydrogen ion equilibrium ratio. table 1, Ursino 2000. units: dimensionless.
p(46)   = 6.0 / 60;     % etaprime_r: reaction velocity, hco3 buffer. table 1, Ursino 2000. see reaction i. units: L^2/s/mmol
p(47)   = 10.0^(-6.1);  % kprime_a: dissociation constant; p209 of Ursino 2000, right column 3rd para.
p(48)   = 6.0 / 60;     % etaprimeprime_r: reaction velocity, protein buffer. table 1, Ursino 2000. see reaction ii. units: L^2/s/mmol
p(49)   = 10.0^(-7.4);  % kprimeprime_a:dissociation constant; p209 of Ursino 2000, right column 3rd para.
p(50)   = 1.2;          % cco2ic: concentration of CO2 in ic compartment; (in eq 20) Ursino 2000, p209, para 3 on right column.
p(51)   = 1.2;          % cco2ex: concentration of CO2 in ex compartment; (in eq 20) Ursino 2000, p209, para 3 on right column.
p(52)   = 4;            % cpic0: basal protein concentration in intracellular compartment; units mmol/L, Table 1 ursino 2000
p(53)   = -5.97;        % pis0: basal pressure in is compartment; units: mmHg; eq. 10, appendix 1. Ursino 2000. table 1.
p(54)   = 0.01;         % La: arterial capillary permeability; units: mL/mmHg/s; table 1. Ursino 2008, table 1
p(55)   = 0.062;        % Lv: venous capillary permeability; units: mL/mmHg/s; table 1. Ursino 2008, table 1
p(56)   = 0.006;        % Rlo: Resistance of left heart outlet; units: (mmHg*s/ml); Table 2 of Heldt 2002
p(57)   = 3.9;          % Rup1: Resistance of upper body (1); units: (mmHg*s/ml); Table 2 of Heldt 2002
p(58)   = 4.1;          % Rkid1: Resistance of kidney (1); units: (mmHg*s/ml); Table 2 of Heldt 2002
p(59)   = 3.0;          % Rsp1: Resistance of splanchic circulation (1); units: (mmHg*s/ml); Table 2 of Heldt 2002
p(60)   = 3.6;          % Rll1: Resistance of legs (1); units: (mmHg*s/ml); Table 2 of Heldt 2002
p(61)   = 0.23;         % Rup2: Resistance of upper body (2); units: (mmHg*s/ml); Table 2 of Heldt 2002
p(62)   = 0.3;          % Rkid2: Resistance of kidneys (2); units: (mmHg*s/ml); Table 2 of Heldt 2002
p(63)   = 0.18;         % Rsp2: Resistance of splanchic circulation (2); units: (mmHg*s/ml); Table 2 of Heldt 2002
p(64)   = 0.3;          % Rll2: Resistance of legs (2); units: (mmHg*s/ml); Table 2 of Heldt 2002
p(65)   = 0.06;         % Rsup: Resistance of superior vena cava; units: (mmHg*s/ml); Table 2 of Heldt 2002
p(66)   = 0.01;         % Rab: Resistance of abdominal vena cava; units: (mmHg*s/ml); Table 2 of Heldt 2002
p(67)   = 0.015;        % Rinf: Resistance of inferior vena cava; units: (mmHg*s/ml); Table 2 of Heldt 2002
p(68)   = 0.003;        % Rro: Resistance of right heart outlet; units: (mmHg*s/ml); Table 2 of Heldt 2002
p(69)   = 0.08;         % Rp: Resistance of pulmonary arteries; units: (mmHg*s/ml); Table 2 of Heldt 2002
p(70)   = 0.01;         % Rpv: Resistance of pulmonary veins; units: (mmHg*s/ml); Table 2 of Heldt 2002
p(71)   = 15.0;         % Ck: Kidney capacitance; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
p(72)   = 55;           % Csp: Splanchnic capacitance; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
p(73)   = 19.0;         % Cll: Legs venous capacitance; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
p(74)   = 25;           % Cab: Capacitance of abdominal veins; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
p(75)   = 2.0;          % Ca: Capacitance of systemic artery, i.e. aorta; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
p(76)   = 8.0;          % Cup: Capacitance of upper body; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
p(77)   = 2.0;          % Cinf: Capacitance of inferior vena cava; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
p(78)   = 15.0;         % Csup: Capacitance of superior vena cava; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
p(79)   = 4.3;          % Cpa: Capacitance of pulmonary arteries; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
p(80)   = 8.4;          % Cpv: Capacitance of pulmonary veins; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
p(81)   = 0.5;          % Csys_l: left ventricular systolic capacitance units: ml/mmHg; Heldt 2002, page 1243, paragraph 1. % new V6 added 0.64 to account for SNA gain
p(82)   = 10;           % Cdias_l: left ventricular diastolic capacitance units: ml/mmHg; Heldt 2002, page 1243, paragraph 1. % new V6 added 8.45 to account for SNA gain
p(83)   = 1.2;          % Csys_r: right ventricular systolic capacitance units: ml/mmHg; Heldt 2002, page 1243, paragraph 1.
p(84)   = 20;           % Cdias_r: right ventricular diastolic capacitance units: ml/mmHg; Heldt 2002, page 1243, paragraph 1.
p(85)   = 5;            % Q_B: bulk blood flow through the dialyzer; units: ml/s; ref: ursino 1997
p(86)   = DELTAT;       % Simulation time-step; units: sec
% New V6 Below parameters are for Lin Baroreflex
p(87)   = 0;            % t_cardCycleInit: initiation time of the current cardiac cycle
p(88)   = 0;            % SNA
p(89)   = 0;            % PNA
p(90)   = 30;           % Pbco2; units: mmHg
p(91)   = 87;           % Pb02; units: mmHg
p(92)   = 0;            % deltaHR_S
p(93)   = 0;            % deltaHR_V
p(94)   = 0;            % sigma_lv; units: mmHg/ml
p(95)   = 0;            % sigma_rv; units: mmHg/ml
p(96)   = 0;            % delta_sigma_V; units: ml
p(97)   = 0;            % sigma_R; units: mmHg/(ml/s)
p(98)   = 650;          % ZPFV_up: zero-pressure filling volume upper body veins; units: ml;
p(99)   = 150;          % ZPFV_kid;
p(100)  = 1300;         % ZPFV_sp;
p(101)  = 350;          % ZPFV_ll
p(102)  = 250;          % ZPFV_ab
p(103)  = 75;           % ZPFV_inf
p(104)  = 10;           % ZPFV_sup
% Below are parameters of the Baroreflex integrator
p(105)  = 0.001;        % tau_aff: units: seconds
p(106)  = 1;            % G_aff
p(107)  = 0.001;        % tau_c: Same as tau_aff,  units: s
p(108)  = 0.0205;       % S_p
p(109)  = 6;            % PNA_max, units: Hz
p(110)  = 0.6;          % PNA_min, units: Hz
p(111)  = -0.0138;      % S_s
p(112)  = 4;            % SNA_max, units: Hz
p(113)  = 1.12;         % SNA_min, units: Hz
p(114)  = 13.8;         % k1
p(115)  = 0.182;        % k2
p(116)  = 828;          % k3
p(117)  = 1;            % k4
p(118)  = -18.118;      % k5
% Below are parameters of the Baroreflex effector
p(119)  = 90;           % G_k_s0, units: beats/min/Hz
p(120)  = 0.28;         % k_k_s0
p(121)  = 3;            % T_s, units: s
p(122)  = 60;           % G_v0, units: beats/min/Hz 45
p(123)  = 0.4;          % k_v0
p(124)  = 0.5;          % T_v, units: s
p(125)  = 1.5;          % tau_sigma_lv, units: s
p(126)  = 2;            % T_e_lv, units: s
p(127)  = 0.45*1.5;     % G_eff_lv, units: mmHg/ml/Hz % EDIT: original value was 0.45
p(128)  = 1.5;          % tau_sigma_rv, units: s
p(129)  = 2;            % T_e_rv, units: s
p(130)  = 0.282*0.935;  % G_eff_rv, units: mmHg/ml/Hz
p(131)  = 10;           % tau_sigma_V, units: s
p(132)  = 5;            % T_e_V, units: s
p(133)  = 275*13.0;     % G_eff_V, units: ml/Hz EDIT original value was -275
p(134)  = 1.5;          % tau_sigma_R, units: s
p(135)  = 3;            % T_e_R, units: s
p(136)  = 0.20*0.55;    % G_eff_R, units: mmHg/(ml/s)/Hz % EDIT G_eff_R increased from 0.2 to 0.21
p(137)  = 25;           % tau_s, simplified from equation 23, units: s
p(138)  = 0.8;          % tau_v, simplified from eq. 28; units: s
% state vectors for baroreflex
p(139)  = 1;            % x1_P_aff
p(140)  = 1;            % x1_temp1
p(141) = 1;             % x1_deltaHR_s;
p(142) = 1;             % x1_deltaHR_v;
p(143) = 1;             % x1_sigma_lv;
p(144) = 1;             % x1_sigma_rv;
p(145) = 1;             % x1_sigma_V;
p(146) = 1;             % x1_sigma_R;
p(147) = 0;             % sigma_V
p(148) = 85;            % HR0 units: bpm


saveFile = sprintf('output%05d.dat',n);
simTime = 30; % seconds
control = true;
% tic;

%
% initial conditions.
% variables 1 to 5 correspond to odes in Appendix 1, page 584 in Lim et al. paper.
% initial values for volume and k, na, cl were taken from https://www.kumc.edu/AMA-MSS/Study/phys1.htm
%
% Sept 5, 2020. Where is the ANS/sympathetic system/baroreflex component here? (if Ursino did not do, dont add, Jermiah doing already).
% What part of this is the dialyser unit?
y0(1)   = 25;       % Vic, units: L
y0(2)   = 11;       % Vis, units: L
y0(3)   = 3.25;     % Vpl, units: L
y0(4)   = 8.97;     % Pup, units: mmHg
y0(5)   = 0;        % 11.32;  % Pk, units: mmHg
y0(6)   = 10.11*.9; % Psp, units: mmHg
y0(7)   = 12.24*.9; % Pll, units: mmHg
y0(8)   = 3.81*.9;  % Pab, units: mmHg
y0(9)   = -5;       % Pth, units: mmHg. % is pth a variable, or a constant? 5 Sept. 2020.
y0(10)  = 10.0;     % Cl, units: ml/mmHg
y0(11)  = 20.0;     % Cr, units: ml/mmHg
y0(12)  = 12.78*.9; % Pl, left ventricular pressure, units: mmHg
y0(13)  = 91*.9;    % Pa, aortic pressure, units: mmHg
y0(14)  = 3.47*.9;  % Psup, units: mmHg
y0(15)  = 3.31*.9;  % Pinf, units: mmHg
y0(16)  = 2.60*.9;  % Pr, units: mmHg
y0(17)  = 15.66*.9; % Ppa, units: mmHg
y0(18)  = 12.99*.9; % Ppv, units: mmHg
y0(19)  = 100;      % Muic, units: mmol
y0(20)  = 250.0;    % Mnaic, units: mmol
y0(21)  = 3535.0;   % Mkic, units: mmol
y0(22)  = 84.0;     % Mclic, units: mmol
y0(23)  = 10.0;     % MHco3ic, units: mmol
y0(24)  = 100.0;    % Mhic, units: mmol
y0(25)  = 0.0;      % Mpic, units: mmol
y0(26)  = 55.0;     % Muex, units: mmol
y0(27)  = 2130.0;   % Mnaex, units: mmol
y0(28)  = 75.0;     % Mkex, units: mmol
y0(29)  = 1470.0;   % Mclex, units: mmol
y0(30)  = 100.0;    % Mhco3ex, units: mmol
y0(31)  = 100.0;    % Mhex, units: mmol
y0(32)  = 0.0;      % Mpex, units: mmol
y0(33)  = 11.32;    % Pk, units: mmHg
y0(34)  = 11.32;    % Pk, units: mmHg
y0(35)  = 11.32;    % Pk, units: mmHg
y0(36)  = 11.32;    % Pk, units: mmHg
y0(37)  = 11.32;    % Pk, units: mmHg
y0(38)  = 11.32;    % Pk, units: mmHg
y0(39)  = 11.32;    % Pk, units: mmHg
y0(40)  = 11.32;    % Pk, units: mmHg
y0(41)  = 11.32;    % Pk, units: mmHg
y0(42)  = 11.32;    % Pk, units: mmHg
y0(43)  = 11.32;    % Pk, units: mmHg
y0(44)  = 11.32;    % Pk, units: mmHg


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% code for one time step. ================================================================
% single time step iterator. this is what forms the for loop.
fileID = fopen(saveFile,'w'); % This may be changed to 'w'
% start time and end time. You do this for each time step, as y(:,end) is what you put into output.
tspan = [0; DELTAT];
% options need defining once. 5 sept. 2020.
% the basic tolerances used by the PM3 solver are absolute and relative tolerances. These two tolerances have standard definitions.
options = odeset('RelTol',1e-6,'AbsTol',1e-6);

% new
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HR_0 = p(148); % units: beats/min

SNA_buffer = zeros(p(132)/DELTAT,1) + 3;
PNA_buffer = zeros(p(124)/DELTAT,1) + 3;

% f = waitbar(0, 'Simulation progress');
% s = sprintf("Time remaining %.0f s\n", toc/(0*DELTAT/simTime)*(1 - (0*DELTAT/simTime)));
% f = waitbar(0*DELTAT/simTime, s);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i_iter = 1:1:simTime/DELTAT
for i_iter = 1:1:simTime/DELTAT

	% set up the time to pass to RHS.
    p(1) = p(1) + DELTAT; % in seconds.
    if ~control
        p = modParam(p);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% new
    %%%%% Baroreflex calculations %%%%%
    [p] = pm30D_baroreflex_integration_FIP2021(y0, p);

    % update SNA and PNA buffers
    for i = 1:length(SNA_buffer) - 1
        SNA_buffer(length(SNA_buffer) + 1 - i) = SNA_buffer(length(SNA_buffer) - i);
    end
    SNA_buffer(1) = p(88);

    for i = 1:length(PNA_buffer) - 1
        PNA_buffer(length(PNA_buffer) + 1 - i) = PNA_buffer(length(PNA_buffer) - i);
    end
    PNA_buffer(1) = p(89);
    % done updating SNA and PNA buffers

    % calculate effector site ceofficients
    [p] = pm30D_baroreflex_effector_FIP2021(p, SNA_buffer, PNA_buffer);

    % calculate heart rate modification
    deltaHR = 19.64*p(92) - 17.95*p(93) - 1.225 * p(92)^2 + 1.357 * p(93)^2 - 1.523 * p(92) * p(93);% page 468 Levy and Zeiske
    HR = HR_0 + deltaHR; % eq. 18

    p(11) = 60 / HR;

    % once cardiac cycle completes, update init time p(87)
    if p(1) >= p(87) + p(11)
        p(87) = p(1);
    end
    %%%%% end of baroreflex %%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % s = sprintf("Time remaining %.0f s\n", toc/(i_iter*DELTAT/simTime)*(1 - (i_iter*DELTAT/simTime)));
    % waitbar(i_iter*DELTAT/simTime, f, s);
    % tspan is always 0 to DELTAT. There are forcing terms in the RHS, and you have to pass time as a parameter.
    [t y] = ode15s(@(t,y) pm30D_dialysis_RHS_FIP2021(t,y,p),tspan,y0,options);
    y0 = y(end,1:NEQ); % init conditions and for output.

    % write to file.
	if(mod(i_iter*DELTAT, WRITE_INTERVAL) == 0) % deltat is 10 ms.
		fprintf(fileID,'%f\t', p(1));
		for nNum=1:1:NEQ
            fprintf(fileID,'%f\t',y0(nNum));
        end;
		for nNum=2:1:NP
            fprintf(fileID,'%f\t',p(nNum));
        end;
		fprintf(fileID,"\n");
	end; % end of writing.
    %
end; % end of time loop. Use semi-colons for consistency regardless of what MATLAB documentation says.
fclose(fileID);
% close(f);
%
% ==================================================================================
% -------------------------end of pm30D_dialysis main driver function.-------------------------------------------------------------------------.
%
end

function dydt = pm30D_dialysis_RHS_FIP2021(t,y,p)
%
% 29 January, 2021
% Author: Timothy Hunter
% Lead Developer: Sanjay Kharche
% Part of PM3 Lab
%


dydt = zeros(size(y)); % derivative vector is same size as state vector.

%%% ALL PASSED PARAMETERS ARE LISTED BELOW %%%
%
% the t in the function dydt = ... is something between 0 and DELTAT. It is used by the implicit solver.
%
% this is the physical time that goes into elastances and all time dependent RHS formulas.
% On p1241 of Heldt 2002, the time used in the elastances starts at LV contraction. For development,
% I am assuming that the period of the cardiac cycle is 1 s. This needs revision.
% DELTAT      = 0.01; 		% units: seconds. integration time step. Assign DELTAT in ONLY ONE PLACE, pass from p or make global
T_n_1       = p(11); %  may be passed as a parameter. Why it must be passed as par? Sept 5, 2020.
% new
loc_t       = p(1) + t - p(87); % 1.0 is the assumed heart period at development, revise this.
%
% Amount of (other) solute in compartments, to be measured in patients
% ref: table 2 Ursino & Innocenti 2008
M_eqic      = p(2);
M_eqex      = p(3);
QF          = p(4);     % Ultrafiltration rate. units:mL/s
Qinfused    = p(5);     % infusion rate (Qinfused is used in place of Qinf; which is already flow in inferior vena cava.) units:mL/s
Q_B      = p(85);
%
% from Ursino 2000, p207
% Dialysate concentrations of solutes TEMPORARY (may be passed from p)
cud         = p(6);     % urea concentration in dialysate. placeholder. p207 of Ursino 2000.
cnad        = p(7);     % sodium concentration in dialysate. placeholder. p207 of Ursino 2000.
ckd         = p(8);     % potassium concentration in dialysate. placeholder. p207 of Ursino 2000.
ccld        = p(9);     % chlorine concentration in dialysate. placeholder. p207 of Ursino 2000.
chco3d      = p(10);    % carbonate concentration in dialysate. placeholder. p207 of Ursino 2000.
respRate    = p(12);    % respiratory rate; units: (breaths/s)
%
% from ursino 2008
cuinf       = p(13);    % urea concentration in infusion, units: mmol/L
cnainf      = p(14);    % sodium concentration in infusion, units: mmol/L
ckinf       = p(15);    % potassium concentration in infusion, units: mmol/L
cclinf      = p(16);    % chlorine concentration in infusion, units: mmol/L

%%% ALL CONSTANT PARAMETERS ARE LISTED BELOW %%%
% Parameters are listed in the order they are called by the ODEs
%
% Fluid and Solute Kinetics
% Table 2, Ursino 2008 (unless stated otherwise).
D_s 		= p(17);        % Na, Cl, K clearance/dialysance. Table 2, Ursino 2008. units: mL/s.
D_u 		= p(18);        % urea clearance/dialysance. Table 2, Ursino 2008. units: mL/s.
D_hco3      = p(19);    % carbonate clearance/dialysance. Table 1, Ursino 2000. units: L/s.
F_p         = p(20);    % F_p: plasma water fraction
F_R         = p(21);    % F_R: RBC water fraction
gamma_U     = p(22);    % gamma_U: fraction of red blood cell water that participates in the transfer through the dialyzer
R_DU        = p(23);    % R_DU: Donnan ratio for Urea.
gamma_Na    = p(24);    % gamma_Na: fraction of red blood cell water that participates in the transfer through the dialyzer
gamma_K     = p(25);    % gamma_K: fraction of red blood cell water that participates in the transfer through the dialyzer
gamma_Cl    = p(26);    % gamma_Cl: fraction of red blood cell water that participates in the transfer through the dialyzer
gamma_HCO3  = p(27);    % gamma_HCO3: fraction of red blood cell water that participates in the transfer through the dialyzer
k_Na        = p(28);    % k_Na: mass transfer coefficient for Na; units (ml/s)
beta_Na     = p(29);    % beta_Na: mass transfer coefficient for Na; units (n/a)
k_K         = p(30);    % k_K: mass transfer coefficient for K; units (ml/s)
beta_K      = p(31);    % beta_K: mass transfer coefficient for K; units (n/a)
k_U         = p(32);    % k_U: mass transfer coefficient for U; units (ml/s)
beta_U      = p(33);    % beta_U: mass transfer coefficient for U; units (n/a)
k_f         = p(34);    % k_f: water exchange coefficient; units (L^2 s^-1 mmol^-1)
E_is        = p(35);    % E_is: elastance of the interstitial space; units (mmHg/L)
V_isn       = p(36);    % V_isn: basal volume of interstitial compartments; units (L)
V_pln       = p(37);    % V_pln: basal volume of blood plasma; units (L)
c_ppln      = p(38);    % c_ppln: basal protein concentration in plasma; units (g/dl)
c_pisn      = p(39);    % c_pisn: basal protein concentration in interstitial compartment; units (g/dl)

alphaa      = p(40);    % Gibbs Donnan ratio for anions. Ursino 2000, table 1.
alphac      = p(41);    % Gibbs Donnan ratio for cations. Ursino 2000, table 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Table 1 from Ursino 2000 p 205
eta_hco3    = p(42);    % eta_hco3: bicarbonate mass transfer coefficient. table 1, Ursino 2000. units: L/s.
eta_h       = p(43);    % eta_h: hydrogen ion mass transfer coefficient. table 1, Ursino 2000. units: L/s.
g_hco3  	  = p(44);    % g_hco3: bicarbonate equilibrium ratio. table 1, Ursino 2000. units: dimensionless.
g_h 		    = p(45);    % g_h: hydrogen ion equilibrium ratio. table 1, Ursino 2000. units: dimensionless.
%
% TEMPORARY insert references
etaprime_r 	= p(46);    % etaprime_r: reaction velocity, hco3 buffer. table 1, Ursino 2000. see reaction i. units: L^2/s/mmol
kprime_a 		= p(47);    % kprime_a: dissociation constant; p209 of Ursino 2000, right column 3rd para.
etaprimeprime_r = p(48);% etaprimeprime_r: reaction velocity, protein buffer. table 1, Ursino 2000. see reaction ii. units: L^2/s/mmol
kprimeprime_a 	= p(49);% kprimeprime_a:dissociation constant; p209 of Ursino 2000, right column 3rd para.

cco2ic 			= p(50);    % cco2ic: concentration of CO2 in ic compartment; (in eq 20) Ursino 2000, p209, para 3 on right column.
cco2ex 			= p(51);    % cco2ex: concentration of CO2 in ex compartment; (in eq 20) Ursino 2000, p209, para 3 on right column.

cpic0       = p(52);    % cpic0: basal protein concentration in intracellular compartment; units mmol/L, Table 1 ursino 2000

% Capillary pressure alues from Lim 2008 p584 paragraph 1
% units (mmHg)
pac1    = y(4)        ; pac2    = y(6)       ;  pac3     = y(7);
pvc1    = y(4) - y(14); pvc2    = y(6) - y(8);  pvc3     = y(7) - y(8);
%
pis0        = p(53);    % pis0: basal pressure in is compartment; units: mmHg; eq. 10, appendix 1. Ursino 2000. table 1.
La          = p(54);    % La: arterial capillary permeability; units: mL/mmHg/s; table 1. Ursino 2008
Lv          = p(55);    % Lv: venous capillary permeability; units: mL/mmHg/s; table 1. Ursino 2008, table 1
%
% new V6
sigma_lv = p(94);       % units: mmHg/ml
sigma_rv = p(95);       % units: mmHg/ml
delta_sigma_V = p(96);  % units: ml
V_gain = delta_sigma_V / 2785; % This is the change in venous volume / total venous volume
sigma_R = p(97);        % units: mmHg/(ml/s)
ZPFV_up = p(98);        % zero-pressure filling volume upper body veins; units: ml
ZPFV_kid = p(99);       % zero-pressure filling volume kidneys ;units: ml
ZPFV_sp = p(100);       % zero-pressure filling volume splanchnic ;units: ml
ZPFV_ll = p(101);       % zero-pressure filling volume lower legs ;units: ml
ZPFV_ab = p(102);       % zero-pressure filling volume abdominal vena cava ;units: ml
ZPFV_inf = p(103);      % zero-pressure filling volume inf. vena cava ;units: ml
ZPFV_sup = p(104);      % zero-pressure filling volume sup. vena cava ;units: ml


% Table 2 of Heldt 2002. AJP. Resistances. I used same notation as Heldt
% Resistances for inlet and outlet of each vascular compartment (mmHg*s/ml)
Rlo         = p(56);    % Rlo: Resistance of left heart outlet; units: (mmHg*s/ml); Table 2 of Heldt 2002
% new V6 : arteriolar resistances are multiplied by gain value
Rup1        = p(57) * sigma_R; % Rup1: Resistance of upper body (1); units: (mmHg*s/ml); Table 2 of Heldt 2002
Rkid1       = p(58) * sigma_R; % Rkid1: Resistance of kidney (1); units: (mmHg*s/ml); Table 2 of Heldt 2002
Rsp1        = p(59) * sigma_R; % Rsp1: Resistance of splanchic circulation (1); units: (mmHg*s/ml); Table 2 of Heldt 2002
Rll1        = p(60) * sigma_R; % Rll1: Resistance of legs (1); units: (mmHg*s/ml); Table 2 of Heldt 2002
Rup2        = p(61);    % Rup2: Resistance of upper body (2); units: (mmHg*s/ml); Table 2 of Heldt 2002
Rkid2       = p(62);    % Rkid2: Resistance of kidneys (2); units: (mmHg*s/ml); Table 2 of Heldt 2002
Rsp2        = p(63);    % Rsp2: Resistance of splanchic circulation (2); units: (mmHg*s/ml); Table 2 of Heldt 2002
Rll2        = p(64);    % Rll2: Resistance of legs (2); units: (mmHg*s/ml); Table 2 of Heldt 2002
Rsup        = p(65);    % Rsup: Resistance of superior vena cava; units: (mmHg*s/ml); Table 2 of Heldt 2002
Rab         = p(66);    % Rab: Resistance of abdominal vena cava; units: (mmHg*s/ml); Table 2 of Heldt 2002
Rinf        = p(67);    % Rinf: Resistance of inferior vena cava; units: (mmHg*s/ml); Table 2 of Heldt 2002
Rro         = p(68);    % Rro: Resistance of right heart outlet; units: (mmHg*s/ml); Table 2 of Heldt 2002
Rp          = p(69);    % Rp: Resistance of pulmonary arteries; units: (mmHg*s/ml); Table 2 of Heldt 2002
Rpv         = p(70);    % Rpv: Resistance of pulmonary veins; units: (mmHg*s/ml); Table 2 of Heldt 2002
%
% Appendix 2. of Heldt 2002 AJP. Table 3 and text on p. 1242-1244.
% Capacitances of each vascular compartment
% (Capacitances are with a capital C)
% units: mL/mmHg
Ck          = p(71);    % Ck: Kidney capacitance; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
Csp         = p(72);    % Csp: Splanchnic capacitance; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
Cll         = p(73);    % Cll: Legs venous capacitance; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
Cab         = p(74);    % Cab: Capacitance of abdominal veins; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
Ca          = p(75);    % Ca: Capacitance of systemic artery, i.e. aorta; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
Cup         = p(76);    % Cup: Capacitance of upper body; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
Cinf        = p(77);    % Cinf: Capacitance of inferior vena cava; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
Csup        = p(78);    % Csup: Capacitance of superior vena cava; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
Cpa         = p(79);    % Cpa: Capacitance of pulmonary arteries; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
Cpv         = p(80);    % Cpv: Capacitance of pulmonary veins; units: mL/mmHg; Heldt 2002 AJP. p. 1242.
%
% Heldt 2002, page 1243, paragraph 1
% left ventrical capacitances.
Csys_l      = p(81);    % Csys_l: left ventricular systolic capacitance units: ml/mmHg; Heldt 2002, page 1243, paragraph 1.
Cdias_l     = p(82);    % Csys_l: left ventricular diastolic capacitance units: ml/mmHg; Heldt 2002, page 1243, paragraph 1.
%
% Heldt 2002, page 1243, paragraph 1
% right ventrical capacitances.
Csys_r      = p(83);    % Csys_r: right ventricular systolic capacitance units: ml/mmHg; Heldt 2002, page 1243, paragraph 1.
Cdias_r     = p(84);    % Cdias_r: right ventricular diastolic capacitance units: ml/mmHg; Heldt 2002, page 1243, paragraph 1.
%
DELTAT = p(86);





%
%%% ALL ALGEBRAIC EQNS ARE LISTED BELOW %%%
%
% eq. 16 of Ursino 2000 p207
% solute concentration in Intracellular compartment
% Check if abs(y(1)) < 10^(-16)
cuic        = y(19) 	/ y(1); % urea concentration. 		Muic is a state variable. eq. 16 of Ursino 2000.
cnaic       = y(20) 	/ y(1); % sodium concentration. 	Mnaic is a state variable. eq. 16 of Ursino 2000.
ckic        = y(21) 	/ y(1); % potassium concentration. 	Mkic is a state variable. eq. 16 of Ursino 2000.
cclic       = y(22) 	/ y(1); % chlorine concentration. 	Mclic is a state variable. eq. 16 of Ursino 2000.
chco3ic 	  = y(23) 	/ y(1); % carbonate concentration. 	Mhco3lic is a state variable. eq. 16 of Ursino 2000.
chic        = y(24) 	/ y(1); % hydrogen concentration. 	Mhlic is a state variable. eq. 16 of Ursino 2000.
cpic        = y(25) 	/ y(1); % protein concentration. 	Mpic is a state variable. eq. 16 of Ursino 2000.
%

% eq. 9 of Ursino 2000 p206.
% Volume of extracellular fluid (l)
Vex 		= y(2) + y(3);  % for equation 5 in Lim et al. Then Urisino 2000 gives formula eq. 9.
%
% This approximation (csex = cspl = csis) proposed by ursino and innocenti
% 2008 eq. 24
cuex        = y(26) / Vex; cupl     = cuex  ; cuis      = cuex;
cnaex       = y(27) / Vex; cnapl    = cnaex ; cnais     = cnaex;
ckex        = y(28) / Vex; ckpl     = ckex  ; ckis      = ckex;
cclex       = y(29) / Vex; cclpl    = cclex ; cclis     = cclex;
chco3ex     = y(30) / Vex; chco3pl  = chco3ex; chco3is  = chco3ex;
chex        = y(31) / Vex; chpl     = chex  ; chis      = chex;
cpex        = y(32) / Vex; cppl     = cpex  ; cpis      = cpex;

%
% eq. 2-3 of Ursino 2000 p205.
% Osmotic activity of fluid compartment (mmol/l)
cic 		= (y(19) + y(20) + y(21) + M_eqic) / y(1) ; 	% eq 2, page 205 Ursino 2000; Intracellular osmotic activity
cis 		= (y(26) + y(27) + y(28) + M_eqex) / Vex ; 		% eq 3, page 205 Ursino 2000; Interstitial osmotic activity
%
c_ppl = c_ppln * (V_pln / y(3)); % Ursino & Innocenti 2008 eq. 39
c_pis = c_pisn * (V_isn / y(2)); % Ursino & Innocenti 2008 eq. 38

% eq. 10, appendix 1 of Lim 2008 p 584
% pressure in the interstitial compartment (mmHg)
pis 		= E_is * ( y(2) - V_isn ) + pis0;
%
% eq. 11, appendix 1 of Lim 2008 p 584
% osmotic pressure in the blood plasma (mmHg)
pipl 		= 2.1 * c_ppl + 0.16 * c_ppl*c_ppl + 0.009 * c_ppl * c_ppl * c_ppl; 		% eq. 11, appendix 1.
%
% eq. 12, appendix 1 of Lim 2008 p 584
% osmotic pressure of the interstitial space
piis 		= 2.8 * c_pis + 0.18 * c_pis*c_pis + 0.012 * c_pis * c_pis * c_pis; 		% eq. 12, appendix 1.
%
% eq. 8, appendix 1 of Lim 2008 p 584
Fa1 		= La * (pac1 - pis - pipl + piis ); 							% eq. 8, appendix 1, for n = 1.
Fa2 		= La * (pac2 - pis - pipl + piis ); 							% eq. 8, appendix 1, for n = 2.
Fa3 		= La * (pac3 - pis - pipl + piis ); 							% eq. 8, appendix 1, for n = 3.
% eq. 9, appendix 1 of Lim 2008 p 584
Rv1 		= Lv * (pis - pvc1 + pipl - piis ); 							% eq. 9, appendix 1, for n = 1.
Rv2 		= Lv * (pis - pvc2 + pipl - piis ); 							% eq. 9, appendix 1, for n = 2.
Rv3 		= Lv * (pis - pvc3 + pipl - piis ); 							% eq. 9, appendix 1, for n = 3.

% eq. 6, appendix 1. of Ursino 2000.
% Fa is the rate at which fluid is filtered at arterial capillaries.
% Units: ml/s
Fa 			= Fa1 + Fa2 + Fa3; 									% units (ml/s), eq. 6, appendix 1. of Ursino 2000.
% eq. 7, appendix 1. of Ursino 2000.
% Rv is the rate at which fluid is filtered at venous capillaries.
% Units: ml/s
Rv 			= Rv1 + Rv2 + Rv3; 									% units (ml/s), eq 7, appendix 1. of Ursino 2000.
%
%
% from Lim 2008, unless it is from Heldt 2002. Notation as close to the two papers as possible.
% Flows into and out of each vascular compartment (l/s?)
% This structure works as a valve; allowing positive flow but not negative
if (y(12) < y(18))
    Qli 	= (y(18) - y(12))	/Rpv	;
else
    Qli     = 0;
end;
if (y(13) < y(12))
    Qlo 	= (y(12) - y(13))	/Rlo	;
else
    Qlo     = 0;
end;
Qupi 	= (y(13) - y(4))	/Rup1	;
if (y(14) < y(4))
    Qupo 	= (y(4)   - y(14))	/Rup2	;
else
    Qupo    = 0;
end;
if (y(16) < y(14))
    Qsup 	= (y(14) - y(16))	/Rsup	;
else
    Qsup    = 0;
end;
Qsp1 	= (y(13) - y(6))	/Rsp1	;
Qsp2 	= (y(6)   - y(8))	/Rsp2	;
% Qk1 	= (y(13) - y(5))	/Rkid1	;
% Qk2 	= (y(5)   - y(8))	/Rkid2	;
Qll1 	= (y(13) - y(7))	/Rll1	;
if (y(8) < y(7))
    Qll2	= (y(7)   - y(8))	/Rll2	;
else
    Qll2    = 0;
end;
Qab 	= (y(8) - y(15))	/Rab	;
if (y(16) < y(15))
    Qinf 	= (y(15) - y(16))	/Rinf	;
else
    Qinf    = 0;
end;
if (y(17) < y(16))
    Qro 	= (y(16) - y(17))	/Rro	;
else
    Qro     = 0;
end;
Qpa 	= (y(17) - y(18))	/Rp		; % TEMPORARY there is no Qpa as printed by Lim et al. double check this, and work out where it goes.
%
%
%================== Pumping function for the heart ==================
%
Ts          = 0.3 * sqrt(T_n_1); % systolic time interval.
%
% Heldt 2002 p.1241, paragraph 4
% Left ventricular elastances
Edias_l     = 1/Cdias_l;  % Formula from Heldt 2002 p.1241, paragraph 4.
Esys_l      = (1/Csys_l) * sigma_lv; % Formula from Heldt 2002 p.1241, paragraph 4.
%
if(0<= loc_t && loc_t <= Ts)
    El      = Edias_l + ((Esys_l - Edias_l)/ 2.0) * (1.0 - cos(pi * loc_t / Ts));
elseif(Ts < loc_t) && (loc_t <= 1.5 *Ts)
    El      = Edias_l + ((Esys_l - Edias_l)/ 2.0) * (1.0 + cos(2 * pi * (loc_t - Ts) / Ts));
else
    El = Edias_l;
end;
Cl          = 1 / El;
%
% Heldt 2002 p.1241, paragraph 4
% Right ventricular elastances
Edias_r     = 1/Cdias_r;  % Formula from Heldt 2002 p.1241, paragraph 4.
Esys_r      = (1/Csys_r) * sigma_rv; % Formula from Heldt 2002 p.1241, paragraph 4.
%
if (0<= loc_t && loc_t <= Ts)
    Er      = Edias_r + (Esys_r - Edias_r)/ 2.0 * (1.0 - cos(pi * loc_t / Ts));
elseif (Ts < loc_t) && (loc_t <= 1.5 *Ts)
    Er      = Edias_r + (Esys_r - Edias_r)/ 2.0 * (1.0 + cos(2 * pi * (loc_t - Ts) / Ts));
else
    Er      = Edias_r;
end;
Cr          = 1 / Er;
%
%===================================================================
%
% eq. 12 from Ursino 2000 p206 (unless stated otherwise)
% mass transfer rate from ic compartment to is compartment
phi_k       = - k_K * ( ckic - beta_K * ckis ) ;            % Ursino 2008, eq 26. amount of solute exchanged at cellular membrane per unit time. units:
phi_na      = - k_Na * ( cnaic - beta_Na * cnais ) ;        % Ursino 2008, eq 26. amount of solute exchanged at cellular membrane per unit time. units:
phi_u       = - k_U * ( cuic - beta_U * cuis ) ;            % Ursino 2008, eq 26. amount of solute exchanged at cellular membrane per unit time. units:
phi_hco3    = - eta_hco3 * ( chco3ic - g_hco3 * chco3is ) ; % Ursino 2000, eq 12. amount of solute exchanged at cellular membrane per unit time. units:
phi_h       = - eta_h * ( chic - g_h * chis ) ;             % Ursino 2000, eq 12. amount of solute exchanged at cellular membrane per unit time. units:
phi_p       = 0.0;                                          % Ursino 2000. p207, second paragraph.
phi_cl      = phi_na + phi_k + phi_h - phi_hco3;            % Ursino 2000, eq 13.
%
% eq 11 from Ursino 2000. These are eq 20 to 22 p207 of Ursino 2000. Algebriac equations.
Rhco3ic 	  = etaprime_r * ( kprime_a * cco2ic - chco3ic * chic );       % eq 20
chpic       = cpic0 - cpic;                                              % ref ursino 2000 p. 207 last paragraph
Rpic 		    = etaprimeprime_r * (kprimeprime_a * chpic - cpic * chic );  % eq 21. Dont know chpic yet (4 June 2020).
Rhic 		    = Rhco3ic + Rpic;                                            % eq 22.
%
% Rsex values are analogous to Rsic values. Ursino 2000 p.
Rhco3is 	  = etaprime_r * ( kprime_a * cco2ex - chco3is * chis );       % eq 20
%
% for chpis in eq 21, see bottom para p207 of Ursino 2000.
% Protein activity is assumed to be mainly in plasma
chpis       = c_pis - cpis;                                              % ref ursino 2000 p. 207 last paragraph
Rpis 		    = etaprimeprime_r * (kprimeprime_a * chpis - cpis * chis );  % etaprimeprime_r * (kprimeprime_a * chpis - cpex * chex ); % eq 21. Dont know chpic yet (4 June 2020).
Rhis 		    = Rhco3is + Rpis;                                            % eq 22.
Rhco3pl 	  = etaprime_r * ( kprime_a * cco2ex - chco3pl * chpl );       % eq 20
%
% for chpic in eq 21, see bottom para p207 of Ursino 2000.
chppl       = c_ppl - cppl;                                              % ref ursino 2000 p. 207 last paragraph
Rppl 		    = etaprimeprime_r * (kprimeprime_a * chppl - cppl * chpl );  % eq 21. Dont know chpic yet (4 June 2020).
Rhpl 		    = Rhco3pl + Rppl;                                            % eq 22.
%

Rhco3ex     = Rhco3pl + Rhco3is; Rpex = Rppl + Rpis; Rhex = Rhpl + Rhis;
%
V_rc        = 1.3;                                                       % red blood cell volume, units (L) ursino + innocenti 2008
V           = V_rc + y(3);                                               % whole blood volume
Hct         = V_rc / V;

% eq 29 from Ursino and Innocenti 2008 p888
Q_eK        = Q_B*(F_p * (1 - Hct) + F_R * gamma_K * alphac);
Q_eNa       = Q_B*(F_p * (1 - Hct) + F_R * gamma_Na * alphac);
Q_eU        = Q_B*(F_p * (1 - Hct) + F_R * gamma_U * R_DU);
Q_eHCO3     = Q_B*(F_p * (1 - Hct) + F_R * gamma_HCO3 * alphac);
Q_eCl       = Q_B*(F_p * (1 - Hct) + F_R * gamma_Cl * alphaa);

% eq 28 from Ursino and Innocenti 2008 p888
% convective and diffusive transport to dialyzer
if (QF ~= 0 && Q_B ~= 0)
    J_k     = (D_s * (1 - QF / Q_eK) + QF) * ckex - D_s * (1 - QF / Q_eK) * ckd ;
    J_na    = (D_s * (1 - QF / Q_eNa) + QF) * cnaex - D_s * (1 - QF / Q_eNa) * cnad ;
    J_u     = (D_u * (1 - QF / Q_eU) + QF) * cuex - D_u * (1 - QF / Q_eU) * cud ;
    J_hco3  = (D_hco3 * (1 - QF / Q_eHCO3) + QF) * chco3ex - D_hco3 * (1 - QF / Q_eHCO3) * chco3d ;
    J_cl    = (D_s * (1 - QF / Q_eCl) + QF) * cclex - D_s * (1 - QF / Q_eCl) * ccld ;
else
    J_k = 0; J_na = 0; J_u = 0; J_hco3 = 0; J_cl = 0;
end;
%
J_h     = 0.0; % for hydrogen J is 0.
J_p     = 0.0; % for proteins, J is 0.
%
%
%%% ALL ODEs ARE LISTED BELOW %%%
%
[dydtKidney, Qki, Qko] = pm30D_detailedKidney_FIP2021(y,Rkid1,Rkid2,Ck); % new V5; Implementation of r+l kidneys, 6 lobes each

% in Urisino eq 2, Oic is the intracellular osmotic activity. units: probably mEq/L or mmol/L unless the 0.93 has
% units itself. Ois is the interstitial osmotic activity.
dydt(1) 	=   k_f*(cic - cis);                        % equation 1, appendix 1, Vic, variable. Vic is intracellular fluid volume. fluid potentially means plasma. According to Ursino 2000. units in Ursino 2000: L (liters).
dydt(2) 	= - k_f*(cic - cis) + (Fa - Rv)/1000;       % eq. 2, appendix 1, Vis. Vis is interstial fluid volume. Ursino 2000. units: L.
dydt(3) 	= ( -(Fa - Rv) - QF + Qinfused)/1000; 	    % eq. 3, appedix 1, Vpl. Vpl is plasma volume. Ursino 2000. unit: L.
%
% y(4) is Pup
dydt(4)     = (Qupi - Qupo + V_gain * ZPFV_up) / Cup; % y(4) is Pup, Cl is non-zero. eq. 22, appendix 2, p585.

% y(5) is Pk, pressure at kidneys (not in use)

% y(6) is Psp, splanchic pressure at inlet.
dydt(6)     = (Qsp1 - Qsp2 + V_gain * ZPFV_sp) / Csp; % y(6) is Psp, Csp is non-zero. eq. 24, appendix 2, p585.
% y(7) is Pll, legs pressure.
dydt(7)     = (Qll1 - Qll2 + V_gain * ZPFV_ll) / Cll; % y(7) is Pll, Cll is non-zero. eq. 25, appendix 2.
% y(8) is Pab.
dydt(8)     = (Qko + Qsp2 + Qll2 - Qab + V_gain * ZPFV_ab) / Cab;  % y(8) is Pab. eq. 26, appendix 2. TEMPORARY Fa was added to account for volume filtered by capilleries
% y(9) is Pth: thoracic pressure.
dydt(9)     = (2*pi*respRate)*cos((2*pi*respRate)*p(1) + t);       % This varies with respiratory variation, mean value taken from Heldt thesis 2.3 p.47 resprate units are (breath/s)
%
%
dydt(10)    = (Cl - y(10))/DELTAT;   %  Cl is only a time dependent parameter, therefore dydt is calculated with a pseudo-euler method
dydt(11)    = (Cr - y(11))/DELTAT;   %  Cr is only a time dependent parameter, therefore dydt is calculated with a pseudo-euler method
%
% y(12) is Pl, LV pressure.
dydt(12)    = (1.0/y(10)) * ( (y(9) - y(12)) * dydt(10) + Qli - Qlo ) + dydt(9);  % eq 20, appendix 2. y(12) is Pl.
% y(13) is Pa.
dydt(13)    = (1.0/Ca) * (Qlo + ( -(Fa - Rv) - QF + Qinfused) - (Qupi + Qsp1 + Qll1 + Qki)) + 0.333*dydt(9);    % eq 21, appendix 2. y(13) is Pa. I am guessing that the Qup in appendix 2 is actually Qupi.
% y(14) is Psup.
dydt(14)    = ((Qupo) - Qsup + V_gain * ZPFV_sup) / Csup + dydt(9);  % y(14) is Psup, eq. 27 appendix 2.
% y(15) is Pinf.
dydt(15)    = (Qab - Qinf + V_gain * ZPFV_inf) / Cinf + dydt(9);     % y(15) is Pinf, eq. 28 appendix 2.

% y(16) is Pr (Prv).
dydt(16)    = (1.0/y(11)) * ( (y(9) - y(16)) * dydt(11) + Qinf + Qsup - Qro ) + dydt(9);      % eq 29, appendix 2. y(16) is Pr. I dont think there is any Cri, only Cr. Double check all the same.

% y(17) is Ppa.
dydt(17)    = (Qro - Qpa) / Cpa + dydt(9);  % y(17) is Ppa, eq. 30 appendix 2.
% y(18) is Ppv.
dydt(18)    = (Qpa - Qli) / Cpv + dydt(9);  % y(18) is Ppv, eq. 31 appendix 2.
%
% for eq 11, Ursino 2000
% State variables 19-25 are Intracellular Solute mass
dydt(19) = phi_u 		;                % urea 		   eq 26 Ursino 2000.
dydt(20) = phi_na 	;                % sodium 	   eq 26 Ursino 2000.
dydt(21) = phi_k 		;                % potassium   eq 26 Ursino 2000.
dydt(22) = phi_cl 	;                % chlorine 	 eq 11 Ursino 2000.
dydt(23) = phi_hco3 + Rhco3ic ;      % bicarbonate eq 11 Ursino 2000.
dydt(24) = phi_h 		+ Rhic 	;        % hydrogen 	 eq 11 Ursino 2000.
dydt(25) = phi_p 		+ Rpic 	;        % urea 		   eq 11 Ursino 2000.
%
% for eq 14, Ursino 2000
% State variables 28-34 are Extracellular Solute mass
dydt(26)   = -phi_u 	- J_u/1000 		+ (Qinfused/1000)*cuinf;        % urea 		eq 14, p207, Ursino 2000. combined with eqn 27, Ursino 2008
dydt(27)   = -phi_na 	- J_na/1000       + (Qinfused/1000)*cnainf;   % sodium 	eq 14, p207, Ursino 2000.combined with eqn 27, Ursino 2008
dydt(28)   = -phi_k 	- J_k/1000 		+ (Qinfused/1000)*ckinf;        % potassium  eq 14, p207, Ursino 2000.combined with eqn 27, Ursino 2008
dydt(29)   = -phi_cl 	- J_cl/1000       + (Qinfused/1000)*cclinf;   % chlorine 	eq 14, p207, Ursino 2000.
dydt(30)   = -phi_hco3 	- J_hco3/1000 	+ Rhco3ex;                  % bicarbonate eq 14, p207, Ursino 2000.
dydt(31)   = -phi_h 	- J_h/1000 		+ Rhex	;                       % hydrogen 	eq 14, p207, Ursino 2000.
dydt(32)   = -phi_p 	- J_p/1000 		+ Rpex	;                       % protein 		eq 14, p207, Ursino 2000.

% new V5 Detailed Kidney
for n = 1:12
    dydt(n + 32) = dydtKidney(n);
end;

% -------------------------end of f_pm30D_dialysis RHS function.-------------------------------------------------------------------------.
end

function [p] = pm30D_baroreflex_effector_FIP2021(p, SNA, PNA)
% This function controlls the action on the effector sites of the
% baroreflex. Effector sites include HR (sympathetic and parasympathetic),
% venous tone (zero pressure filling volume), Arterial tone (resistance),
% and right and left ventricular contractility.
%
% references:
% Lin et al. DOI: 10.1177/0954411912451823
%
% "Autonomic control of cardiac pacemaker activity and atrioventricular
% transmission" Levy and Zeiske 1969
%
% all values and equations are from Lin unless otherwise stated
%

DELTAT = p(86);

% State vectors
x0_deltaHR_s    = p(141);
x0_deltaHR_v    = p(142);
x0_sigma_lv     = p(143);
x0_sigma_rv     = p(144);
x0_sigma_V      = p(145);
x0_sigma_R      = p(146);

% Table 3 baroreflex control parameters
G_k_s0 = p(119);    % units: beats/min/Hz
k_k_s0 = p(120);
T_s = p(121);       % units: s

G_v0 = p(122);      % units: beats/min/Hz 45
k_v0 = p(123);
T_v = p(124);       % units: s

% Myocardial constriction
tau_sigma_lv = p(125);  % units: s
T_e_lv = p(126);        % units: s
G_eff_lv = p(127);      % units: mmHg/ml/Hz

tau_sigma_rv = p(128);  % units: s
T_e_rv = p(129);        % units: s
G_eff_rv = p(130);      % units: mmHg/ml/Hz

% Veins
tau_sigma_V = p(131);   % units: s
T_e_V = p(132);         % units: s
G_eff_V = p(133);       % units: ml/Hz

% Arterioles
tau_sigma_R = p(134);   % units: s
T_e_R = p(135);         % units: s
G_eff_R = p(136);       % units: mmHg/(ml/s)/Hz % EDIT G_eff_R increased from 0.2 to 0.21

tau_s = p(137);         % simplified from equation 23, units: s

k_s = G_k_s0 / k_k_s0;  % eq. 22
G_s = 1 * (1 - exp(-k_s * SNA(1)));   % eq. 21 EDIT gain value changed to 1
% sys3 = c2d(tf(G_s,[tau_s 1]),DELTAT);
% [A3, B3, C3, D3] = tf2ss(sys3.Numerator{1},sys3.Denominator{1});
A3 = exp((-1/tau_s) * DELTAT);
B3 = (1-A3)*tau_s;
C3 = G_s/tau_s;
D3 = 0;
x1_deltaHR_s    = A3 * x0_deltaHR_s + B3 * SNA(T_s/DELTAT);
deltaHR_s       = C3 * x0_deltaHR_s + D3 * SNA(T_s/DELTAT);


tau_v = p(138);         % simplified from eq. 28; units: s

k_v = G_v0 / k_v0;      % eq. 27
G_v = 1 * (1 - exp(-k_v * PNA(1)));   % eq. 26, EDIT gain value changed to 1 for more accurate description of HR used in
% sys4 = c2d(tf(G_v,[tau_v 1]),DELTAT);
% [A4, B4, C4, D4] = tf2ss(sys4.Numerator{1},sys4.Denominator{1});
A4 = exp((-1/tau_v) * DELTAT);
B4 = (1-A4)*tau_v;
C4 = G_v/tau_v;
D4 = 0;
x1_deltaHR_v    = A4 * x0_deltaHR_v + B4 * PNA(T_v/DELTAT);
deltaHR_v       = C4 * x0_deltaHR_v + D4 * PNA(T_v/DELTAT);


% % SNA Effector Sites
% again the paper uses inconsistent notation
T_sigma_lv = T_e_lv;
T_sigma_rv = T_e_rv;
T_sigma_V = T_e_V;
T_sigma_R = T_e_R;

% The following are all described by eq. 29
% Left ventricular contractility
% sys5 = c2d(tf(G_eff_lv,[tau_sigma_lv 1]),DELTAT);
% [A5, B5, C5, D5] = tf2ss(sys5.Numerator{1},sys5.Denominator{1});
A5 = exp((-1/tau_sigma_lv) * DELTAT);
B5 = (1-A5)*tau_sigma_lv;
C5 = G_eff_lv/tau_sigma_lv;
D5 = 0;
x1_sigma_lv    = A5 * x0_sigma_lv + B5 * SNA(T_sigma_lv/DELTAT);
sigma_lv       = C5 * x0_sigma_lv + D5 * SNA(T_sigma_lv/DELTAT);

% Right ventricular contractility
% sys6 = c2d(tf(G_eff_rv,[tau_sigma_rv 1]),DELTAT);
% [A6, B6, C6, D6] = tf2ss(sys6.Numerator{1},sys6.Denominator{1});
A6 = exp((-1/tau_sigma_rv) * DELTAT);
B6 = (1-A6)*tau_sigma_rv;
C6 = G_eff_rv/tau_sigma_rv;
D6 = 0;
x1_sigma_rv    = A6 * x0_sigma_rv + B6 * SNA(T_sigma_rv/DELTAT);
sigma_rv       = C6 * x0_sigma_rv + D6 * SNA(T_sigma_rv/DELTAT);

% Venous tone
% sys7 = c2d(tf(G_eff_V,[tau_sigma_V 1]),DELTAT);
% [A7, B7, C7, D7] = tf2ss(sys7.Numerator{1},sys7.Denominator{1});
A7 = exp((-1/tau_sigma_V) * DELTAT);
B7 = (1-A7)*tau_sigma_V;
C7 = G_eff_V/tau_sigma_V;
D7 = 0;
x1_sigma_V    = A7 * x0_sigma_V + B7 * SNA(T_sigma_V/DELTAT);
sigma_V       = C7 * x0_sigma_V + D7 * SNA(T_sigma_V/DELTAT);

% Arterial resistance
% sys8 = c2d(tf(G_eff_R,[tau_sigma_R 1]),DELTAT);
% [A8, B8, C8, D8] = tf2ss(sys8.Numerator{1},sys8.Denominator{1});
A8 = exp((-1/tau_sigma_R) * DELTAT);
B8 = (1-A8)*tau_sigma_R;
C8 = G_eff_R/tau_sigma_R;
D8 = 0;
x1_sigma_R    = A8 * x0_sigma_R + B8 * SNA(T_sigma_R/DELTAT);
sigma_R       = C8 * x0_sigma_R + D8 * SNA(T_sigma_R/DELTAT);

p(141) = x1_deltaHR_s;
p(142) = x1_deltaHR_v;
p(143) = x1_sigma_lv;
p(144) = x1_sigma_rv;
p(145) = x1_sigma_V;
p(146) = x1_sigma_R;

p(92) = deltaHR_s;
p(93) = deltaHR_v;
p(94) = sigma_lv;
p(95) = sigma_rv;
p(96)   = (p(147) - sigma_V)/DELTAT;
p(147) = sigma_V;
p(97) = sigma_R;

end

function p = pm30D_baroreflex_integration_FIP2021(y, p)
%
%
% 12 October, 2020
% Author: Timothy Hunter
% Lead Developer: Sanjay R. Kharche
% Part of PM3
% This code includes implementation of the baroreflex model from Lin et al.
% 2012.
%
% Ref. DOI: 10.1177/0954411912451823
% Constants are from Lin et al. 2012

MAP = y(13);
Pbco2 = p(90);
Pbo2 = p(91);
DELTAT = p(86);

x0_P_aff    = p(139);
x0_temp1    = p(140);

% Table 3 - Baroreflex control parameters
% Afferent compartment:
tau_aff = p(105); % units: seconds
G_aff   = p(106);

% Central Compartment:
tau_c = p(107); % Same as tau_aff,  units: s

% Efferent compartment:
S_p = p(108);
PNA_max = p(109); % units: Hz
PNA_min = p(110); % units: Hz
S_s = p(111);
SNA_max = p(112); % units: Hz
SNA_min = p(113); % units: Hz

% Table A3
k1 = p(114);
k2 = p(115);
k3 = p(116);
k4 = p(117);
k5 = p(118);

% sys1 = c2d(tf(G_aff,[tau_aff 1]),DELTAT);
% [A1, B1, C1, D1] = tf2ss(G_aff,[tau_aff 1]);
A1 = exp((-1/tau_aff) * DELTAT);
B1 = (1-A1)*tau_aff;
C1 = G_aff/tau_aff;
D1 = 0;
x1_P_aff    = A1 * x0_P_aff + B1 * MAP;
P_aff       = C1 * x0_P_aff + D1 * MAP;

% Central Compartment
% Pbco2 = 30; % TEMPORARY
% Pbo2 = 104;
deltaMAP = 0;
if Pbco2 > 40 && Pbo2 < 104
    deltaMAP = k1 + k2 * Pbco2 + k3 / Pbo2;
elseif Pbco2 <= 40 && Pbo2 < 104
    deltaMAP = k1 + k2 * 40 + k3 / Pbo2;
elseif Pbco2 > 40 && Pbo2 >= 104
    deltaMAP = k1 + k2 * Pbco2 + k3 / 104;
end

P_demand = 90 + deltaMAP/100; % Desrcibed in fig. 2 and on page 793

% temp1 is part of eq. 13
% sys2 = c2d(tf(1,[tau_c 1]),DELTAT); % for reference, the tf is of the form 1/(tau_c*s + 1) with s as the laplace variable
% [A2, B2, C2, D2] = tf2ss(sys2.Numerator{1},sys2.Denominator{1});
A2 = exp((-1/tau_c) * DELTAT);
B2 = (1-A2)*tau_c;
C2 = 1/tau_c;
D2 = 0;
x1_temp1    = A2 * x0_temp1 + B2 * P_demand;
temp1       = C2 * x0_temp1 + D2 * P_demand;

P_error = temp1 - P_aff ;% eq. 13


% Efferent Compartment:
% eq. 31 chemoreceptor operating point
if Pbco2 > 40
    deltaG_SNA = k4 * Pbco2 + k5;
else
    deltaG_SNA = k4 * 40 + k5;
end


k_s = (SNA_max - SNA_min) / (4 * S_s); % eq. 15
SNA = ((SNA_max + SNA_min * exp(P_error/k_s))/(1 + exp(P_error/k_s))) * (1 + deltaG_SNA/100); % eq. 14, units: Hz

S_v = S_p; % paper uses inconsistent notation
k_v = (PNA_max - PNA_min) / (4 * S_v); % eq. 17
PNA = (PNA_max - PNA_min * exp(P_error/k_v)) / (1 + exp(P_error/k_v)); % eq. 16

p(139) = x1_P_aff;
p(140) = x1_temp1;
p(88) = SNA;
p(89) = PNA;
end

function [dydt, Qki, Qko] = pm30D_detailedKidney_FIP2021(y,Rkid1, Rkid2, Ck)
% This function describes the hemodynamics of the
neq = 33;
dydt = zeros(12,1);
% Rkid1 = 4.1;
Rkid1 = 12*Rkid1;
% Rkid2 = 0.3;
Rkid2 = Rkid2 * 12;
% Ck = 15;
Ck = Ck / 12;

Qkri1 	= (y(13) - y(0+neq))	/Rkid1	;
Qkri2	= (y(13) - y(1+neq))	/Rkid1	;
Qkri3	= (y(13) - y(2+neq))	/Rkid1	;
Qkri4	= (y(13) - y(3+neq))	/Rkid1	;
Qkri5	= (y(13) - y(4+neq))	/Rkid1	;
Qkri6	= (y(13) - y(5+neq))	/Rkid1	;
Qkli1 	= (y(13) - y(6+neq))	/Rkid1	;
Qkli2	= (y(13) - y(7+neq))	/Rkid1	;
Qkli3	= (y(13) - y(8+neq))	/Rkid1	;
Qkli4	= (y(13) - y(9+neq))	/Rkid1	;
Qkli5	= (y(13) - y(10+neq))	/Rkid1	;
Qkli6	= (y(13) - y(11+neq))	/Rkid1	;

Qkro1 	= (y(0+neq)   - y(8))	/Rkid2	;
Qkro2 	= (y(1+neq)   - y(8))	/Rkid2	;
Qkro3 	= (y(2+neq)   - y(8))	/Rkid2	;
Qkro4 	= (y(3+neq)   - y(8))	/Rkid2	;
Qkro5	= (y(4+neq)   - y(8))	/Rkid2	;
Qkro6 	= (y(5+neq)   - y(8))	/Rkid2	;
Qklo1 	= (y(6+neq)   - y(8))	/Rkid2	;
Qklo2 	= (y(7+neq)   - y(8))	/Rkid2	;
Qklo3 	= (y(8+neq)   - y(8))	/Rkid2	;
Qklo4 	= (y(9+neq)   - y(8))	/Rkid2	;
Qklo5 	= (y(10+neq)   - y(8))	/Rkid2	;
Qklo6 	= (y(11+neq)   - y(8))	/Rkid2	;

Qki = Qkri1 + Qkri2 + Qkri3 + Qkri4 + Qkri5 + Qkri6 + Qkli1 + Qkli2 + Qkli3 + Qkli4 + Qkli5 + Qkli6;
Qko = Qkro1 + Qkro2 + Qkro3 + Qkro4 + Qkro5 + Qkro6 + Qklo1 + Qklo2 + Qklo3 + Qklo4 + Qklo5 + Qklo6;

% y(1) is Pk, pressure at kidneys.
dydt(1)     = (Qkri1 - Qkro1) / Ck; % y(5) is Pk, Ck is non-zero. eq. 23, appendix 2, p585.
dydt(2)     = (Qkri2 - Qkro2) / Ck; % y(5) is Pk, Ck is non-zero. eq. 23, appendix 2, p585.
dydt(3)     = (Qkri3 - Qkro3) / Ck; % y(5) is Pk, Ck is non-zero. eq. 23, appendix 2, p585.
dydt(4)     = (Qkri4 - Qkro4) / Ck; % y(5) is Pk, Ck is non-zero. eq. 23, appendix 2, p585.
dydt(5)     = (Qkri5 - Qkro5) / Ck; % y(5) is Pk, Ck is non-zero. eq. 23, appendix 2, p585.
dydt(6)     = (Qkri6 - Qkro6) / Ck; % y(5) is Pk, Ck is non-zero. eq. 23, appendix 2, p585.
dydt(7)     = (Qkli1 - Qklo1) / Ck; % y(5) is Pk, Ck is non-zero. eq. 23, appendix 2, p585.
dydt(8)     = (Qkli2 - Qklo2) / Ck; % y(5) is Pk, Ck is non-zero. eq. 23, appendix 2, p585.
dydt(9)     = (Qkli3 - Qklo3) / Ck; % y(5) is Pk, Ck is non-zero. eq. 23, appendix 2, p585.
dydt(10)     = (Qkli4 - Qklo4) / Ck; % y(5) is Pk, Ck is non-zero. eq. 23, appendix 2, p585.
dydt(11)     = (Qkli5 - Qklo5) / Ck; % y(5) is Pk, Ck is non-zero. eq. 23, appendix 2, p585.
dydt(12)     = (Qkli6 - Qklo6) / Ck; % y(5) is Pk, Ck is non-zero. eq. 23, appendix 2, p585.

end

function pNew = modParam(p)
% 17 Sept, 2020
% Author: Timothy Hunter
% Lead Developer: Sanjay R. Kharche
% Part of PM3
%
% This program currently applies a time-gradient to the cna_d parameter of
% the PM3 0D Dialysis (driver V1) model. The dialyzer runs for 4 hours then
% shuts off.
%


% DELTAT = p(end);
t = p(1);
pNew = p;


% This gives the Na profile for the dialysate, and turns off dialyzer after
% 4 hours
if t < 240*60 % 4 hours in seconds
    pNew(7) = 141.75 - t.*(4/(240*60));
else
    pNew(7) = 0;
    pNew(4) = 0;
    pNew(85) = 0;
end

end
