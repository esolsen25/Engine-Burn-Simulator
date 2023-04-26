%% INTRODUCTION
% Title: ISOPROPANOL & NITROUS OXIDE BURN SIMULATOR
% Developer: Evan Olsen
% Date: 5-3-2023
% Project Iteration: Nitron I
% Team: RIT Liquid Propulsion
clc; clear; close all;
%% INPUT GLOBAL DECLARATION
global dt AMB_P EXIT_area C_star_EFF NOZ_EFF dP OX_P OX_mass FUEL_mass ...
 OX_density OX_D_orifice OX_N_orifice OX_OD_annulus OX_ID_annulus OX_Discharge_Coeff ...
 FUEL_density FUEL_D_orifice1 FUEL_N_orifice1 FUEL_D_orifice2 FUEL_N_orifice2 ...
 FUEL_OD_annulus FUEL_ID_annulus FUEL_Discharge_Coeff EXPANSION_RATIO THROAT_area

OX_P = []; OX_mass = []; FUEL_mass = []; 
%% OUTPUT GLOBAL DECLARATION
global OVERALL_EFF FUEL_P CHAMBER_P CHAMBER_temp EXIT_temp EXIT_vel EXIT_P...
EXIT_mach gam0 MW R CdA_OX CdA_FUEL MR OX_flowrate FUEL_flowrate NET_flowrate ...
OXGEN_mass FGEN_mass OX_Status Status t CHAMBER_P_PSI DELTA_P DELTA_P_PSI ...
EXIT_P_PSI impulse DELIVERED_impulse thrust thrust_lbs

FUEL_P = []; CHAMBER_P = []; CHAMBER_temp = []; EXIT_temp = []; EXIT_vel = [];
EXIT_P = []; EXIT_mach = []; gam0 = []; MW = []; R = []; MR = []; OX_flowrate = [];
FUEL_flowrate = []; NET_flowrate = []; OXGEN_mass = []; FGEN_mass = []; OX_Status = [];
t = []; CHAMBER_P_PSI = []; DELTA_P = []; DELTA_P_PSI = []; EXIT_P_PSI = []; impulse = [];
DELIVERED_impulse = []; thrust = []; thrust_lbs = [];
%% SETUP
dt=0.025;%[s] - Timestep
AMB_P=101300;%[Pa] - Ambient Pressure
g=9.81;%[m/s^2]
%% NOZZLE
THROAT_dia=0.5;%[in]
    THROAT_area=(pi/4)*THROAT_dia^2/1550;%[m^2]
EXIT_dia=1.1;%[in]
    EXIT_area=(pi/4)*EXIT_dia^2/1550;%[m^2]
EXPANSION_RATIO=EXIT_area/THROAT_area;%[-] - Expansion Ratio
NOZ_EFF=0.98; %[-] - Nozzle Efficiency
%% COMBUSTION CHAMBER
CHAMBER_dia=2;%[in]
CHAMBER_len=4;%[in]
C_star_EFF=0.80; %[-] - C-star Efficiency
%% FUEL SPECIFICATIONS
FUEL_type = "Isopropanol";

FUEL_TANK_OD=1.75;%[in]
FUEL_TANK_ID=1.610;%[in]
FUEL_TANK_len=25;%[in]
FUEL_displacement=20;%[in]
FUEL_TANK_vol=pi*FUEL_TANK_OD^2*FUEL_TANK_len/2;%[in^3]

dP=15;%[psi] - Pressure Difference
FUEL_density=786;%[kg/m^3]

FUEL_vol=pi*FUEL_TANK_ID^2*FUEL_displacement/4;%[in^3]
FUEL_mass(1)=FUEL_density*FUEL_vol/61024;%[kg]
%% OXIDIZER SPECIFICATIONS
OX_type="Nitrous Oxide [N20]";

OX_TANK_OD=3;%[in]
OX_TANK_ID=2.750;%[in]
OX_displacement=26;%[in]
TANK_style="Concentric";
SIPHON_dia=0;%[in]

% could use a lookup table to determine the density at a given temperature
OX_temp=25;%[C]
OX_density=744;%[kg/m^3]
OX_P(1)=5660;%[kPa]
OX_vol=pi*OX_TANK_ID^2*OX_displacement/4-pi*FUEL_TANK_OD^2*FUEL_TANK_len/4-...
    pi*SIPHON_dia^2*(OX_displacement-FUEL_TANK_len)/4;%[in^3]

if(SIPHON_dia==0)
    if(FUEL_TANK_len<OX_displacement/2)
        OX_usable_vol=OX_vol-pi*FUEL_TANK_OD^2*FUEL_displacement/4;%[in^3]
    else
        OX_usable_vol=pi*(OX_TANK_ID^2-FUEL_TANK_OD^2)*OX_displacement/4;%[in^3]
    end
else
    OX_usable_vol=OX_vol;%[in^3]
end

OX_total_mass=OX_vol*OX_density/61024;%[kg]
OX_usable_mass=OX_usable_vol*OX_density/61024;%[kg]
OX_mass(1)=OX_usable_mass;%[kg]
%% OXIDIZER INJECTOR PARAMETERS
OX_D_orifice = 0.067; % [in]
OX_N_orifice = 6; % [-]
OX_OD_annulus = 0; % [in]
OX_ID_annulus = 0; % [in]
OX_Discharge_Coeff = 0.25; % [-]
%% FUEL INJECTOR PARAMETERS
FUEL_D_orifice1 = 0.0465; % [in]
FUEL_N_orifice1 = 6; % [-]
FUEL_D_orifice2 = 0; % [in]
FUEL_N_orifice2 = 0; % [in]
FUEL_OD_annulus = 0; % [in]
FUEL_ID_annulus = 0; % [in]
FUEL_Discharge_Coeff = 0.25; % [-]
%% CALCULATE
CALCULATE_BURN();
%% RESULTS
CHAMBER_P_AVG=mean(CHAMBER_P(7:14));%[kPa]
CHAMBER_P_PSI_AVG=CHAMBER_P_AVG/6.895;%[psi]

THRUST_AVG=mean(thrust(7:14));%[N]
THRUST_LBS_AVG=THRUST_AVG*0.2248090795;%[lbf]

OX_FLOWRATE_AVG=mean(OX_flowrate(7:14));%[kg/s]
FUEL_FLOWRATE_AVG=mean(FUEL_flowrate(7:14));%[kg/s]
FLOWRATE_AVG=mean(NET_flowrate(7:14));%[kg/s]

MIXTURE_RATIO_AVG=mean(MR(7:14));%[-]

SPECIFIC_impulse=THRUST_AVG/(FLOWRATE_AVG*g);
TOTAL_impulse=sum(impulse,"all");

CHAR_len=(((pi*CHAMBER_dia^2)/4*CHAMBER_len)/61024)/THROAT_area;

DELTA_P_AVG=mean(DELTA_P(7:14));%[kPa]
DELTA_P_PSI_AVG=mean(DELTA_P_PSI(7:14));%[psi]
PRESSURE_drop=DELTA_P./CHAMBER_P;%[-]
PRESSURE_drop_AVG=mean(PRESSURE_drop(7:14));%[-]

count=0;
for i=1:length(Status)
    if(Status(i)=="Burning")
        count=count+1;
    end
end
burnTime=count*dt;%[s]
%% COMMAND WINDOW OUTPUT
fprintf('CHAMBER PRESSURE:   %.2f  [psi]\n\n',CHAMBER_P_PSI_AVG)
fprintf('INITIAL THRUST:     %.2f  [lbs]\n                    %.2f [N]\n\nMASS FLOWRATE:      %.3f   [kg/s]\n\n',THRUST_LBS_AVG,THRUST_AVG,FLOWRATE_AVG);
fprintf('MIXTURE RATIO(O/F): %.3f   [-]\n\nSPECIFIC IMPULSE:   %.2f  [s]\n\nPRESSURE Δ:         %.2f  [psi]\n                    %.2f [kPa]\n\n',MIXTURE_RATIO_AVG,SPECIFIC_impulse,DELTA_P_PSI_AVG,DELTA_P_AVG);
fprintf('BURN TIME (FUEL):   %.2f    [s]\n\nTOTAL IMPULSE:      %.2f [Ns]\n',burnTime,TOTAL_impulse);
%% FIGURE WINDOW OUTPUT
figure('Name','Performance Plots vs. Time');
set(gcf,'color','k');
subplot(2,2,1);
plot(t,thrust_lbs,'color','red','LineWidth',1.5);
title('\color{white}THRUST');
xlabel('\color{white}Time [s]');ylabel('\color{white}Thrust [lbs]');
grid on
set(gca,'Color','k');
set(gca,'YColor','w');set(gca,'XColor','w');
subplot(2,2,2);
plot(t,CHAMBER_P_PSI,'color','yellow','LineWidth',1.5);
title('\color{white}CHAMBER_P');
xlabel('\color{white}Time [s]');ylabel('\color{white}Pressure [psi]');
grid on
set(gca,'Color','k');
set(gca,'YColor','w');set(gca,'XColor','w');
subplot(2,2,3);
plot(t,NET_flowrate,'color','green','LineWidth',1.5);
title('\color{white}ṁ_T_O_T_A_L')
xlabel('\color{white}Time [s]');ylabel('\color{white}Mass Flowrate [kg/s]');
grid on
set(gca,'Color','k');
set(gca,'YColor','w');set(gca,'XColor','w');
subplot(2,2,4)
plot(t,DELTA_P_PSI,'color','blue','LineWidth',1.5);
title('\color{white}Δ_P [CHAMBER_P-OX_P]');
xlabel('\color{white}Time [s]');ylabel('\color{white}Pressure [psi]');
grid on
set(gca,'Color','k');
set(gca,'YColor','w');set(gca,'XColor','w');