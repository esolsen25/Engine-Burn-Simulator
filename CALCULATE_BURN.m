function [] = CALCULATE_BURN()
%% INPUT GLOBAL DECLARATION
global dt AMB_P EXIT_area C_star_EFF NOZ_EFF dP OX_P OX_mass FUEL_mass ...
    OX_density OX_D_orifice OX_N_orifice OX_OD_annulus OX_ID_annulus OX_Discharge_Coeff ...
    FUEL_density FUEL_D_orifice1 FUEL_N_orifice1 FUEL_D_orifice2 FUEL_N_orifice2 ...
    FUEL_OD_annulus FUEL_ID_annulus FUEL_Discharge_Coeff EXPANSION_RATIO THROAT_area
%% OUTPUT GLOBAL DECLARATION
global OVERALL_EFF FUEL_P CHAMBER_P CHAMBER_temp EXIT_temp EXIT_vel EXIT_P...
    EXIT_mach gam0 MW R CdA_OX CdA_FUEL MR OX_flowrate FUEL_flowrate NET_flowrate ...
    OXGEN_mass FGEN_mass OX_Status Status t CHAMBER_P_PSI DELTA_P DELTA_P_PSI ...
    EXIT_P_PSI impulse DELIVERED_impulse thrust thrust_lbs
%% VARIABLES - BASED ON INPUTS
OVERALL_EFF=C_star_EFF*NOZ_EFF; %[-] - Overall Efficiency
FUEL_P(1)=OX_P(1)-dP/(6894.76*1000);%[kPa] - Fuel Initial Pressure
%% CONSTANTS
t(1)=0;%[s] - Initial Time
CHAMBER_P(1)=0;%[kPa] - Chamber Initial Pressure
CHAMBER_temp=[500 1000];%[K] - Chamber Initial Temperatures
OX_Status=["Liquid"]; % - Initial Oxidizer State
EXIT_mach=[0.5,1.0]; %[-] - Initial Mach Number at Exit
gam0(1:2)=1.2; % - Gamma
MW(1:2)=28.9647; % - Initial (1-2) Molar Weight
R=8314./MW; % - Initial (1-2) Gas Constant
%% OXIDIZER INJECTOR AREAS AND CdA
% Calculations - Orifice and Annulus Areas, CdA_OX
OX_A_orifice = (pi*(OX_D_orifice/2)^2)/1550.003; % [m^2]
OX_A_annulus = (pi*(OX_OD_annulus/2)^2-pi*(OX_ID_annulus/2)^2)/1550.003; % [m^2]

CdA_OX = (OX_A_orifice*OX_N_orifice+OX_A_annulus)*OX_Discharge_Coeff; % [m^2]
%% FUEL INJECTOR AREAS AND CdA
% Calculations - Orifice, Annulus Areas, and CdA_FUEL
FUEL_A_orifice1 = (pi*(FUEL_D_orifice1/2)^2)/1550.003; % [m^2]
FUEL_A_orifice2 = (pi*(FUEL_D_orifice2/2)^2)/1550.003; % [m^2]
FUEL_A_annulus = (pi*(FUEL_OD_annulus/2)^2-pi*(FUEL_ID_annulus/2)^2)/1550.003; % [m^2]

CdA_FUEL = (FUEL_A_orifice1*FUEL_N_orifice1+FUEL_A_orifice2*FUEL_N_orifice2+FUEL_A_annulus)*FUEL_Discharge_Coeff; % [m^2]
%% INITIAL VALUES
i=1;
% Calculate Initial Flowrates based on CHAMBER_P(1)=0[kPa]
OX_flowrate(i) = CdA_OX*sqrt(2*OX_density*(OX_P(i)-CHAMBER_P(i))*1000); % [kg/s]
FUEL_flowrate(i) = CdA_FUEL*sqrt(2*FUEL_density*(FUEL_P(i)-CHAMBER_P(i))*1000); % [kg/s]
NET_flowrate(i)=OX_flowrate(i)+FUEL_flowrate(i);

MR(i)=OX_flowrate(i)/FUEL_flowrate(i);

% Calculate Initial OXGEN_mass based on OX_flowrate(1)
OXGEN_mass(i)=OX_flowrate(i)*dt;
FGEN_mass(i)=FUEL_flowrate(i)*dt;

% Initalize Loop Control
i=2;
while(OXGEN_mass(i-1)>0 && FGEN_mass(i-1)>0)
    %% TIME
    t(i)=t(i-1)+dt;
    %% gam0, MW, R
    if(i>2)
        % Calculate Constant Values from Lookup
        gam0(i)=gam0_calc(CHAMBER_P(i-1),MR(i-1));
        MW(i)=MW_calc(CHAMBER_P(i-1),MR(i-1));
        R(i)=8314/MW(i);
    end
    %% OX_status
    if((OX_mass(i-1)/OX_mass(1))>0.1)
        OX_Status(i)="Liquid";
    else
        OX_Status(i)="Gas";
    end

    if(OX_mass(i-1)>0 && FUEL_mass(i-1)>0)
        CHAMBER_P(i)=NET_flowrate(i-1)*sqrt(CHAMBER_temp(i-1))/(THROAT_area*sqrt(gam0(i)/R(i))*((gam0(i)+1)/2)^(-(gam0(i)+1)/(2*(gam0(i)-1))))/1000;
    end
    %% CHAMBER_temp (T0)
    if(i>2)
        T0=T0_calc(MR(i-1),CHAMBER_P(i-1));%[K] - Chamber Temperature
        % Calculates by which factor to multiply the result, depending on if
        % the previous pressure is greater than the current, or not.
        if(CHAMBER_P(i)<CHAMBER_P(i-1))
            factor = CHAMBER_P(i)/CHAMBER_P(i-1);
        else
            factor = 1;
        end
        CHAMBER_temp(i)=C_star_EFF*T0*factor;
    end
    %% OX_P
    if(strcmp(OX_Status(i-1),"Liquid"))
        if(OX_mass(i-1)>0 && FUEL_mass(i-1)>0)
            OX_P(i)=OX_P(1)*(OX_mass(i-1)/OX_mass(1)*0.3+0.7); % [kPa]
        end
    elseif(strcmp(OX_Status(i-1),"Gas"))
        OX_P(i)=OX_mass(i-1)/exp(0.05*(t(i)-t(i-1))); % [kPa]
    end
    %% FUEL_P
    if(OX_mass(i-1)>0 && FUEL_mass(i-1)>0)
        if(OX_P(i)-dP/6.895<0)
            FUEL_P(i)=0;%[kPa]
        else
            FUEL_P(i)=OX_P(i)-dP/6.895; % [kPa]
        end
    end
    %% OX,FUEL,NET_flowrate
    % Calculate Step Flowrates based on CHAMBER_P(i)
    if(FUEL_P(i)>0 && OX_P(i)>0)
        OX_flowrate(i) = CdA_OX*sqrt(2*OX_density*(OX_P(i)-CHAMBER_P(i))*1000); % [kg/s]
        FUEL_flowrate(i) = CdA_FUEL*sqrt(2*FUEL_density*(FUEL_P(i)-CHAMBER_P(i))*1000); % [kg/s]
    else
        OX_flowrate(i)=0;
        FUEL_flowrate(i)=0;
    end
    NET_flowrate(i)=OX_flowrate(i)+FUEL_flowrate(i);
    %% MR
    % Calculate Step MR
    if(OX_flowrate(i)==0 || FUEL_flowrate(i)==0)
        if(MR(i-1) ~= 0)
            MR(i)=0; %[-]
        end
    else
        MR(i)=OX_flowrate(i)/FUEL_flowrate(i); %[-]
    end
    %% OX_mass, OXGEN_mass
    % Calculate Step OX_mass and OXGEN_mass
    if((OX_mass(i-1)-OXGEN_mass(i-1))<0)
        OX_mass(i)=0; %[kg]
    else
        OX_mass(i)=OX_mass(i-1)-OXGEN_mass(i-1); %[kg]
    end
    if(OX_mass(i)==0)
        OXGEN_mass(i)=0; %[kg]
    else
        OXGEN_mass(i)=OX_flowrate(i)*dt;%[kg]
    end
    %% FUEL_mass, FGEN_mass
    % Calculate Step FUEL_mass and FGEN_mass
    if((FUEL_mass(i-1)-FGEN_mass(i-1))<0)
        FUEL_mass(i)=0;%[kg]
    else
        FUEL_mass(i)=FUEL_mass(i-1)-FGEN_mass(i-1);%[kg]
    end
    if(FUEL_mass(i)==0)
        FGEN_mass(i)=0;
    else
        FGEN_mass(i)=FUEL_flowrate(i)*dt;
    end
    i=i+1;
end
%% OUTSIDE LOOP CONTROL CALCULATIONS
CHAMBER_P(end)=AMB_P/1000;%[kPa]
CHAMBER_P_PSI=CHAMBER_P/6.895;%[psi]

for i=1:length(CHAMBER_P)
    if(OX_P(i)-CHAMBER_P(i)>0)
        DELTA_P(i)=OX_P(i)-CHAMBER_P(i);%[kPa]
    else
        DELTA_P(i)=0;
    end
end

DELTA_P_PSI=DELTA_P/6.895;%[psi]

% Calculate exit mach matrix using gam0 values
EXIT_mach(3)=OVERALL_EFF*EXIT_mach_calc(EXPANSION_RATIO,gam0(3));
for i=4:length(gam0)
    EXIT_mach(i)=EXIT_mach_calc(EXPANSION_RATIO,gam0(i));%[-]
end

% Calculate Exit Pressure
EXIT_P=[0,CHAMBER_P(2)/10.0];%[kPa]
for i=3:length(CHAMBER_P)
    EXIT_P(i) = CHAMBER_P(i)*(1+(gam0(i)-1)/2*EXIT_mach(i)^2)^(-gam0(i)/(gam0(i)-1));%[kPa]
end
EXIT_P_PSI=EXIT_P/6.895;%[psi]

% Calculate Exit Temperature
for i=1:length(CHAMBER_temp)
    EXIT_temp(i)=CHAMBER_temp(i)/(1+(gam0(i)-1)/2*EXIT_mach(i)^2);%[K]
end

% Calculate Exit Velocity
for i=1:length(EXIT_mach)
    EXIT_vel(i)=EXIT_mach(i)*sqrt(gam0(i)*R(i)*EXIT_temp(i));%[m/s]
end

% Calculate Thrust
thrust=[NET_flowrate(1)*EXIT_vel(1),NET_flowrate(2)*EXIT_vel(2)];
for i=3:length(EXIT_vel)
    if(NET_flowrate(i)*EXIT_vel(i)+(EXIT_P(i)*1000-AMB_P)*EXIT_area > 0)
        thrust(i)=NET_flowrate(i)*EXIT_vel(i)+(EXIT_P(i)*1000-AMB_P)*EXIT_area;%[N]
    else
        thrust(i)=0;
    end
end
thrust_lbs=thrust*0.2248090795;

% Calculate Impulse
impulse=thrust.*dt;

% Calculate Delivered Impulse
DELIVERED_impulse(1)=0;
for i=2:length(impulse)
    DELIVERED_impulse(i)=DELIVERED_impulse(i-1)+impulse(i-1);
end

% Calculate Fuel/OX Status
Status=["Burning"];
for i=2:length(FUEL_mass)
    if(strcmp(Status(i-1),"FD") || strcmp(Status(i-1),"OD"))
        Status(i)=Status(i-1);
    else
        if(FUEL_mass(i)==0)
            Status(i)="FD";
        else
            if(OX_mass(i)==0)
                Status(i)="OD";
            else
                Status(i)="Burning";
            end
        end
    end
end
end
