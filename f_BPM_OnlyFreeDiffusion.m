%% I. Define Parameters
clearvars; close all
load('CoefficientTablesPade.mat')

%Information on the bead and medium
a           =   0.5e-6;              %Base radius                    (m)            (Radius of sphere on which we are going to put protrusions)
eta         =   8.9e-4;              %Dynamic viscocity of fluid     (Pas)          (Assumed to be water)
T           =   298;                 %Temperature                    (K)            (Assumed to be room temperature)
rho_medium  =   1000;                %Density of medium;             (kg/m^3)       (Assumed to be water)
rho_part    =   1800;                %Density of particle            (kg/m^3)       (Silica 0.5 micron 2 gr/cm^3, Myone 0.5 micron 1.8 g/cm^3, Myone 1.4  micron 1.6 g/cm^3)
slip        =   20e-9;               %Slip length at the surface (Shifted boundary model) (m)   (PEG brush is around 10 nm in thickness, slip length around 20 nm)) 
                                     %Assumption of slip length is based on 'Boundary flow on end-grafted brushes', E. Charrault et al. (2016)  

%Experimental settings
framerate           =       60;                      %Used framerate in the experiment          (Hz.)
t_e                 =       5e-3;                    %Exposure time                             (s)
t_meas              =       200;                     %Total measurement time                    (s)

%% II. Generate and check the particle
%Set initial position of particle
r                           =       [0 ; 0; 2*a]; %Initial position: Start at two times the radius of the sphere (unbound)
e_or                        =       eye(3,3);
    
Volume          =   (4/3)*pi*a^3;
%% III. Initialize
%Some physical entities
kBT                 =       T*(1.38*10^(-23));          %Thermal energy      
dT                  =       1/framerate;                %Time between each location acquisition
mu_tt               =       1/(6*pi*eta*a);         %Translational mobility coefficient at infinite distance from the wall 
mu_rr               =       1/(8*pi*eta*a^3);       %Rotational mobility coefficient at infinite distance from the wall
mu_rt               =       1/(8*pi*eta*a^2);       %tr coefficient of mobility matrix

%Prelocate matrices in which we will save the positions and state of the particle
N               =   round(t_meas*framerate,0);      %Number of localizations
r_blurred       =   zeros(3,N);                     %Particle motion with motion blur
r_noblur        =   zeros(3,N);                     %Particle motion without motion blur

%Random number generation
M            =   10^5;                            %Generate random numbers each M iterations
X1           =   sqrt(kBT)*normrnd(0,1,[6,M]);    %Generate random numbres for the random force at the beginning of the simulation
X2           =   sqrt(kBT)*normrnd(0,1,[6,M]);    %Generate random numbres for the random force at the beginning of the simulation

%Determine constant part of external force and torque working on the particle
%Read as [Fx;Fy;Fz;Tx;Ty;Tz]
F_g       =   [0; 0; -Volume*(rho_part-rho_medium)*9.81; 0; 0; 0]; 
%In this case, only gravity and buyancy force in the z-direction.

%Integration time
dt              =       0.1*t_e;                %Thus, dt is equal to dt_f at the start, since we start unbound (At least, if you put the particle high enough at the start of the simulation :) )
                                                %If State_change = -1, a bond will be deleted
t               =       0;                      %We start at time 0, t is used to keep track of time
i               =       1;                      %Variable to keep track of the used random numbers
j               =       1;                      %Variable to keep track of the number of localizations    
r_blurr         =       [0;0;0];                %Here we will save the localizations over which wel will blurr
m               =       0;                      %Variable to save the number of localizations we will average in timeo include blur

%% III. Integration
tic()
while j <= N  %For all acquisitions
    %We make use of Fixman midpoint scheme to avoid calculating the
    %stochastic drift term explicitely
    %Source: 'Brownian Dynamics of Confined Rigid Bodies', Delong et al.
    %(2015) 
    
    if i > M     %In case we are out of random numbers, we will generate new random numbers
        X1           =   sqrt(kBT)*normrnd(0,1,[6,M]);    %Generate the first batch of random numbers
        X2           =   sqrt(kBT)*normrnd(0,1,[6,M]);    %Generate the second batch of random numbers
        i            =   1;                               %Reset i = 1, since we start again at the first random number
    end
    
    [mu,Z]  =       mu_wall_Jones(a,r(3)+slip,mu_tt,mu_rr,mu_rt,A,B,C);    %First, determine the mobility of the sphere at the current position
    ht              =       mu*F_g;                     %Determine the deterministic displacement due to the force and torque working on the particle  
    vrand           =       chol((4/dt)*mu,'lower')*X1(:,i);    %Determine the noise term due to thermal fluctations
    dr_p            =       0.5*dt*(ht + vrand);                %Calculate the predicted displacement
    
    %Apply reflective boundary conditions on the predicted position:
    if (r(3) + dr_p(3)) <= a                            %Use reflective boundary condition 
        dr_p(3) = dr_p(3) + 2*abs(a-(r(3)+dr_p(3)));    
    end
    
    %Determine predicted positions:
    r_p             =       r + dr_p(1:3);              %predicted position;     
    mu_p            =       mu_wall_Jones_red(a,r_p(3)+slip,mu_tt,mu_rr,mu_rt,A,B,C);      %Determine mobility matrix at the predicted positio
    
    %If there are no bonds, there is only gravity
    hp              =       mu_p*F_g;                   %Determine the deterministic displacement due to the force and torque working on the particle
    
    %determined corrected displacement: 
    dr       =       dt*(hp + sqrt(1/dt)*mu_p*chol(Z,'lower')*(X1(:,i) + X2(:,i)));            %determine the correct displacement

    %Apply reflective boundary condition
    if (r(3) + dr(3)) <= a
        dr(3) = dr(3) + 2*abs(a-(r(3)+dr(3))); 
    end

    %Update position and orientation (test, to see if there is overlap of
    %the protrusions with the boundary)
    r                   =       r + dr(1:3);                %Position update
    Rot                 =       RotationMatrix(dr(4:6));    %Determine rotation matrix
    e_or                =       Rot*e_or;           %Update protrusion location w.r.t. C.O.V. bead (rotate)         
    
    %Data acquisition, with and without motion blurring
    if t >= j*dT - t_e/2 && t <= j*dT + t_e/2         %During the acquisition time,
        r_blurr =  r_blurr + r;                       %Add up al localizations
        m       =  m + 1;                             %And store the number of localizations
    end

    if floor(t/dT) == j                 %Every dT (acquisition)
        r_noblur(:,j) = r;              %save the position without motion blurring for comparison.
    end

    if t > j*dT + t_e/2                 %At the end of the acquisition t
        r_blurred(:,j) = r_blurr/m;     %Save the localization of the particle; divide by the number of added up localizations
        j         = j+1;                %set j one higher, in which j is the acquisition number
        if mod(j,N/20) == 0
           disp(['We are at ',num2str((j/N)*100),' percent, t = ',num2str(round(t,1)),' (s)'])     %Show progress to user:
        end
        m         = 0;                  %set m again 0, since we start counting again.
        r_blurr   = [0;0;0];            %We should set r_blurr again equal to zero
    end   

    t               =       t + dt;         %Store at which time the next iteration will happen
    i               =       i + 1;          %Store which random number we used
end
toc()

%% V. Display simulation result
t           = (1:N)*dT; %total t          (s)

%plot plot x-y and z-position in seperate graphs in time
figure(1)
subplot(4,1,1)
plot(t,r_blurred(1,:)*10^(6))
hold on
title('x-position in time (Blurred)')
ylabel('x-position (\mum)')

subplot(4,1,2)
plot(t,r_blurred(2,:)*10^(6))
hold on
title('y-position in time (Blurred)')
ylabel('y-position (\mum)')

%Plot z-position of particle in time
subplot(4,1,3)
plot(t,(r_blurred(3,:)-a)*10^(6))
hold on
title('z-position in time (Blurred)')
ylabel('z-position (\mum)')

%plot the trajectory in a 3D plot
figure(2)
plot3(r_noblur(1,:)*10^6,r_noblur(2,:)*10^6,r_noblur(3,:)*10^6,'.')
title('3D motion of particle (No Blur)')
xlabel('x - position (\mum)')
ylabel('y - position (\mum)')
zlabel('z - position (\mum)')

% Determine and plot the diffusion coefficient in time
% Use weighted average method
MeasurementWindow   =   framerate;    %Sliding window length in number of frames
maxdt               =   20;
DtraceWindow        =   zeros(maxdt,numel(r_blurred(1,:)) - MeasurementWindow + 1);
weight              =   zeros(1,maxdt);

for dt = 1:maxdt
        SDtrace             =   (r_blurred(1,1+dt:end) - r_blurred(1,1:end-dt)).^2 + (r_blurred(2,1+dt:end)-r_blurred(2,1:end-dt)).^2;
        SDtraceWindow       =   movmean(SDtrace,MeasurementWindow-dt,'Endpoints','discard'); %mean squared displacement
        DtraceWindow(dt,:)  =   SDtraceWindow*framerate/(4*dt);
        Vrel                =   dt*(2*dt^2+1)/(MeasurementWindow-dt+1);      %Calculate relative variance
        weight(dt)          =   1/Vrel; %assign weights to data points (rows) DtraceWindow)
end
 
sumweight = sum(weight);
weight = weight/sumweight;
weight = weight';
Dcoef_blurred = sum(DtraceWindow.*weight); %Calculate diffusion coefficient
 
figure(3)
hold on
plot(r_blurred(1,:)*10^6,r_blurred(2,:)*10^6,'k .')
axis equal
title('Motion (blurred) in x-y plane (k = 0, r = 1, g = 2, b = > 2)')
xlabel('x-position (\mum)')
ylabel('y-position (\mum)')
  

% Plot xy position in time (blurred as measured in the experimental set-up)
figure(4)
hold on
plot(t,r_blurred(1,:)*10^6)
plot(t,r_blurred(2,:)*10^6)
xlim([0 max(t)])
ylabel('Position (\mum)')
xlabel('Time (s)')
legend('x-position','y-position')
title('Motion in time (blurred)')

figure(5)
hold on
plot(t(1:size(Dcoef_blurred,2)),10^(12)*Dcoef_blurred);
ylabel('Diffusion Coefficient (\mum^2/s)')
title('Diffusion Coefficient in time')
xlabel('t (s)')
xlim([0 max(t)])
ylim([0 1.1*max(10^(12)*Dcoef_blurred)])
yyaxis right
ylabel('State')
ylim([-2.1 0.1])
yticks([-2 -1 0])

figure(6)
hold on
histogram(Dcoef_blurred*10^(12),75)
xlim([0, 1.1*max(Dcoef_blurred*10^(12))])
xlabel('Diffusion Coefficien(\mum^2/s)')
ylabel('Counts (-)')
title('Measured Diffusion Coefficient')
