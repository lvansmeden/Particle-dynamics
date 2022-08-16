%% I. Define Parameters
clearvars; close all
load('CoefficientTablesPade.mat')

%Information on the bead and medium
a_base      =   0.5e-6;              %Base radius                    (m)            (Radius of sphere on which we are going to put protrusions)
eta         =   8.9e-4;              %Dynamic viscocity of fluid     (Pas)          (Assumed to be water)
T           =   298;                 %Temperature                    (K)            (Assumed to be room temperature)
rho_medium  =   1000;                %Density of medium;             (kg/m^3)       (Assumed to be water)
rho_part    =   1800;                %Density of particle            (kg/m^3)       (Silica 0.5 micron 2 gr/cm^3, Myone 0.5 micron 1.8 g/cm^3, Myone 1.4  micron 1.6 g/cm^3)
slip        =   20e-9;               %Slip length at the surface (Shifted boundary model) (m)   (PEG brush is around 10 nm in thickness, slip length around 20 nm)) 
                                     %Assumption of slip length is based on 'Boundary flow on end-grafted brushes', E. Charrault et al. (2016)  
                                     
%Set the protrusions:
n_protrusions       =    100;              %Max Number of protrusions on the particle       (If this number does not fit without overlap, the algorithm will use less)
d_min               =    20e-9;            %Minimum radius of a protrusion  (m)
d_max               =    60e-9;            %Maximum radius of a protrusion  (m)
offset_max          =    .3;              %Number between 0 and 1, which indicates by how much a protrusion can stick out of or can be embedded in the bead
                                  
l_0                 =    100e-9;           %Distance between active binders on the surface     (m) 
%Note that altough we only plot a limited number of binders on the surface,
%that in principle we have an 'infinite' number of binders on the surface.

Rc                  =    25e-9;        %Cut-off interaction radius (nm)
%This is the maximum distance between binders between which a bond can be formed 
%note, it would be unphysical to let this distance be much longer than the
%contour length of your bond (0.34 nm/bp for ds-DNA)
%(NOTE: When Rc >= l_0, this model is no longer valid and will give
%errorenous results)

%Experimental settings
framerate           =       60;                      %Used framerate in the experiment          (Hz.)
t_e                 =       5e-3;                    %Exposure time                             (s)
t_meas              =       60;                     %Total measurement time                    (s)

%Information on the bonds and bond formation
k_bond              =       10e-5;               %Bond stifness                 (N/m) (Currently set for ds-DNA with a contour length of 50 bp)
%You can look up which value of k to use via the lookuptable here: \\physstor\mbx-common\02 People and Projects\02 MSc\Rik van Haaften

%When simulating ds-DNA, You can look up which value of k to use via the lookuptable here:
%\\physstor\mbx-common\02 People and Projects\02 MSc\Rik van Haaften
%Tables are based on measuring the same Diffusion coefficient OR motion amplitude OR max
%radial discplacement as a simplified WLC model.

%Determines probability to bind and/or unbind when patches are close to one another
k_off          =    1/8;               %Dissocation rate   (1/s) 
k_on_1         =    1/8;               %Association rate: Free -> Single bound  (1/s)
k_on_2         =    1/8;               %Association rate: Single bound -> Double bound (1/s)

Prob_S_D       =    k_on_2/(k_on_2 + k_off);    %Probability to go from single bound to double bound
%Note that no real chemical interactions are modelled in this code, yet
%that association times and dissociation times are generated after each
%event after which the code will wait to let the next event happen.

%Association is after such a waiting governed by 'hit-and-stick'

%For the association, this might mean that there is a difference between
%the set k_on and the output k_on, since the particle might be too far away
%from the surface when an binding event is due.

%Especially going from a freely moving state to a bound state, this method
%might give huge errors in association times.
%Yet, since we simulate free diffusion at a larger time step than the actual chemical
%reaction, we can not include real adhesion.

%% II. Generate and check the particle
%Set initial position of particle
r                           =       [0 ; 0; 2*a_base]; %Initial position: Start at two times the radius of the sphere (unbound)

%Generate protrusions and binders on the surface of the particle:
[d, e_protrusion, Volume]   =   GenerateParticle_nobinders(d_min, d_max, n_protrusions, a_base, offset_max);
n_protrusions               =   size(e_protrusion,2);

%You can put in and edit the code to have any type of configuration
a_eff                       =   (3*Volume/(4*pi))^(1/3);  %Effective radius (Bead + protrusions) of the bead
%Volume of the total particle is determined and can be used to determine
%some kind of effective radius

%Warn user for possible errors and possible fixes:
if a_base + slip < a_eff
    disp('----------------------')
    disp('These settings might give an error')
    disp('If that happens, resolve by doing one of the following:')
    disp('-  In case of requiring only small number of protrusions, use smaller protrusions')
    disp('-  Or, use more protrusions and set d_min to a smaller value (larger coverage of the sphere)')
    disp('-  Or, use a_eff at the reflective boundary (Unphysical at a small coverage of protrusions)')
    disp('-  Or, increase the slip length (Large values of slip are unphysical as well)')
    disp('----------------------')
end

%Plot the particle you are going to use to check
close all
[X,Y,Z] = sphere;   %Generate a sphere
figure(1)
surf(10^(6)*a_base*X,10^(6)*a_base*Y,10^(6)*(a_base*1.5  + a_base*Z),'facecolor',[0.5 0.5 0.5]) %Plot all the protrusions
for i = 1:size(e_protrusion,2)
figure(1)
hold on
surf(10^(6)*(d(i)*X + e_protrusion(1,i)), 10^(6)*(d(i)*Y + e_protrusion(2,i)), 10^(6)*(d(i)*Z + a_base*1.5 + e_protrusion(3,i)),'facecolor',[0.5 0.5 0.5])
end
    
figure(1)
xlabel('\mum')
ylabel('\mum')
zlabel('\mum')
title('Is this particle okay? (yes = [space], no = [ctrl] + [c])')
hold on
surf(linspace(-a_base*10^6,a_base*10^6,2),linspace(-a_base*10^6,a_base*10^6,2),zeros(2),'facecolor',[0.1 0.3 0.4])
surf(linspace(-a_base*10^6,a_base*10^6,2),linspace(-a_base*10^6,a_base*10^6,2),zeros(2)-(10^6)*slip,'facecolor',[0.1 0.2 0.1])
scatter3((-2*l_0:l_0:2*l_0)*10^6,-2*l_0*ones(1,5)*10^6,zeros(1,5),'r','filled')
scatter3((-2*l_0:l_0:2*l_0)*10^6,-1*l_0*ones(1,5)*10^6,zeros(1,5),'r','filled')
scatter3((-2*l_0:l_0:2*l_0)*10^6,zeros(1,5),zeros(1,5),'r','filled')
scatter3((-2*l_0:l_0:2*l_0)*10^6,1*l_0*ones(1,5)*10^6,zeros(1,5),'r','filled')
scatter3((-2*l_0:l_0:2*l_0)*10^6,2*l_0*ones(1,5)*10^6,zeros(1,5),'r','filled')
pause()

%% III. Initialize
%Some physical entities
kBT                 =       T*(1.38*10^(-23));          %Thermal energy      
dT                  =       1/framerate;                %Time between each location acquisition
mu_tt               =       1/(6*pi*eta*a_eff);         %Translational mobility coefficient at infinite distance from the wall 
mu_rr               =       1/(8*pi*eta*a_eff^3);       %Rotational mobility coefficient at infinite distance from the wall
mu_rt               =       1/(8*pi*eta*a_eff^2);       %tr coefficient of mobility matrix

%Prelocate matrices in which we will save the positions and state of the particle
N               =   round(t_meas*framerate,0);      %Number of localizations
r_blurred       =   zeros(3,N);                     %Particle motion with motion blur
r_noblur        =   zeros(3,N);                     %Particle motion without motion blur
Bonds           =   zeros(1,N);                     %Number of bonds at the moment of particle location acquisition

%Random number generation
M            =   10^5;                            %Generate random numbers each M iterations
X1           =   sqrt(kBT)*normrnd(0,1,[6,M]);    %Generate random numbres for the random force at the beginning of the simulation
X2           =   sqrt(kBT)*normrnd(0,1,[6,M]);    %Generate random numbres for the random force at the beginning of the simulation

%Determine constant part of external force and torque working on the particle
%Read as [Fx;Fy;Fz;Tx;Ty;Tz]
F_g       =   [0; 0; -Volume*(rho_part-rho_medium)*9.81; 0; 0; 0]; 
%In this case, only gravity and buyancy force in the z-direction.

%Integration time
dt_f            =       0.1*t_e;                %Simulation integration time while unbound (s)
dt              =       dt_f;                   %Thus, dt is equal to dt_f at the start, since we start unbound (At least, if you put the particle high enough at the start of the simulation :) )
dt_b1           =       min([dt_f 1/(mu_tt*k_bond)]);   %Integration time while having one bond (s)
dt_b2           =       min([dt_f 0.5/(mu_tt*k_bond)]); %Integration time while having two bonds (s)

%Create objects, used throuhgout the simulation
anch_point      =       zeros(3,1);             %In this matrix, we will store the positions of ligands on the surface to which the particle has bound
e_binder        =       zeros(3,1);             %Matrix, in which we are going to store binder orientation vectors when a bond is present
State           =       0;                      %We start in the freely moving state (State = 0)
State_change    =       1;                      %State_change indicates whether a bond will be formed or deleted. Throughout a simulation
                                                    %If State_change = 1, a bond will be added 
                                                    %If State_change = -1, a bond will be deleted
t               =       0;                      %We start at time 0, t is used to keep track of time
i               =       1;                      %Variable to keep track of the used random numbers
j               =       1;                      %Variable to keep track of the number of localizations    
r_blurr         =       [0;0;0];                %Here we will save the localizations over which wel will blurr
m               =       0;                      %Variable to save the number of localizations we will average in timeo include blur

failure         =       0;                      %Use failure to keep track of the number of failed iterations (Exclusion of the protrusions)
                                                %If too many successive failures occur (more than 10000 after another) a bond will be disregarded
t_wait          =       t + GenerateWaitTime(k_on_1); %Time, at which we wait at minimum before the next state
                                                %We will update t_wait after each event
Q               =       1;                      %Q is used (among others) to keep track of whether to calculate the forces again at the start of each iteration
deltaT          =       0;                      %Use deltaT to keep track of the errors in the association time;
AssociationError=       zeros(1,round(t_meas*(k_on_1+k_off),0)); %Generate a zero array to save the errors in association time in
WW              =       1;                      %Use WW to keep track of the number of errors in the association time;

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
    
    if Q == 1 %When we are at a new position (the previous iteration was successfull), it makes sence to calculate the forces working on the particle again
        [mu,Z]  =       mu_wall_Jones(a_eff,r(3)+slip,mu_tt,mu_rr,mu_rt,A,B,C);    %First, determine the mobility of the sphere at the current position

        if State > 0            %If there are bonds, calculate corresponding forces and torques
            %Direct influence of the bond:
            f_bond          =       -k_bond*(r + e_binder - anch_point);    %First, determine all forces working on the particle    
            F_bond          =       [sum(f_bond,2); sum(cross(e_binder,f_bond),2)];  %Determine the torque and add up al forces and torques due to the bonds

            %Determine the volume exclusion effect: governed by Eq. 6 by the
            %paper of Segall et al. (2006) "Volume-Exclusion Effects in
            %Tethered-Particle Experiments: Bead Size Matters"
            f_ve            =       [0;0;1].*(kBT./(r(3) + e_binder(3,:)));
            
            %check whether particle is sufficiently close (Untick these
            %lines in case of bonds which are two times larger than the
            %radius of the sphere; this is not the case for f-BPM most of
            %the times).
%             for Q = 1:State                                         %Check all the bonds 
%                 if r(3) + e_binder(3,bond_index(Q)) > 2*a_eff       %If distance between particle and surface is large, there is no longer a effect
%                     f_ve(:,Q) = [0;0;0];                            %So set the corresponding force equal to zero
%                 end
%             end            
            
            F_ve            =       [sum(f_ve,2); sum(cross(e_binder,f_ve),2)]; %Determine the torque and add up al forces and torques due to the bonds
            ht              =       mu*(F_g + F_bond + F_ve);   %Determine the deterministic displacement due to the force and torque working on the particle
        else                                                    %If there are no bonds, there is only gravity
            ht              =       mu*F_g;                     %Determine the deterministic displacement due to the force and torque working on the particle
        end

    end
    
    vrand           =       chol((4/dt)*mu,'lower')*X1(:,i);    %Determine the noise term due to thermal fluctations
    dr_p            =       0.5*dt*(ht + vrand);                %Calculate the predicted displacement
    Rot             =       RotationMatrix(dr_p(4:6));          %Determine the corresponding rotation matrix
    
    %Apply reflective boundary conditions on the predicted position:
    if (r(3) + dr_p(3)) <= a_base                            %Use reflective boundary condition 
        dr_p(3) = dr_p(3) + 2*abs(a_base-(r(3)+dr_p(3)));    
    end
    
    %Determine predicted positions:
    r_p                 =       r + dr_p(1:3);              %predicted position;
    e_binder_p          =       Rot*e_binder;               %predicted orientation of the binder
    e_protrusion_p      =       Rot*e_protrusion;           %predicted orientation of the protrusions
    
    test = r_p(3) + e_protrusion_p(3,:) - d < 0;            %Test contains which protrusions cross the boundary (Corresponding index value = 1)
       
    if sum(test) == 0 || n_protrusions == 0    %If there was no overlap or there are no protrusions
        
        mu_p            =       mu_wall_Jones_red(a_eff,r_p(3)+slip,mu_tt,mu_rr,mu_rt,A,B,C);      %Determine mobility matrix at the predicted position

        %Determine forces and torques due to bonds at predicted position
        if State > 0            %If there are bonds, calculate corresponding forces and torques
            %Direct influence of the bond:
            f_bond          =       -k_bond*(r_p(1:3) + e_binder_p - anch_point);    %First, determine all forces working on the particle    
            F_bond          =       [sum(f_bond,2); sum(cross(e_binder_p,f_bond),2)];                               %Add up al forces and torques due to the bonds

            f_ve            =       [0;0;1].*(kBT./(r_p(3) + e_binder_p(3,:)));
            
            %check whether particle is sufficiently close (Untick these
            %lines in case of bonds which are two times larger than the
            %radius of the sphere; this is not the case for f-BPM most of
            %the times).
            %for Q = 1:State          
            %    if r_p(3) + e_binder_p(3,bond_index(Q)) > 2*a_eff
            %        f_ve(:,Q) = [0;0;0];
            %    end
            %end
            
            F_ve           =       [sum(f_ve,2); sum(cross(e_binder_p,f_ve),2)];       %Put all in the same vector
            hp              =       mu_p*(F_g + F_bond + F_ve);                                   %Determine the deterministic displacement due to the force and torque working on the particle
        else                                                    %If there are no bonds, there is only gravity
            hp              =       mu_p*F_g;                   %Determine the deterministic displacement due to the force and torque working on the particle
        end

        %determined corrected displacement: 
        dr       =       dt*(hp + sqrt(1/dt)*mu_p*chol(Z,'lower')*(X1(:,i) + X2(:,i)));            %determine the correct displacement

        %Apply reflective boundary condition
        if (r(3) + dr(3)) <= a_base
            dr(3) = dr(3) + 2*abs(a_base-(r(3)+dr(3))); 
        end

        %Update position and orientation (test, to see if there is overlap of
        %the protrusions with the boundary)
        r_test              =       r + dr(1:3);                %Position update
        Rot                 =       RotationMatrix(dr(4:6));    %Determine rotation matrix
        e_protrusion_test   =       Rot*e_protrusion;           %Update protrusion location w.r.t. C.O.V. bead (rotate)

        %Now, do some testing on whether the new position is accepted (If
        %protrusion(s) cross the boundary, we will need to redo the iteration)
        test = r_test(3) + e_protrusion_test(3,:) - d < 0;                                            %Test contains which protrusions cross the boundary (Corresponding index value = 1)
    end 
    
    if sum(test) == 0 || n_protrusions == 0        %If there is no overlap or there are no protrusions, continue
        failure = 0;            %Reset the failure counter

        %Set the new locations to the accepted location
        r               =   r_test;               
        e_binder        =   Rot*e_binder;               
        e_protrusion    =   e_protrusion_test;

        if State_change == 1 && t > t_wait && r(3) - a_base < Rc + 2*d_max %If we are prone to a new bond in the next iteration
            
            %We have waited long enough (t>t_wait_min)
            %The particle is sufficiently close to the surface  (r(3) <= e_max + Rc)
            
            Q   = 0; %Use Q to keep track of the number of points that possibly might give a bond
            
            %Check which binders on surface are sufficiently close to the
            %bead while using the periodicity of the latice
            r_surface_potential = zeros(3,1);
            e_binder_potential  = zeros(3,1);
            %For all possible points on the surface
            for nx = (floor((r(1) - Rc - a_base)/l_0)-1):1:(floor((r(1) + Rc + a_base)/l_0)+1) %Check all points surrouning the bead
                for ny = (floor((r(2) - Rc - a_base)/l_0)-1):1:(floor((r(2) + Rc + a_base)/l_0)+1)
                    if vecnorm(r - [nx*l_0;ny*l_0;0]) <= (a_base + Rc) %If the distance between the sphere and the point on the surface is sufficiently small
                        %Check overlap of this point with any protrusion
                        e_bind  =   a_base*(([nx*l_0;ny*l_0;0] - r)/vecnorm([nx*l_0;ny*l_0;0] - r)); %Determine the corresponding
                        
                        if n_protrusions == 0 %If there are no protrusions
                            Q = Q + 1;  %Remember that another point has potential to form a bond
                            r_surface_potential(:,Q)    =   [nx*l_0;ny*l_0;0];
                            e_binder_potential(:,Q)     =   e_bind;
                        else
                            if sum(vecnorm(e_bind - e_protrusion)<d) == 0 %If there is no overlap at all with any of the protrusions
                                Q = Q + 1;  %Remember that another point has potential to form a bond
                                r_surface_potential(:,Q)    =   [nx*l_0;ny*l_0;0];
                                e_binder_potential(:,Q)     =   e_bind;
                            end
                        end
                    end
                end
            end
            
            %Check which binders on surface are sufficiently close to
            %protrusions on the bead while using the periodicity of the latice
            r_protrusion    =   r + e_protrusion;                   %Determine position of center of protrusion in the stationary frame of reference
            for nx = (floor((r(1) - Rc - a_base - 2*d_max)/l_0)-1):1:(floor((r(1) + Rc + a_base + 2*d_max)/l_0)+1) %Check all points surrouning the bead
                for ny = (floor((r(2) - Rc - a_base - 2*d_max)/l_0)-1):1:(floor((r(2) + Rc + a_base + 2*d_max)/l_0)+1)
%                     find((vecnorm([nx*l_0;ny*l_0;0] - r_protrusion)<(d+Rc))==1)
                    for QQ = find((vecnorm([nx*l_0;ny*l_0;0] - r_protrusion)<(d+Rc))==1)
                        Q = Q + 1;  %Remember that another point has potential to form a bond
                        r_surface_potential(:,Q)    =  [nx*l_0;ny*l_0;0];
                        vec                         =  r_protrusion(:,QQ) + d(QQ)*(r_surface_potential(:,Q) - r_protrusion(:,QQ))/vecnorm((r_surface_potential(:,Q) - r_protrusion(:,QQ)));
                        e_binder_potential(:,Q)     =  vec - r;
                    end
                end
            end
            
            %check whether binder locations on the surface are occupied
            del_index = zeros(1,Q);
            if State ~= 0 %If there is already a bond present
                QQ = Q;
                %Check overlap of binders at the surface of earlier bonds
                for c = 1:QQ %For all potential binder locations on the surface
                    for n = 1:State %For all currently locations already present
                        if sum(r_surface_potential(1:2,c) == anch_point(1:2,n)) == 2 %if one of the potential anchoring points is the same as a present anchoring point
                           Q = Q - 1;       %Update that we have less points less
                           del_index(c) = c; %Save which potential binders and anchoring points we need to delete
                        end
                    end
                end
            end
            
            r_surface_potential(:,nonzeros(del_index))  =   [];
            e_binder_potential(:,nonzeros(del_index))   =   [];

            if Q ~= 0 %If there are potential bonds left
                QQ = randi([1 Q]); %Select one at random
                
                State = State + State_change; %Update the state we are currently in
                
                %Set the next state_change and waiting time accordingly
                if State == 1 %If we are currently in the single-bound state
                    Q = rand(); %Draw a random number
                    if Q < Prob_S_D
                        State_change = 1;
                        t_wait = t + GenerateWaitTime(k_on_2 + k_off);
                    else
                        State_change = -1;
                        t_wait = t + GenerateWaitTime(k_on_2 + k_off);
                    end
                elseif State == 2
                    State_change = -1;
                    t_wait = t + GenerateWaitTime(2*k_off);
                elseif State == 0
                    State_change = 1;
                    t_wait = t + GenerateWaitTime(k_on_1);
                end
                
                %Set time step for the next iteration accordingly
                if State >= 2
                    dt = dt_b2;                 %If now two bonds or more are present, set dt = dt_b2 (Small time step)
                else 
                    dt = dt_b1;                 %If now one bond is present, set dt = dt_b1 (Medium sized time step)
                end
                
                %Save the binder orientation vector anchoring point
                e_binder(:,State)     =   e_binder_potential(:,QQ);
                anch_point(:,State)   =   r_surface_potential(:,QQ);
                
                %update the user on what has happened
                disp(['A bond has formed (State = ',num2str(State),' bonds).'])
                
                %Store how large the error was in association time:
                AssociationError(WW)     =  deltaT;   %Store how large the error is
                deltaT                   =  0;        %Restart measuring the error
                WW                       =  WW +1;    %Store that we will add another number next iteration
            else
                deltaT = deltaT + dt; %Store how long the error in association is
            end
        end

    %check bond disociation of currently formed bonds
    if State_change == -1 && t > t_wait %when we are due for dissociation of a bond
            Q   = randi([1 State]);                    %Select at random one of the formed bonds
            anch_point(:,Q)    =   [];                 %Delete the corresponding binder surface location 
            e_binder(:,Q)      =   [];                 %Delte the corresponding binder location on the particle                   

            State = State + State_change;
            %Set time step for the next section accordingly
            if  State >= 2
                dt = dt_b2;                 %If two bonds or more are present, set dt = dt_b2 (Small time step)
            elseif State == 1 
                dt = dt_b1;                 %If one bond is present, set dt = dt_b1 (Medium sized time step)
            else
                dt = dt_f;                  %If no bonds are present, set dt = dt_f (Just enough to average over 10 displacements in the motion blur)
            end
            
            %Set the next state_change and waiting time accordingly
            if State == 1 %If we are currently in the single-bound state
                Q = rand(); %Draw a random number
                if Q < Prob_S_D
                    State_change = 1;
                    t_wait = t + GenerateWaitTime(k_on_2 + k_off);
                else
                    State_change = -1;
                    t_wait = t + GenerateWaitTime(k_on_2 + k_off);
                end
            elseif State == 2
                State_change = -1;
                t_wait = t + GenerateWaitTime(2*k_off);
            elseif State == 0
                State_change = 1;
                t_wait = t + GenerateWaitTime(k_on_1);
            end
            
            %Finally, update user with what has happened:
            disp(['A bond has dissociated (State = ',num2str(State),' bonds).'])
    end
    
    %Data acquisition, with and without motion blurring
    if t >= j*dT - t_e/2 && t <= j*dT + t_e/2         %During the acquisition time,
        r_blurr =  r_blurr + r;                       %Add up al localizations
        m       =  m + 1;                             %And store the number of localizations
    end

    if floor(t/dT) == j                 %Every dT (acquisition)
        r_noblur(:,j) = r;              %save the position without motion blurring for comparison.
        Bonds(j)      = -State;   %Save the number of bonds
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
    Q               =       1;              %Save that we ended at a new iteration, so at the start we need to calculate forces again.
    
    else  %If overlap between protrusion and the surface, redo the iteration
        failure     =       failure + 1;    %Set the failure counter
        Q           =       0;              %Save that we start again at the same iteration, such that we do not need to calculate the forces again at the start of the iteration.
        i           =       i + 1;          %Use now a different random number
        
        if failure == 1e5 && State >=1     %If failed 10000 times and at least one bond
            %Delete at random one of the bonds
            Q                              =    randi([1 State]);        %Select at random one of the already formed bonds
            anch_point(:,Q)    =    [];                                             %Delete corresponing ligand location
            e_binder(:,Q)      =    [];                      %Delete corresponding bond index

            State = State - 1; %Save the current amount of states, briefly
            
            %Set the next state_change and waiting time accordingly
            if State == 1 %If we are currently in the single-bound state
                Q = rand(); %Draw a random number
                if Q < Prob_S_D
                    State_change = 1;
                    t_wait = t + GenerateWaitTime(k_on_2 + k_off);
                else
                    State_change = -1;
                    t_wait = t + GenerateWaitTime(k_on_2 + k_off);
                end
            elseif State == 2
                State_change = -1;
                t_wait = t + GenerateWaitTime(2*k_off);
            elseif State == 0
                State_change = 1;
                t_wait = t + GenerateWaitTime(k_on_1);
            end

            
            %Finally, update user with what has happened:
            disp(['A bond has dissociated (Fail safe) (b = ',num2str(State),' ).'])
            disp('If this happens a lot, you might want to consider to decrease the integration time of your simulation')
            failure                        =   0;                             %And reset the failure counter       
            Q                              =   1;                             %And save that we should recalculate the forces again at the start of the iteration.
        end
    end     
end
toc()
disp(['The mean error in association time is equal to ',num2str(mean(AssociationError(1:WW-1))),' seconds'])
disp(['This is equal to ',num2str(mean(AssociationError(1:WW-1))/dT),' frames'])

%% V. Display simulation result
t           = (1:N)*dT; %total t          (s)

%plot plot x-y and z-position in seperate graphs in time
figure(2)
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
plot(t,(r_blurred(3,:)-a_base)*10^(6))
hold on
title('z-position in time (Blurred)')
ylabel('z-position (\mum)')

%Plot the state of particle in time
subplot(4,1,4)
plot(t,Bonds);
ylim([min(Bonds)-0.1 0.1])
title('Particle state in time')
xlabel('Time (s)')
ylabel('State')

%plot the trajectory in a 3D plot
figure(3)
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
 
J = 1;
n = 1;
while n < N
    if Bonds(n+1) ~= Bonds(n) 
        if Bonds(n) == 0
            disp('Press [Space] to continue')
            figure(4)
            hold on
            plot(r_blurred(1,J:1:n)*10^6,r_blurred(2,J:1:n)*10^6,'k .')
            axis equal
            pause()
            J = n+1;
        elseif Bonds(n) == -1
            disp('Press [Space] to continue')
            figure(4)
            hold on
            plot(r_blurred(1,J:1:n)*10^6,r_blurred(2,J:1:n)*10^6,'r .')
            axis equal
            J = n+1;
            pause()
        elseif Bonds(n) == -2
            disp('Press [Space] to continue')
            figure(4)
            hold on
            plot(r_blurred(1,J:1:n)*10^6,r_blurred(2,J:1:n)*10^6,'g .')
            axis equal
            xlabel('x-position (\mum)')
            ylabel('y-position (\mum)')
            J = n+1;
            pause()
        elseif Bonds(n) < -2
            disp('Press [Space] to continue')
            figure(4)
            axis equal
            hold on
            plot(r_blurred(1,J:1:n)*10^6,r_blurred(2,J:1:n)*10^6,'b .')
            J = n+1;
            pause()
        end
    end
    n = n + 1;
end

disp('Press [Space] to continue')
       if Bonds(n) == 0
            figure(4)
            hold on
            plot(r_blurred(1,J:1:n)*10^6,r_blurred(2,J:1:n)*10^6,'k .')
            axis equal
            title('Motion (blurred) in x-y plane (k = 0, r = 1, g = 2, b = > 2)')
            xlabel('x-position (\mum)')
            ylabel('y-position (\mum)')
            pause()
       elseif Bonds(n) == -1
            figure(4)
            hold on
            plot(r_blurred(1,J:1:n)*10^6,r_blurred(2,J:1:n)*10^6,'r .')
            axis equal
            title('Motion (blurred) in x-y plane (k = 0, r = 1, g = 2, b = > 2)')
            xlabel('x-position (\mum)')
            ylabel('y-position (\mum)')
            pause()
       elseif Bonds(n) == -2
            figure(4)
            hold on
            plot(r_blurred(1,J:1:n)*10^6,r_blurred(2,J:1:n)*10^6,'g .')
            axis equal
            title('Motion (blurred) in x-y plane (k = 0, r = 1, g = 2, b = > 2)')
            xlabel('x-position (\mum)')
            ylabel('y-position (\mum)')
            pause()
       elseif Bonds(n) < -2
            figure(4)
            hold on
            plot(r_blurred(1,J:1:n)*10^6,r_blurred(2,J:1:n)*10^6,'b .')
            axis equal
            title('Motion (blurred) in x-y plane (k = 0, r = 1, g = 2 , b = > 2)')
            xlabel('x-position (\mum)')
            ylabel('y-position (\mum)')
            pause()
        end

% Plot xy position in time (blurred as measured in the experimental set-up)
figure(5)
subplot(2,1,1)
hold on
plot(t,r_blurred(1,:)*10^6)
plot(t,r_blurred(2,:)*10^6)
xlim([0 max(t)])
ylabel('Position (\mum)')
xlabel('Time (s)')
legend('x-position','y-position')
title('Motion in time (blurred)')

subplot(2,1,2)
plot(t,Bonds)
xlabel('Time (s)')
ylabel('State')
ylim([-2.1 0.1])
yticks([-2 -1 0])
title('State in time')

figure(6)
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
plot(t-dT*MeasurementWindow/2,Bonds)

figure(7)
hold on
histogram(Dcoef_blurred*10^(12),75)
xlim([0, 1.1*max(Dcoef_blurred*10^(12))])
xlabel('Diffusion Coefficien(\mum^2/s)')
ylabel('Counts (-)')
title('Measured Diffusion Coefficient')

%% VI Save the result in a certain format
x_position_nm = r_blurred(1,:)';
y_position_nm = r_blurred(2,:)';
N_Bonds = -Bonds';
T = table(x_position_nm, y_position_nm, N_Bonds);
writetable(T,'Result_fBPM_Simulation.csv') %Save in excel sheet 
%column 1 x-position, column 2 y-position, column 3 State (Number of bonds)