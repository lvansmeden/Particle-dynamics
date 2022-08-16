function [d, e_protrusions, V_p] = GenerateParticle_nobinders(d_min, d_max, n_protrusions, a, offset_max)
%Input:
%   d_min:         Minimum radius of a protrusion
%   d_max:         Maximum radius of a protrusion
%   n_protrusions: The number of protrusions on the particle
%   n_binders:     The number of binders on the particle (On the top of a
%                  protrusion or on the surface of the particle)
%   a:             Radius of the sphere/bead on which we will put
%                  protrusions and binders
%   offset_max     %Number between 0 and 1, which indicates by how much the
%                  protrusion can stick out or can be embedded in/on top of the bead

%Output:
%   d:             An array which contains the raddi of the protrusions which are put
%                  on the bead
%   e_protrusions: Matrix which contains in it columns the orientation
%                  vectors of the protrusion. An orientation vector points from the center
%                  of the sphere to the center of a protrusion
%   V_p:           %Volume of the bead + the non-shared volume of the
%                  protrusions. This is called the effective volume and is
%                  used to calculate the effective radius of the bead and
%                  the gravity working on the bead.

NP              =     1;                             %USe NP to keep track of the number of protrusions we have generated
d               =     zeros(1,n_protrusions);        %Use d to store the raddii of the protrusions we create
offset          =     zeros(1,n_protrusions);        %Use offset to store the distance from the center of the bead to the center of the protrusions
e_protrusions   =     zeros(3,n_protrusions);        %Prelocate for faster computation

failure         =     0;                             %Keep track of the number of failures (If too many failures happen, we start generating new protrusions)

%Algorithm is partly based on the BEP Report of Jakob Limpens
%'The effect of surface roughness on motion patterns in TPM' (2020)


%% 1) First, determine the locations of the protrusions (we do not allow overlap of the spheres)
while NP <= n_protrusions                                               %While still having protrusions to generate
    d(NP)      =   d_min + rand*(d_max-d_min);                          %Generate a random radius for the protrusion
    offset(NP) =   a + offset_max*(- d(NP) + 2*d(NP)*rand());           %Generate at random the offset of the protrusion
    phi        =   2*pi*rand();                                         %Create a random phi value, to assign the protrusion to the surface of the sphere
    theta      =   acos(2*rand()-1);                                    %Create a random theta value, to assign the protrusion to the surface of the sphere

    e_protrusions(:,NP)   =     offset(NP)*[sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]; %This is the location of the center of the protrusion

    Q = 0; %Set Q = 0 at first, which is, assume that there is no overlap
    %Check the candidate, so that it does not overlap with earlier protrusions
    for i = 1:1:NP-1
        if abs(vecnorm(e_protrusions(:,NP)-e_protrusions(:,i))) <= d(NP) + d(i) %If the distance between the center is smaller than the combined radius
            Q = 1;  %There is overlap, so set Q = 1
            break   %Stop the loop, since there is no reason for checking the other protrusions
        end
    end

    if Q == 0 %If the protrusion was okay
        NP = NP + 1;
        failure = 0; %Also, reset the fail-counter
    else %If the protrusion overlapped with other protrusions,redo this protrusion, and add up the failure counter
        failure = failure + 1;
    end

    if failure == 1e5 %Stop after many attempts
        disp(['Due to overlap, we stopped at ',num2str(NP),' protrusions, (You wanted ',num2str(n_protrusions),' protrusions.)'])
        if n_protrusions ~= 0
            e_protrusions(:,NP:end)     =   [];             %Throw away the zero entries and single wrong entry
            d                           =   d(1:NP-1);      %Do the same for the raddii
        end
        break
    end
end

%% 2) Determine total volume of the sphere and the protrusions
V_s         =       (4/3)*pi*a^3;   %Volume of the sphere
V_add       =       zeros(1,NP-1);  %Create array to save the non-shared volumes

%Calculation based on formulas 6 7 and 8 of the BEP Report of Jakob Limpens
%'The effect of surface roughness on motion patterns in TPM' (2020)
for i = 1:NP-1 %For all protrusions
    if d(i) > 0 && offset(i) > 0
        beta = acos((a^2 + offset(i)^2 - d(i)^2)/(2*a*offset(i))); %theta
        alpha  = acos((d(i)^2 + offset(i)^2 - a^2)/(2*offset(i)*d(i))); %alpha

        V_shared        =   (pi/3)*(a^3)*(2 - 3*cos(beta) + cos(beta)^3) + (pi/3)*(d(i))^3*(2-3*cos(alpha) + cos(alpha)^3); %calculate the shared volume of the bead and the protrusion
        V_add(i)        =   (4/3)*pi*d(i)^3 - V_shared;                                                                     %Determine the volume we thus need to add on top of that of the sphere
    end
end

V_p         =       V_s    + sum(V_add);  %Total volume of the particle including protrusions
