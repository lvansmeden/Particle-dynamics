function t_wait = GenerateWaitTime(k)
    %k: rate at which association/dissociation occurs (k_on or k_off)
    %State1: state of the particle before that the event happened
    %State2: State of the particle afer that the event happened
    Q = 0; %Use Q to keep track of whether the new waiting time is okay
    while Q ==0 %While the waiting time is not okay yet
        t_wait = 100*(1/k)*rand; %Generate a random waiting time
        y = k*rand; %Generate a random number
        
        if y <= k*exp(-k*t_wait) %Sample
            Q = 1; 
        end
    end
end