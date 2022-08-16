function mu = mu_wall_Jones_red(a,h,mu_tt,mu_rr,mu_rt,A,B,C)
    t = a/h;
    T = [1 t t^2 t^3 t^4 t^5 t^6 t^7 t^8 t^9];
    F1 = C(:,1)*(t/(1-t)) + C(:,2)*log(1-t) + C(:,4)*((1-t)/t)*log(1-t);
    phi_tt = F1(1) + dot(T,A(1,:))/dot(T,B(1,:));
    psi_tt = F1(2) + dot(T,A(2,:))/dot(T,B(2,:));
    phi_rr = F1(3) + dot(T,A(3,:))/dot(T,B(3,:));
    psi_rr = F1(4) + dot(T,A(4,:))/dot(T,B(4,:));
    psi_tr = F1(5) + dot(T,A(5,:))/dot(T,B(5,:));
    
    alpha_tt    = mu_tt*(1/phi_tt);
    alpha_rr    = mu_rr*(1/phi_rr);
    beta_tt     = mu_tt*(psi_rr/(psi_tt*psi_rr - (4/3)*psi_tr^2));
    beta_rr     = mu_rr*(psi_tt/(psi_tt*psi_rr - (4/3)*psi_tr^2));
    beta_tr     = mu_rt*(- psi_tr/((3/4)*psi_tt*psi_rr - psi_tr^2));
    
    mu = [[beta_tt 0 0 0 beta_tr 0] ; [0 beta_tt 0 -beta_tr 0 0] ; [0 0 alpha_tt 0 0 0] ; [0 -beta_tr 0 beta_rr 0 0] ; [beta_tr 0 0 0 beta_rr 0] ; [ 0 0 0 0 0 alpha_rr]];
end