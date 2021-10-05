function [Wprecod, E] = avgCSI_MRT(M,N,H,H_los,beta,delta,epsilon)
%Beamforming design (and corresponding receive energy) for an average
%CSI precoding using MRT exploiting average CSI [REF, Algorithm 1]. 

%% Input parameters:
% M    -> Number of Tx Antennas in the PB  (scalar)
% N    -> Number of EH devices
% H    -> NxM normalized complex instantaneous channel matrix 
% H_los-> NxM normalized complex LOS (average) channel matrix 
% The M-dimentional channel vector corresponding to the i-th device appears
% in the i-th row of H or H_los
% beta -> N-dimentional path loss vector. beta(i) corresponds to the large-
%         scale path loss associated to the i-th device

%% Output parameters:
% Wprecod-> rxM precoding matrix. r<=M and it is found here as a by-product 
% E      -> minimum RF energy available for the set of EH devices

   
%% Main code

% Set Q as defined after [REF, eq.7]
for j=1:N
    for i=1:N
        Q(j,i) = abs(conj(H_los(j,:))*H_los(i,:).')^2;
        Q(j,i) = Q(j,i)/norm(H_los(j,:))^2;
    end
end
Q=Q./max(Q,[],'all');

% Initialization of Algorithm 1 in [REF]
B = diag(beta);
A = [ones(N,1) -B*Q' eye(N,N); 0 ones(1,N) zeros(1,N)];
mu = [-1; zeros(2*N,1)];
p_0 = 1./(beta*sum(1./beta));  % initial MRT power allocation [REF, eq.10]
xi_0 = min(B*Q'*p_0);          % initialization of xi_0
z = [xi_0; p_0; B*Q'*p_0-xi_0];% initialization of z


iterations=0;
err = 1; % accuracy initialization
%Main body of Algorithm 1 in [REF], codelines: 5-10
while err>epsilon && iterations < 100 %Stop only when reaching the required
                                      %accuracy epsilon or maximum number of
                                      %iterations=100
    %i
    Z = diag(z);  % Algorithm 1, code line: 6
    lambda = (A*Z^2*A')\(A*Z^2*mu);   % Algorithm 1, code line: 6
    r = mu - A'*lambda;               % Algorithm 1, code line: 7
    if min(r) >= -1e-4                %accuracy update Algorithm 1, code line: 7 and 10. 
                                      %Instead of a comparison with 0, we
                                      %adopt a relatively small negative
                                      %number as suggested in [REF, footnote 3]
        err = ones(1,2*N+1)*Z*abs(r); 
    end

    z = z - delta*Z^2*r/norm(Z*r);    % update of z. Algorithm 1, code line 8
    iterations = iterations + 1;
end
p = z(2:N+1);    % Solution of Algorithm 1 (code line 11)
    
% Set the MRT precoder using [REF, eq.6]
for i=1:N
    Wprecod(i,:) = sqrt(p(i)).*conj(H_los(i,:))./norm(H_los(i,:));
end

% the following code computes the RF energy available at each device using
% the instantaneus channel realizations
E=zeros(1,length(beta));
for i = 1:length(beta)  
    E(i) = beta(i)*norm(H(i,:)*Wprecod.')^2;
end
% Minimum RF energy at the EH devices.
E=min(E);


%% References:
%[REF]    - O. L. A. López, F. A. Monteiro, H. Alves, R. Zhang and M. Latva-Aho,
%           "A Low-Complexity Beamforming Design for Multiuser Wireless Energy 
%           Transfer," in IEEE Wireless Communications Letters, vol. 10,
%           no. 1, pp. 58-62, Jan. 2021, doi: 10.1109/LWC.2020.3020576.
end

