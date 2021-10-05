function [Wprecod, E] = CSIBeamf_SDP(M, N, H, beta)
%Beamforming design (and corresponding receive energy) for a full or average
%CSI precoding [REF]. 

%% Input parameters:
% M    -> Number of Tx Antennas in the PB  (scalar)
% N    -> Number of EH devices
% H    -> NxM normalized complex instantaneous/average channel matrix 
%        (the M-dimentional channel vector corresponding to the i-th device
%         appears in the i-th row of H)
% beta -> N-dimentional path loss vector. beta(i) corresponds to the large-
%         scale path loss associated to the i-th device

%% Output parameters:
% Wprecod-> rxM precoding matrix. r<=M and it is found here as a by-product 
% E      -> minimum RF energy available for the set of EH devices

   
%% Main code

% Code for creating {H_i}_i in [REF, eq.4b] as 3D matrix
for i=1:N
   Hi(:,:,i) = (H(i,:).')*(H(i,:).')';
end

% SDB solution. This requires CVX solver (http://cvxr.com/cvx/)
% The problem is formulated as it appears in [REF, eq.4]
cvx_begin sdp quiet
    variable W(M,M) hermitian semidefinite
    variable xi
    minimize( - xi )  
    for i=1:N           
        beta(i)*real(trace((Hi(:,:,i)*W))) >= xi
    end
    trace(W) == 1;    
cvx_end

r = rank(W);      % the rank determines the number of precoders
% the precoders match the eigenvectors of W and are computed as follows
[U,S] = eig(W);
for j = 1:r
    Wprecod(j,:) = sqrt(S(j,j))*U(:,j)';
end
Wprecod = Wprecod./norm(Wprecod,'fro');

% the following code computes the RF energy available at each device
E = zeros(1,N);
for i = 1:N
  for j = 1:r
       E(i) = E(i) + beta(i)*abs(sum(Wprecod(j,:).*H(i,:)))^2;
  end
end

% the minimum RF available energy at the devices is given next. 
E = min(E);
% E should match the value of xi from the solution of the SDP


%% References:
%[REF]    - O. L. A. López, F. A. Monteiro, H. Alves, R. Zhang and M. Latva-Aho,
%           "A Low-Complexity Beamforming Design for Multiuser Wireless Energy 
%           Transfer," in IEEE Wireless Communications Letters, vol. 10,
%           no. 1, pp. 58-62, Jan. 2021, doi: 10.1109/LWC.2020.3020576.
