%This script allows generating a figure similar to [REF, Fig.2]
%The performance curves corresponding to the CSI-free schemes are not
%displayed with this code. 

%% System parameters
M = 8;        %Number of PB antennas
N = 8;        %Number of EH devices
radius = 10;  %The EH devices are assumed to be randomly and uniformly 
              %distributed around the PB in a circular area of radius 10m.
samples = 250;%Number of Monte Carlo runs 

%Algorithm 1 parameters [REF]:
delta = 0.9;
epsilon = 1e-6;

% Results are drawn as a function of the Rician LOS factor kappa
Kappa_dB = -10:5:30;            % dB scale
Kappa = 10.^(Kappa_dB./10);     % linear scale

avgCSI_mrt = zeros(samples,length(Kappa));
fullCSI_sdp = zeros(samples,length(Kappa));
avgCSI_sdp = zeros(samples,length(Kappa));
for k = 1:length(Kappa)
    kappa = Kappa(k);
    for j=1:samples   
        [Kappa_dB(k) j]        
        rand('seed',j)
        
        d_i = radius*sqrt(rand(N,1));  %distance between the i-th device and the PB
        theta_i = 2*pi*rand(N,1);      %azimuth angle of the i-th device relative to
                                       %the boresight of the PB's antenna array
        beta = 1e-2*max(d_i,1).^(-2.7);%path loss model
        beta = sort(beta,'descend');

        %Rician channel samples
        h_los = zeros(N,M);
        h_nlos = zeros(N,M);
        for i=1:N
            phi = -(0:(M-1))*pi*sin(theta_i(i));
            h_los(i,:) = sqrt(kappa/(1+kappa))*exp(1i*pi/4)*exp(1i*phi);   %LOS component
            h_nlos(i,:) = sqrt(1/(2*(1+kappa)))*(randn(1,M)+1i*randn(1,M));%NLOS component
        end
        h = h_los + h_nlos;

        %MRT precoder using average CSI
        [~, E] = avgCSI_MRT(M,N,h,h_los,beta,delta,epsilon);  %Precoder computed using average CSI
        avgCSI_mrt(j,k) = E;  %instantaneous minimum RF energy available at the EH devices

        %average-CSI based  SDP 
        [W, ~] = CSIBeamf_SDP(M,N,h_los,beta); %Precoder computed using average CSI
        [r,~]=size(W);
        E = zeros(1,N);
        for ii = 1:N
            for jj = 1:r
                E(ii) = E(ii) + beta(ii)*abs(sum(W(jj,:).*h(ii,:)))^2;
            end
        end
        avgCSI_sdp(j,k) = min(E);  %instantaneous minimum RF energy available at the EH devices

        %full-CSI based  SDP   
        [~, E] = CSIBeamf_SDP(M,N,h,beta);  %Precoder computed using instantaneous CSI
        fullCSI_sdp(j,k) = min(E);          %instantaneous minimum RF energy available at the EH devices
      
    end
        
end

%Results in dBm
AVG_CSI_SDP = 10*log10(nanmean(avgCSI_sdp))+30;
FULL_CSI_SDP = 10*log10(nanmean(fullCSI_sdp))+30;
AVG_CSI_MRT = 10*log10(nanmean(avgCSI_mrt))+30;

% Plot of the results
figure
set(gcf, 'Units', 'centimeters'); 
axesFontSize = 16;
legendFontSize = 16;
afFigurePosition = [2 7 19 12]; 
set(gcf, 'Position', afFigurePosition,'PaperSize',[18 8],'PaperPositionMode','auto'); % [left bottom width height], setting printing properties 

%%
h1=plot(Kappa_dB,FULL_CSI_SDP,'-o','MarkerSize',6,'LineWidth',2.5); hold on
h2=plot(Kappa_dB,AVG_CSI_SDP,'-s','MarkerSize',8,'LineWidth',2.5); hold on
plot(Kappa_dB,AVG_CSI_MRT,'-x','MarkerSize',12,'LineWidth',2.5); hold on

hl=legend('Full-CSI (optimum)','Average-CSI (optimum)', 'Average-CSI (Alg.~1)');
set(hl,'interpreter','latex','FontSize',legendFontSize);
set(gca,'FontSize',12)

box on
grid on
ylim([-25 -5])

xlabel('$\kappa$ (dB)','Interpreter', 'latex','fontsize',axesFontSize)
ylabel('average $\inf\{E_i\}$ (dBm)','Interpreter', 'latex','fontsize',axesFontSize)

%% References:
%[REF]    - O. L. A. López, F. A. Monteiro, h. Alves, R. Zhang and M. Latva-Aho,
%           "A Low-Complexity Beamforming Design for Multiuser Wireless Energy 
%           Transfer," in IEEE Wireless Communications Letters, vol. 10,
%           no. 1, pp. 58-62, Jan. 2021, doi: 10.1109/LWC.2020.3020576.