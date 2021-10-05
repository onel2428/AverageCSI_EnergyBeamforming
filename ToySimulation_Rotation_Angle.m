%This script allows generating a figure similar to [REF, Fig.4]
%The performance curves correspond only to Scenario B

%% System parameters
M = 8;        %Number of PB antennas
N = 8;        %Number of EH devices
kappa = 10;   %LOS factor
samples = 250;%Number of Monte Carlo runs 

% Scenario B configuration
d_i = (2:9)';                %distance between the i-th device and the PB
theta_i = (80:-10:10)*pi/180;%azimuth angle of the i-th device relative to
                             %the boresight of the PB's antenna array
beta = 1e-2*d_i.^(-2.7);     %path loss

%Algorithm 1 parameters [REF]:
delta = 0.9;
epsilon = 1e-6;

% Results are drawn as a function of the angle rotation
A = -pi/2:pi/50:pi/2;
avgCSI_mrt = zeros(samples,length(A));
fullCSI_sdp = zeros(samples,length(A));
avgCSI_sdp = zeros(samples,length(A));
for  a = 1:length(A)
    alpha = A(a);    
    for j=1:samples
        [a j]
        
        %Rician channel samples
        h_los = zeros(N,M);
        h_nlos = zeros(N,M);
        for i=1:N
            phi = -(0:(M-1))*pi*sin(theta_i(i) + alpha);
            h_los(i,:) = sqrt(kappa/(1+kappa))*exp(1i*pi/4)*exp(1i*phi);     %LOS component
            h_nlos(i,:) = sqrt(1/(2*(1+kappa)))*(randn(1,M) + 1i*randn(1,M));%NLOS component
        end
        h=h_los+h_nlos;
        
        %MRT precoder using average CSI
        [~, E] = avgCSI_MRT(M,N,h,h_los,beta,delta,epsilon);  %Precoder computed using average CSI
        avgCSI_mrt(j,a) = E;  %instantaneous minimum RF energy available at the EH devices

        %average-CSI based  SDP 
        [W, ~] = CSIBeamf_SDP(M,N,h_los,beta); %Precoder computed using average/LOS CSI
        [r,~]=size(W);
        E = zeros(1,N);
        for ii = 1:N
            for jj = 1:r
                E(ii) = E(ii) + beta(ii)*abs(sum(W(jj,:).*h(ii,:)))^2;
            end
        end
        avgCSI_sdp(j,a) = min(E);  %instantaneous minimum RF energy available at the EH devices    
    end
end

%Results in dBm
AVG_CSI_MRT = 10*log10(nanmean(avgCSI_mrt))+30;
AVG_CSI_SDP = 10*log10(nanmean(avgCSI_sdp))+30;

% Plot of the results
figure
set(gcf, 'Units', 'centimeters'); % set units 
axesFontSize = 16;
legendFontSize = 16;
% setting size & position
afFigurePosition = [2 7 19 12]; % [pos_x pos_y width_x width_y]   19x12
set(gcf, 'Position', afFigurePosition,'PaperSize',[19 12],'PaperPositionMode','auto'); % [left bottom width height], setting printing properties 

h1=plot(A*180/pi,AVG_CSI_SDP,'-','MarkerSize',8,'LineWidth',2.5); hold on
h2=plot(A*180/pi,AVG_CSI_MRT,':','MarkerSize',8,'LineWidth',2.5); hold on

[p,m]=max(AVG_CSI_SDP);
h3=plot([A(m)*180/pi A(m)*180/pi+180],[p p],'d','MarkerSize',8,'LineWidth',2,'MarkerEdgeColor','r');

grid on
box on
xlim([-90 90])
xticks(-90:30:90)
ylim([-16 -6])

set(gca,'FontSize',12)

hl=legend([h1 h2 h3],'Average-CSI (optimum)','Average-CSI (Alg. 1)','Optimum orientation');
set(hl,'interpreter','latex','FontSize',legendFontSize);

xlabel('$\alpha^\circ$','Interpreter', 'latex','fontsize',axesFontSize)
ylabel('average $\inf\{E_i\}$ (dBm)','Interpreter', 'latex','fontsize',axesFontSize)