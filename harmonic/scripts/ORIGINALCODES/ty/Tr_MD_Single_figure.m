%clear

%%%%%%%% Do convolution on data series %%%%%%%%%%%%%%%%%%%
% dwin=0.3e12*(3*pi); % THz change this value to adjust the convolution width
% 
% win=round(dwin/dom);
% g = gausswin_my(win); % <-- this value determines the width of the smoothing window
% g = g/sum(g);
% Tr1=conv(Jom_ave,g,'same');
% Dos_g=conv(Dos_ave, g, 'same'); % convolution for Dos using same width as Tr
% % Dos_L_g=conv(Dos_L_avg, g, 'same');
% % Dos_R_g=conv(Dos_R_avg, g, 'same');
% 
% dwin=0.1e12*(2*pi); % THz
% win=round(dwin/dom);
% g = gausswin_my(win); % <-- this value determines the width of the smoothing window
% g = g/sum(g);      
% Tr2=conv(Jom_ave,g,'same');
% 
% OM = oms_fft/2/pi*1e-12;
% data=[OM; Jom_ave; Tr1; Tr2; Dos_ave];
% data = real(data);
% fid = fopen(['Tr_', CS],'w');
% fprintf(fid,'%10.5f\t%f\t%f\t%f\t%f\n',data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dossp = Dossp./trapz(OM,Dossp);
Dossv = Dossv./trapz(OM,Dossv);
Doslp = Doslp./trapz(OM,Doslp);
Doslv = Doslv./trapz(OM,Doslv);

%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%
figure

%%%%%%%%%%% Smoothed Tr %%%%%%%%%%%%%%%%%%%%%%%%%

% ax1 = subplot(2,1,1);
plot(OM, datasp(3,:), 'r', OM, datasv(3,:), 'r:', OM, datalp(3,:), 'b', OM, datalv(3,:), 'b:'); %smoothed Tr
xlabel('Frequency, \omega (THz)');
ylabel('Transmission');
title('Carbon-diamond membranes', 'FontSize', 10);
legend(['1.1 nm' newline 'prist.'], ['1.1 nm' newline '1% vacs.'], ['72.3 nm' newline 'prist.'],['72.3 nm' newline '1% vacs.']);
axis([0 60 0 120]);
ax1.XGrid='on';
ax1.GridLineStyle='--';
ax1.GridAlpha = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Smoothed Dos %%%%%%%%%%%%%%%%%%%%%%%%%%%

% ax2 = subplot(2,1,2);
% plot(ax2, OM, Dossp, 'r', OM, Dossv, 'r:', OM, Doslp, 'b', OM, Doslv, 'b:'); %Plot DOS
% xlabel('Frequency, \omega (THz)');
% ylabel('DoS (a.u.)');
% % yticks([]);
% % title('Density of States', 'FontSize', 10);
% axis([0 60 0 0.06]);
% ax2.XGrid='on';
% ax2.GridLineStyle='--';
% ax2.GridAlpha = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Other Useful Plots %%%%%%%%%%%%%%%%%%%%%%%%

%plot(OM, Dos_L_g, 'k'); %smoothed Left Dos_ave
%plot(OM, Dos_L, 'k'); %Left Dos_ave
%plot(OM, Dos_R_g, 'k'); %smoothed Right Dos_ave
%plot(OM, Dos_R, 'k'); %Right Dos_ave
%plot(OM, Dos_g, 'k'); %smoothed Dos_ave
%plot(OM, data(5,:), 'k'); %raw Dos_ave
%plot(OM, data(4,:), 'k'); %raw Tr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%