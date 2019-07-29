%%% INTEGRATED TR DATA vs. MD HEATFLUX %%%%%%%%%%
%length = [2; 25; 204; 500; 1000] %nm
%anh_pr = [3.48; 40.8; 95.9; 97.8; 106] %anharmonic pristine
%anh_pr_err = [0.32; 5.1; 3.8; 2.8; 3.6];
%anh_ints = [2.22; 5.1; 6.44; 7.96; NaN];
%anh_ints_err[0.11; 1.2; 0.5; 1.4; NaN];
%harm_pr = [3.7; 41.97; 99.41; 94.08; 100.75];
%harm_ints = [2.35; 5.5; 6.18; 8.78; NaN];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 

load('MFP-FinalVersion');

figure;

% ax3=subplot(1,1,1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(ax3,OM,prist(1,:),'b-+',OM,prist(2,:),'m-+',OM,prist(5,:),'k-+',...
%     'MarkerSize',3);
% hold on;
% plot(ax3,OM,ints(1,:)...
%     ,'b:',OM,ints(2,:),...
%     'm:',OM,ints(5,:),'k:','MarkerSize',2);
% lg3=legend('2nm-Si','25nm-Si','500nm-Si',...
%     '2nm-Si_{1.56% Ge}','25nm-Si_{1.56% Ge}','500nm--Si_{1.56% Ge}');
% lg3.FontSize = 14;
% % title('Spectral Transmission');
% axis([0 16 0 55]);
% % xlabel('Frequency, \omega (THz)','fontsize',20);
% ylabel('Tr(\omega)','fontsize',20);
% xticks(0:2:20);
% ax3.XGrid='on';
% ax3.GridLineStyle='--';
% ax3.GridAlpha = 0.25;
% ax3.XAxis.FontSize = 16;
% ax3.YAxis.FontSize = 16;



ax3=subplot(2,2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(ax3,OM,prist(1,:),'b+',OM,prist(2,:),'m+','MarkerSize',2);
hold on;
plot(OM,prist(3,:), '+', 'color', [.3 0.1 0.5],'MarkerSize',2)%'purple'
plot(OM,prist(4,:),'+', 'color', [0.1 0.5 0.1],'MarkerSize',2);%'g'
plot(OM,prist(5,:),'k+',OM,prist(6,:),'r+','MarkerSize',2)
% lg3=legend('2nm','25nm','41nm','204nm','500nm','1000nm');
% lg3.FontSize = 16;
axis([0 20 0 150]);
% xlabel('Frequency, \omega (THz)','fontsize',20);
ylabel('Tr(\omega)','fontsize',20);
xticks(0:2:20);
ax3.XGrid='on';
ax3.GridLineStyle='--';
ax3.GridAlpha = 0.25;
ax3.XAxis.FontSize = 16;
ax3.YAxis.FontSize = 16;

ax4=subplot(2,2,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(ax4,OM,ints(1,:),'b+',OM,ints(2,:),'m+','MarkerSize',2);
hold on;
plot(OM,ints(3,:), '+','color', [0.3 0.1 0.5],'MarkerSize',2);%'c'
plot(OM,ints(4,:), '+','color', [0.1 0.5 0.1],'MarkerSize',2);%'g'
plot(OM,ints(5,:),'k+','MarkerSize',2);
% lg4=legend('2nm','25nm','41nm','204nm','500nm');
% lg4.FontSize = 16;
axis([0 20 0 120]);
% xlabel('Frequency, \omega (THz)','fontsize',20);
ylabel('Tr(\omega)','fontsize',20);
xticks(0:2:20);
ax4.XGrid='on';
ax4.GridLineStyle='--';
ax4.GridAlpha = 0.25;
ax4.XAxis.FontSize = 16;
ax4.YAxis.FontSize = 16; 

ax1=subplot(2,2,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
semilogy(ax1,OM, 1e9*pMFP,'r+',OM, 6e3*om_fit1,'r--',... 
    OM,1e9*iMFP,'b+',OM,1.5e3*om_fit2,'b--','MarkerSize',2);
axis([0 20 0.2 5e4]);
xlabel('Frequency, \omega (THz)','fontsize',20);
ylabel('\Lambda(\omega) (nm)','fontsize',20);
xticks(0:2:20);
yticks([1e0 1e1 1e2 1e3 1e4])
ax1.XGrid='on';
ax1.GridLineStyle='--';
ax1.GridAlpha = 0.25;
ax1.XAxis.FontSize = 16;
ax1.YAxis.FontSize = 16;

ax2=subplot(2,2,4);
% ax2=subplot(1,1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(ax2,OM,pkom(2,:),'m-+',OM,pkom(5,:),'k-+'...
    ,OM,ikom(2,:),'m-.',OM,ikom(5,:),'k-.','MarkerSize',2);
axis([0 20 0 55]);
% lg2=legend('2nm-Si','25nm-Si','500nm-Si','2nm-Si_{1.56% Ge}',...
%     '25nm-Si_{1.56% Ge}','500nm-Si_{1.56% Ge}');
% lg2.FontSize = 16;
% title('Spectral Thermal Conductivity');
xlabel('Frequency, \omega (THz)','fontsize',20);
ylabel('\kappa(\omega) (W/m-K)','fontsize',20);
xticks(0:2:20);
% yticks([1e-1 1e0 1e1 1e2]);
ax2.XGrid='on';
ax2.GridAlpha=0.25;
ax2.GridLineStyle='--';
ax2.XAxis.FontSize = 16;
ax2.YAxis.FontSize = 16;


% close all;