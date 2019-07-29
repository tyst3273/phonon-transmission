% clear all;
load('allTR');

dom = 1.570796326794897e+10;

ints = zeros(6,6521);
prist = zeros(6,6521);
len = zeros(6,1);
OM = OM(1:6521);

ints(1,:) = -i2(3,1:6521);
ints(2,:) = -i25(3,1:6521);
ints(3,:) = -i41(3,1:6521);
ints(4,:) = -i204(3,1:6521);
ints(5,:) = -i500(3,1:6521);

prist(1,:) = -p2(3,1:6521);
prist(2,:) = -p25(3,1:6521);
prist(3,:) = -p41(3,1:6521);
prist(4,:) = -p204(3,1:6521);
prist(5,:) = -p500(3,1:6521);
prist(6,:) = -p999(3,1:6521);

len(1) = 1.1e-9; %nm to m
len(2) = 25e-9;
len(3) = 41.4e-9; 
len(4) = 204.8e-9;
len(5) = 501.1e-9;
len(6) = 1003e-9;

lamstr = ['2vs25','2vs41','2vs204','2vs500','2vs999',...
    '25vs41','25vs204','25vs500','25vs999',....
    '41vs204','41vs500','41vs999',...
    '204vs500','204vs999',...
    '500vs999'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lam_p = zeros(15,6521);
lam_i = zeros(15,6521);

for j = 1:6521
    
    %%%% prist
    
    lam_p(1,j) = ((len(2)-len(1))/... #2 vs 25
        ((prist(1,j)/prist(2,j))-1))...
        -len(1);
    lam_p(2,j) = ((len(3)-len(1))/... #2 vs 41
        ((prist(1,j)/prist(3,j))-1))...
        -len(1);
    lam_p(3,j) = ((len(4)-len(1))/... #2 vs 204
        ((prist(1,j)/prist(4,j))-1))...
        -len(1);
    lam_p(4,j) = ((len(5)-len(1))/... #2 vs 500
        ((prist(1,j)/prist(5,j))-1))...
        -len(1);
    lam_p(5,j) = ((len(6)-len(1))/... #2 vs 999
        ((prist(1,j)/prist(6,j))-1))...
        -len(1);
    lam_p(6,j) = ((len(3)-len(2))/... #25 vs 41
        ((prist(2,j)/prist(3,j))-1))...
        -len(2);
    lam_p(7,j) = ((len(4)-len(2))/... #25 vs 204
        ((prist(2,j)/prist(4,j))-1))...
        -len(2);
    lam_p(8,j) = ((len(5)-len(2))/... #25 vs 500
        ((prist(2,j)/prist(5,j))-1))...
        -len(2);
    lam_p(9,j) = ((len(6)-len(2))/... #25 vs 999
        ((prist(2,j)/prist(6,j))-1))...
        -len(2);
    lam_p(10,j) = ((len(4)-len(3))/... #41 vs 204
        ((prist(3,j)/prist(4,j))-1))...
        -len(3);
    lam_p(11,j) = ((len(5)-len(3))/... #41 vs 500
        ((prist(3,j)/prist(5,j))-1))...
        -len(3);
    lam_p(12,j) = ((len(6)-len(3))/... #41 vs 999
        ((prist(3,j)/prist(6,j))-1))...
        -len(3);
    lam_p(13,j) = ((len(5)-len(4))/... %204 vs 500
        ((prist(4,j)/prist(5,j))-1))...
        -len(4);
    lam_p(14,j) = ((len(6)-len(4))/... %204 vs 999
        ((prist(4,j)/prist(6,j))-1))...
        -len(4);
    lam_p(15,j) = ((len(6)-len(5))/... %500 vs 999
        ((prist(5,j)/prist(6,j))-1))...
        -len(5);
    
    %%%%% ints
    
    lam_i(1,j) = ((len(2)-len(1))/... #2 vs 25
        ((ints(1,j)/ints(2,j))-1))...
        -len(1);
    lam_i(2,j) = ((len(3)-len(1))/... #2 vs 41
        ((ints(1,j)/ints(3,j))-1))...
        -len(1);
    lam_i(3,j) = ((len(4)-len(1))/... #2 vs 204
        ((ints(1,j)/ints(4,j))-1))...
        -len(1);
    lam_i(4,j) = ((len(5)-len(1))/... #2 vs 500
        ((ints(1,j)/ints(5,j))-1))...
        -len(1);
    lam_i(5,j) = ((len(6)-len(1))/... #2 vs 999
        ((ints(1,j)/ints(6,j))-1))...
        -len(1);
    lam_i(6,j) = ((len(3)-len(2))/... #25 vs 41
        ((ints(2,j)/ints(3,j))-1))...
        -len(2);
    lam_i(7,j) = ((len(4)-len(2))/... #25 vs 204
        ((ints(2,j)/ints(4,j))-1))...
        -len(2);
    lam_i(8,j) = ((len(5)-len(2))/... #25 vs 500
        ((ints(2,j)/ints(5,j))-1))...
        -len(2);
    lam_i(9,j) = ((len(6)-len(2))/... #25 vs 999
        ((ints(2,j)/ints(6,j))-1))...
        -len(2);
    lam_i(10,j) = ((len(4)-len(3))/... #41 vs 204
        ((ints(3,j)/ints(4,j))-1))...
        -len(3);
    lam_i(11,j) = ((len(5)-len(3))/... #41 vs 500
        ((ints(3,j)/ints(5,j))-1))...
        -len(3);
    lam_i(12,j) = ((len(6)-len(3))/... #41 vs 999
        ((ints(3,j)/ints(6,j))-1))...
        -len(3);
    lam_i(13,j) = ((len(5)-len(4))/... #204 vs 500
        ((ints(4,j)/ints(5,j))-1))...
        -len(4);
    lam_i(14,j) = ((len(6)-len(4))/... #204 vs 999
        ((ints(4,j)/ints(6,j))-1))...
        -len(4);
    lam_i(15,j) = ((len(6)-len(5))/... #500 vs 999
        ((ints(5,j)/ints(6,j))-1))...
        -len(5);
    
end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
om_fit1 = zeros(1,6521);
om_fit2 = zeros(1,6521);
for i = 1:6521
    om_fit1(i) = 1/(OM(i)^2);
    om_fit2(i) = 1/(OM(i)^4);
end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dwin=0.3e12*(1.9*pi); % THz change this value to adjust the convolution width
% win=round(dwin/dom);
% g = gausswin_my(win); % <-- this value determines the width of the smoothing window
% g = g/sum(g);
% lam1(3500:6521)=conv(lam1(3500:6521),g,'same');

% dwin=0.3e12*(1*pi); % THz change this value to adjust the convolution width
% win=round(dwin/dom);
% g = gausswin_my(win); % <-- this value determines the width of the smoothing window
% g = g/sum(g);
% lam2(2501:6521)=conv(lam2(2501:6521),g,'same');


lam_p(:,3935:6521) = 0;
lam_i(:,3950:6521) = 0;

iMFP = lam_i(2,:); %2/41
iMFP(1:626) = lam_i(11,1:626); %41/500
pMFP = lam_p(5,:); %2/500
pMFP(3000:5000) = lam_p(2,3000:5000); %2/41

pMFP(1:125) = 0;
iMFP(1:125) = 0;

% 
% %%%%%%%%% PLOT ALL FITS %%%%%%%%%%%%%%%%%%
% figure;
% ax1 = subplot(3,5,1);ax2 = subplot(3,5,2);ax3 = subplot(3,5,3);ax4 = subplot(3,5,4);
% ax5 = subplot(3,5,5);ax6 = subplot(3,5,6);ax7 = subplot(3,5,7);ax8 = subplot(3,5,8);
% ax9 = subplot(3,5,9);ax10 = subplot(3,5,10);ax11 = subplot(3,5,11);ax12 = subplot(3,5,12);
% ax13 = subplot(3,5,13);ax14 = subplot(3,5,14);ax15 = subplot(3,5,15);
% 
% semilogy(ax1,OM, 1e9*lam_p(1,:),'Color',[1 0 0]);
% hold on;
% semilogy(ax2,OM, 1e9*lam_p(2,:),'Color',[1 0.25 0]);
% semilogy(ax3,OM, 1e9*lam_p(3,:),'Color',[1 0.5 0]);
% semilogy(ax4,OM, 1e9*lam_p(4,:),'Color',[1 0.75 0]);
% semilogy(ax5,OM, 1e9*lam_p(5,:),'Color',[0.75 0.75 0]);
% semilogy(ax6,OM, 1e9*lam_p(6,:),'Color',[0.5 0.75 0]);
% semilogy(ax7,OM, 1e9*lam_p(7,:),'Color',[0.25 0.75 0]);
% semilogy(ax8,OM, 1e9*lam_p(8,:),'Color',[0 0.75 0]);
% semilogy(ax9,OM, 1e9*lam_p(9,:),'Color',[0 0.75 0.25]);
% semilogy(ax10,OM, 1e9*lam_p(10,:),'Color',[0 0.75 0.5]);
% semilogy(ax11,OM, 1e9*lam_p(11,:),'Color',[0 0.75 0.75]);
% semilogy(ax12,OM, 1e9*lam_p(12,:),'Color',[0 0.5 0.75]);
% semilogy(ax13,OM, 1e9*lam_p(13,:),'Color',[0 0.25 0.75]);
% semilogy(ax14,OM, 1e9*lam_p(14,:),'Color',[0 0 0.75]);
% semilogy(ax15,OM, 1e9*lam_p(15,:),'Color',[0.25 0 0.75]);
% legend(ax1,'2/25');
% legend(ax2,'2/41');
% legend(ax3,'2/204');
% legend(ax4,'2/500');
% legend(ax5,'2/1000');
% legend(ax6,'25/41');
% legend(ax7,'25/204');
% legend(ax8,'25/500');
% legend(ax9,'25/1000');
% legend(ax10,'41/204');
% legend(ax11,'41/500');
% legend(ax12,'41/1000');
% legend(ax13,'204/500');
% legend(ax14,'204/1000');
% legend(ax15,'500/1000');
% ylabel(ax1,'MFP, nm');
% ylabel(ax6,'MFP, nm');
% ylabel(ax11,'MFP, nm');
% xlabel(ax11,'Frequency, THz');
% xlabel(ax12,'Frequency, THz');
% xlabel(ax13,'Frequency, THz');
% xlabel(ax14,'Frequency, THz');
% xlabel(ax15,'Frequency, THz');

% % PLOT MFPS %%%%%%%%
for i = 1:6251
    if iMFP(i) <= 0
        iMFP(i) = NaN;
    end
    if pMFP(i) <= 0
        pMFP(i) = NaN;
    end
end


plot(OM,1e9*pMFP,'r+','MarkerSize',1);
hold on;
plot(OM,1e9*iMFP,'+-','Color',[0.2 0.2 0.8],'MarkerSize',1);
plot(OM,om_fit1*7e3,'r:');
plot(OM,om_fit2*1.5e3,':','Color',[0.2 0.2 0.8]);
lgd = legend('(Bulk Si)_{prist.}',...
    '(Bulk Si)_{1.56% ints.}','\omega^{-2}','\omega^{-4}');
lgd.FontSize = 12;
axis([0 16 -500 1e4]);
xlabel('Frequency (THz)','fontsize',16);
ylabel('Mean Free Path, nm','fontsize',16);
title('Phonon Mean Free Paths','fontweight','bold','fontsize',18);

% %%%% PLOT VS VITALYS %%%%%%
% semilogy(omBTE,BTE,'+k','MarkerSize',4);
% hold on;
% semilogy(omDFT,DFT,'+','MarkerSize',4,'Color',[0.7 0 0.7]);
% semilogy(OM,1e9*pMFP,'r+-','MarkerSize',4,'LineWidth',4);
% 
% 
% lgd = legend('BTE-Tersoff','DFT-AlmaBTE','NEMD');
% lgd.FontSize = 16;
% axis([0 16 0.2 5e4]);
% xlabel('Frequency (THz)','fontsize',16);
% ylabel('Mean Free Path, nm','fontsize',16);
% title('Phonon Mean Free Paths','fontweight','bold','fontsize',18);
% clear i i2 i204 i25 i41 i500 j lam_i lam_p lamstr len lgd OM p2 ...
%     p204 p25 p41 p500 dom 
