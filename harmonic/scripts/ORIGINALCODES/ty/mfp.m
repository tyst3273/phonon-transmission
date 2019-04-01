%clear all;
load('membranes.mat');

tr1 = data(2,:); 
l1 = 1.1e-9; %1.1 nm

tr2 = data_l2(2,:);
l2 = 108.8e-9; %108.8 nm


dom = 1.570796326794897e+10;
dwin=0.3e12*(3*pi); % THz change this value to adjust the convolution width
win=round(dwin/dom);
g = gausswin_my(win); % <-- this value determines the width of the smoothing window
g = g/sum(g);
tr1=conv(tr1,g,'same');
tr2=conv(tr2,g,'same');
 
% for i = 1:size(tr1,2)
%     if tr1(i) < 0
%         tr1(i) = 0;
%     end
%     if tr2(i) < 0
%         tr2(i) = 0;
%     end
% end

lam = zeros(1,50000); %mfp

for i = 1:size(data,2)
    lam(i) = ((l2-l1)/((tr1(i)/tr2(i))-1))-l1;
end

% dwin=0.3e12*(3*pi); % THz change this value to adjust the convolution width
% win=round(dwin/dom);
% g = gausswin_my(win); % <-- this value determines the width of the smoothing window
% g = g/sum(g);
% lam=conv(lam,g,'same');

save('mfp.mat');

% plot(OM, 1e9*lam);
semilogy(OM, 1e9*lam, ref_om, ref);
legend('ty','ref');
axis([0 16 1 10e5]);
% plot(OM, tr1, OM, tr2);
