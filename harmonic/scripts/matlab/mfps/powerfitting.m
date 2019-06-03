% best fit a*x^b

% prist, a = 5.98e-6, b = -1.8
% fit from 3 THz to 14 THz

% ints, a= 3.27e-6, b= -3.807
% fit from 1.2 THz to 14 THz

% for 2 to 5, a=2.7565e-11, b = -3.83


load('MFP-FinalVersion');

lam = transpose(pMFP); %change to math the mfp function 
OM = transpose(OM);
lam(3974:6521) = NaN;

[fit1,r1] = fit(OM(500:1250),lam(500:1250),'power1'); %2 to 5
[fit2,r2] = fit(OM(500:2500),lam(500:2500),'power1'); %2 to 10
[fit3,r3] = fit(OM(500:3750),lam(500:3750),'power1'); %2 to 15
[fit4,r4] = fit(OM(1250:2500),lam(1250:2500),'power1'); %5 to 10
[fit5,r5] = fit(OM(1250:3750),lam(1250:3750),'power1'); %5 to 15
[fit6,r6] = fit(OM(2500:3750),lam(2500:3750),'power1'); %10 to 15

[fit7,r7] = fit(OM(500:3500),lam(500:3500),'power1'); %1.2 to 14

% figure;
% ax1 = subplot(2,3,1);
% ax2 = subplot(2,3,2);
% ax3 = subplot(2,3,3);
% ax4 = subplot(2,3,4);
% ax5 = subplot(2,3,5);
% ax6 = subplot(2,3,6);

a = 9.08e-7;
b = 0.2669;

ex = zeros(6521);
for i = 1:6521
    ex(i) = OM(i)^-b;
end
ex = ex*a;

lam(1:554) = NaN;
semilogy(OM,lam);
hold on;
plot(fit7);
% axis([0 20 10 1e4]);
legend('MFP, Prist','\omega^{-4.659}')