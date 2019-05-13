
dT= 80;
grad_T= 0.00652/1e-10;
A= 1065e-20;

data=prist(6,:);
kb=1.3806e-23;


kom=data.*(kb*dT)/A/grad_T*1e12;

kappa=trapz(OM,kom);

