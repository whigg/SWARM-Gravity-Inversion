clear;clc;close all;

% d = '../work/2006-05/2006-05-05/2006-05-05';
d = '../../osu7539/output/bias/eigen_5p4cEIGEN/results/2011-10-23/2011-10-23';
rawdata = load ([d, '-FIT-KBRon.sim']);
t0 = rawdata(1,1);
tr = (rawdata(:,1) - t0) / 86400.0;
emv1 = rawdata(:,23);
emv2 = rawdata(:,27);
epw1 = rawdata(:,24);
epw2 = rawdata(:,28);
resmvr = rawdata(:,10);
respwr = rawdata(:,11);
resflr = rawdata(:,12);
csrr = rawdata(:,16);
gfzr = rawdata(:,17);
jplr = rawdata(:,18);

%ind = emv1 < 0.005 & emv2 < 0.005 & epw1 < 0.005 & epw2 < 0.005;
%ind = emv1 < 0.005 & emv2 < 0.005;
%ind = emv1 < 0.003; 
%ind = emv1 < 0.005 & emv2 < 0.005;
ind =  resmvr-csrr < 0.005 & resmvr-csrr > -0.005 & emv1 < 0.01;
%ind =  emv1 < 0.005 & resmvr-csrr < 0.005 & resmvr-csrr > -0.005; 


simdata = rawdata(:,:);

t = (simdata(:,1) - t0) / 86400.0;
c = (simdata(:,1) - t0) / 5400.0;
l1c = simdata(:,8);
resmv = simdata(:,10);
respw = simdata(:,11);
resfl = simdata(:,12);
csr = simdata(:,16);
gfz = simdata(:,17);
jpl = simdata(:,18);

figure;plot(t, l1c-mean(l1c));

figure;plot(t, resmv, '.',t, csr,'r');

errdata = rawdata(~ind,:);
te = (errdata(:,1) - t0) / 86400.0;
ce = (errdata(:,1) - t0) / 5400.0;
l1ce = errdata(:,8);
resmve = errdata(:,10);
respwe = errdata(:,11);
resfle = errdata(:,12);
csre = errdata(:,16);
gfze = errdata(:,17);
jple = errdata(:,18);
figure;plot(tr, resmvr, '.',te, resmve,'or',tr, csrr,'r');



ressm = smooth(resmv,20);
figure;plot(t,ressm, t,csr, t, jpl, t, gfz);
legend('res','csr','jpl','gfz');
corrcoef(jpl,csr)
corrcoef(gfz,csr)
corrcoef(ressm,csr)
size(t)/size(emv1)
