clear;clc;close all;
rawdata = load ('l1c.asc');
t0 = rawdata(1,1);

%emv1 = rawdata(:,23);
%emv2 = rawdata(:,27);
%epw1 = rawdata(:,24);
%epw2 = rawdata(:,28);

%ind = emv1 < 0.004 & emv2 < 0.004 & epw1 < 0.004 & epw2 < 0.004;
%ind = emv1 < 0.004 & emv2 < 0.004;

%simdata = rawdata(ind,:);
simdata = rawdata;

t = (simdata(:,1) - t0) / 86400.0;
c = (simdata(:,1) - t0) / 5400.0;
resmv = simdata(:,8);
csr = simdata(:,9);
gfz = simdata(:,10);
jpl = simdata(:,11);

rmv = smooth(resmv,20);

per = length(t)/17280.0;


j2c=corrcoef(jpl,csr);
g2c=corrcoef(gfz,csr);
g2j=corrcoef(gfz,jpl);
l2g=corrcoef(rmv,gfz);
l2c=corrcoef(rmv,csr);
l2j=corrcoef(rmv,jpl);
fileID = fopen('corr.txt','w');
fprintf(fileID,'%.2f%% JPL2CSR %f GFZ2CSR %f GFZ2JPL %f L1C2CSR %f L1C2JPL %f L1C2GFZ %f\n',...
    per*100, j2c(1,2),g2c(1,2),g2j(1,2), l2c(1,2), l2j(1,2), l2g(1,2));
fclose(fileID);
