% clear;clc;close all;

% rawdata = load ('../work/2006-05/2006-05-05/2006-05-05.l1c');

% rawdata = load ('../../osu7539/output/bias/eigen_5p4cEIGEN/results/2003-05-22/2003-05-22.mvf.l1c');
rawdata = load ('../../osu7539/output/L1C_my/results/2003-05-22/2003-05-22.l1c');


t0 = rawdata(1,1);
t = (rawdata(:,1) - t0) / 86400.0;
lat = rawdata(:,2);
lon = rawdata(:,3);
l1c = rawdata(:,8);
csr = rawdata(:,9);
gfz = rawdata(:,10);
jpl = rawdata(:,11);


figure;plot(t, l1c, '.',t, csr,'r');


