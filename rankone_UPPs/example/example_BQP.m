%% ***************************************************************
% filename: example_BQP
%% ***************************************************************

clear all

load('G1');  

Q{1} = full(W);  n = size(W,1);

c{1} = zeros(n,1);  a{1} = 0;

q = 1;   m = 50;  

%% *********************** DCP_QAfactor ***************************

[fobj1,xsol1,infeas1,time1] = dcFAC_start(Q,c,a,q,m);


%% *********************** DCP_QASNCG ***************************

[fobj,xsol,infeas,time,rho] = dcSNCG_start(Q,c,a,q);


