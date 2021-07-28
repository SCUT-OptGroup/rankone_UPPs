%% ***************************************************************
% filename: example_Poly
%% ***************************************************************

randstate = 100;
randn('state',double(randstate));
rand('state',double(randstate));

state = rng;

kappa = 2*100;

Q{1} = randn(kappa/2);

rng(state);

Q{1} = (Q{1} + Q{1}')/2;

Q{1} = Q{1}/norm(Q{1},2);

q1 = randn(1);

Qe = sum(Q{1},2);

c{1} = 0.5*Qe;

a{1} = 0.25*sum(Qe) + q1;

Q{1} = 0.25*Q{1};

randstate=1000;
randn('state',double(randstate));
rand('state',double(randstate));

state = rng;

Q{2} = randn(kappa/2);

rng(state);

Q{2} = (Q{2} + Q{2}')/2;

Q{2} = Q{2}/norm(Q{2},2);

q2 = randn(1);

Qe = sum(Q{2},2);

c{2} = 0.5*Qe;

a{2} = 0.25*sum(Qe) + q2;

Q{2} = 0.25*Q{2};

n = size(Q{1},1);

q = 2;  m = max(min(50,round((q*n+1)/2)),2);

%% *********************** DCP_QAfactor ***************************

[fobj1,xsol1,infeas1,time1] = dcFAC_start(Q,c,a,q,m);

%% *********************** DCP_QAfactor ***************************

[fobj2,xsol2,infeas2,time2] = dcFACls_start(Q,c,a,q,m);


%% *********************** DCP_QASNCG ***************************

[fobj,xsol,infeas,time,rho] = dcSNCG_start(Q,c,a,q);

