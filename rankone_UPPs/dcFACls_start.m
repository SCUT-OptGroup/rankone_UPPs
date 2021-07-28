%% ***************************************************************
%%
%% dcFAC_ls to the following UBPPs:
%
%   min_{x_i\in{-1,1}^n} \Pi_{i=1}^q(<x_i,Q_ix_i>+<c_i,x_i>+a_i)
%
%%
%% ***************************************************************
function [fobj,xsol,infeas,time,rho] = dcFACls_start(Q,c,a,q,m)

%% Q, c and a are cells where Q{i}=Q_i, c{i}=c_i and a{i}=a_i, i=1,...,q.
%% q = 1 or q = 2.
%% m is the number of rows of the variable in factorized problem

randstate = 100;
randn('state',double(randstate));
rand('state',double(randstate));

OPTIONS_VMM.gaptol = 1e-8;

OPTIONS_VMM.objtol = 1e-8;

OPTIONS_VMM.printyes = 1;

OPTIONS_VMM.maxiter = 10000;

n = size(Q{1},1);

if q == 1
    
    tempC = generate_C(Q,c,a,q,n);
    
    C0 = tempC{1};
    
    p = size(C0,1);
    
    seig = eigs(C0,1,'smallestreal');
    
    minus_seig = -min(seig,0);
    
    C = C0 + minus_seig*eye(p);  %% C is necessarily PSD
    
    OPTIONS_VMM.CFnorm = norm(C,'fro');

    rho = 1e-3;  sigma = 1.005;
    
    tstart = clock;
    
    tempV = randn(m,p);
    
    tempV_cnorm = sqrt(dot(tempV,tempV,1));
    
    tempV = bsxfun(@rdivide,tempV,tempV_cnorm);
    
    V = tempV;  VC = V*C;
    
    [xsol,rho] = DCP_QAfactor_ls(C,V,VC,OPTIONS_VMM,rho,sigma);
    
    time = etime(clock,tstart)
    
    infeas = max(abs((xsol.*xsol).^(1/2)-1))
    
    fobj = xsol'*C0*xsol
    
    
else
    
    tempC = generate_C(Q,c,a,q,n);
    
    C1 = tempC{1};
    
    C2 = tempC{2};
    
    p = size(C1,1);
    
    maxFnorm = max(norm(C1,'fro'),norm(C2,'fro'));
    
    LL = abs(eigs(C1,1))*norm(C2,'fro')+abs(eigs(C2,1))*norm(C1,'fro');
    
    OPTIONS_VMM.CFnorm = maxFnorm;

    OPTIONS_VMM.Lip = 6*p*LL; 
    
    rho = 1e-3;  sigma = 1.005; 
    
    tstart = clock;
    
    tempV = randn(m,p);
    
    tempV_cnorm = sqrt(dot(tempV,tempV,1));
    
    V = bsxfun(@rdivide,tempV,tempV_cnorm);
    
    VC1 = V*C1; VC2 = V*C2;
    
    [xsol,rho] = DCP_PolyQAfac_ls(C1,C2,V,VC1,VC2,OPTIONS_VMM,rho,p,sigma);
    
    time = etime(clock,tstart)
    
    infeas = max(abs((xsol.*xsol).^(1/2)-1))
    
    fobj = (xsol'*C1*xsol)*(xsol'*C2*xsol)

end
