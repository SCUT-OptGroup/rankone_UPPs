%% ***************************************************************
%%
%% dcSNCG to the following UBPPs:
%
%   min_{x_i\in{-1,1}^n} \Pi_{i=1}^q(<x_i,Q_ix_i>+<c_i,x_i>+a_i)
%
%%
%% ***************************************************************
function [fobj,xsol,infeas,time,rho] = dcSNCG_start(Q,c,a,q)

%% Q, c and a are cells where Q{i}=Q_i, c{i}=c_i and a{i}=a_i, i=1,...,q.
%% q = 1 or q = 2.

randstate = 100;
randn('state',double(randstate));
rand('state',double(randstate));

OPTIONS_XMM.gaptol = 1e-8;

OPTIONS_XMM.objtol = 1e-5;

OPTIONS_XMM.printyes = 1;

OPTIONS_XMM.maxiter = 10000;

n = size(Q{1},1);

if q == 1
    
    tempC = generate_C(Q,c,a,q,n);
    
    C0 = tempC{1};
    
    p = size(C0,1);
    
    seig = eigs(C0,1,'smallestreal');
    
    minus_seig = -min(seig,0);
    
    C = C0 + minus_seig*eye(p);  %% C is necessarily PSD
    
    CFnorm = norm(C,'fro');
    
    OPTIONS_XMM.lowgam = 1e-2;
    
    OPTIONS_XMM.uppgam = 1e-3;
    
    OPTIONS_XMM.CFnorm = CFnorm;
    
    if CFnorm>=1.0e+5
        
        Lip = 10;
        
    elseif CFnorm>=5.0e+3
        
        Lip = 1;
        
    else
        
        Lip = 0.1;
    end
    
    tstart = clock;
    
    rho = 1e-1;  sigma = 1.05;
    
    [xsol,rho] = DCP_QASNCG(C,OPTIONS_XMM,p,rho,sigma,Lip);
    
    time = etime(clock,tstart)
    
    infeas = max(abs((xsol.*xsol).^(1/2)-1))
    
    fobj = xsol'*C0*xsol
    
    
else
    
    tempC = generate_C(Q,c,a,q,n);
    
    C1 = tempC{1};
    
    C2 = tempC{2};
    
    p = size(C1,1);
    
    maxFnorm = max(norm(C1,'fro'),norm(C2,'fro'));
    
    Lip = abs(eigs(C1,1))*norm(C2,'fro')+abs(eigs(C2,1))*norm(C1,'fro');
    
    OPTIONS_XMM.CFnorm = maxFnorm;
    
    OPTIONS_XMM.fingam = Lip;
    
    tstart = clock;
    
    rho = 1e-1;  sigma = 1.05;
    
    m = max(min(50,round(p/2)),2);
    
    tempV = randn(m,p);
    
    tempV_cnorm = sqrt(dot(tempV,tempV,1));
    
    V = bsxfun(@rdivide,tempV,tempV_cnorm);
    
    Xstart = V'*V;
    
    [xsol,rho] = DCP_PolySNCG(C1,C2,Xstart,OPTIONS_XMM,p,rho,sigma,Lip);
    
    time = etime(clock,tstart)
    
    infeas = max(abs((xsol.*xsol).^(1/2)-1))
    
    fobj = (xsol'*C1*xsol)*(xsol'*C2*xsol)

end
