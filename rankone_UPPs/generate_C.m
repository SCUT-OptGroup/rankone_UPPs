function C = generate_C(Q,c,a,q,n)

if q == 1
    
    c1 = c{1};
    
    if isempty(c1(c1>0))
        
        C{1} = Q{1} + diag(a{1}*ones(n,1)/n);
        
    else
        
        Q1 = Q{1};
        
        C{1} = [a{1} 0.5*c1';0.5*c1 Q1];
        
    end
    
else
    
    C{1} = [a{1} 0.5*c{1}' zeros(1,n);0.5*c{1} Q{1} zeros(n);zeros(n,1) zeros(n) zeros(n)];
    
    C{2} = [a{2} zeros(1,n) 0.5*c{2}';zeros(n,1) zeros(n) zeros(n);0.5*c{2} zeros(n) Q{2}];
   
end