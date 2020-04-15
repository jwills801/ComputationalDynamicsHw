

function [y,i] = exp01(x)
% y = exp(x)

y = 1;


%%

i = 0;

flag = 0;

tol = 1e-4;

max_iter = 300;

newterm = 1;




while flag == 0
    i = i+1;
    
    newterm = newterm*x/i;
    
    y = y + newterm;
    
    rel_err = abs(newterm);
    
    if rel_err <= tol
        flag = 1;
        
        
    elseif i >= max_iter
        flag = -1;
        
    end
end
