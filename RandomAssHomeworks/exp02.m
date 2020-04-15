function [y,i] = exp02(x)

%%check if x<0
flag_xIsNeg = 0;
if x < 0 
    flag_xIsNeg = 1;
    x = -x;
end



%%

% y = exp(x) = exp(n+z)


n = floor(x);
z = x-n;



e = 2.718281828459046;
en = 1;
for i = 1:n
    en = en*e;
end



y = 1;


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
% this gave exp(z)
y = y*en;


if flag_xIsNeg == 1
    y = 1/y
end
