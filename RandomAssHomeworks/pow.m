function p = pow(x,y)

%returns  z = x^y
%...I hope....




% p = ln(x)

%x^y = exp(y*ln(x))

% x^y = exp(w)
% w = y*ln(x) = y*p
% OK! Lets find p!!!

% ln(x+1) = x - x^2/2 + x^3/3 - x^4/4 + x^5/5













% p = ln(x)

%x^y = exp(y*ln(x))

% x^y = exp(w)
% w = y*ln(x) = y*p
w = y*p;


%%

%%check if w<0
flag_wIsNeg = 0;
if w < 0 
    flag_wIsNeg = 1;
    w = -w;
end





% y = exp(w) = exp(n+z)


n = floor(w);
z = w-n;



e = 2.718281828459046;
en = 1;
for i = 1:n
    en = en*e;
end



p = 1;


i = 0;

flag = 0;

tol = 1e-4;

max_iter = 300;

newterm = 1;




while flag == 0
    i = i+1;
    
    newterm = newterm*z/i;
    
    p = p + newterm;
    
    rel_err = abs(newterm);
    
    if rel_err <= tol
        flag = 1;
        
        
    elseif i >= max_iter
        flag = -1;
        
    end
end
% this gave exp(w)
p = p*en


if flag_wIsNeg == 1
    p = 1/p
end

    