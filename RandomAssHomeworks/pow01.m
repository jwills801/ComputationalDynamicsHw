function p = pow01(x,y)

%returns  p = x^y
%...I hope....




% b = ln(x)

%x^y = exp(y*ln(x))

% x^y = exp(w)
% w = y*ln(x) = y*b
% OK! Lets find b!!!


%%
% b = ln(x) = ln(exp(q)*u) = q + ln(u)
% c = ln(u)

% b = q + c;


q = floor(x);

% u = x/exp(q)
% d = exp(q)
% u = x/d

%% WE NEED d = EXP(q)
%%check if q<0
flag_qIsNeg = 0;
if q < 0 
    flag_qIsNeg = 1;
    q = -q;
end





% y = exp(w) = exp(n+z)


n = floor(q);
z = q-n;



e = 2.718281828459046;
en = 1;
for i = 1:n
    en = en*e;
end



d = 1;


i = 0;

flag = 0;

tol = 1e-4;

max_iter = 300;

newterm = 1;




while flag == 0
    i = i+1;
    
    newterm = newterm*z/i;
    
    d = d + newterm;
    
    rel_err = abs(newterm);
    
    if rel_err <= tol
        flag = 1;
        
        
    elseif i >= max_iter
        flag = -1;
        
    end
end
% this gave exp(w)
d = d*en;


if flag_qIsNeg == 1
    d = 1/d;
end

% OKAY! WE FOUND d AND IT WORKS









%% NOW WE ARE GOING TO FIND THE LN(u)
% c = ln(u)

u = x/d;

% ln(u) = (u-1) - (u-1)^2/2 + (u-1)^3/3 - (u-1)^4/4 + (u-1)^5/5 + ....

%%check if u=0
if u == 0 
    b = q;
end
    counter = 1;
    c = u-1;

    
    MaxIterations = 300;
    Tolerance = 1e-5;
    Flag = 0;
    
    
    while Flag == 0
        
        counter = counter + 1;
        
        NewCTerm = (-1)^(counter-1)*((u-1)^counter)/counter;
        
        c = c + NewCTerm;
    
        error = abs(NewCTerm);
        
    if counter >= MaxIterations
       Flag = -1;
    end
        
    if error < Tolerance
        Flag = 1;
    end
        
    end
    
 % OKAY! C IS NOW THE LN(u) AS LONG AS u<0   
    
 
 
 
 
 
 %%
    b = q + c;




% b = ln(x)

%x^y = exp(y*ln(x))

% x^y = exp(w)
% w = y*ln(x) = y*b
w = y*b;


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
p = p*en;


if flag_wIsNeg == 1
    p = 1/p;
end

%Ucheck = x/exp(q)
%u

%Bcheck = log(x)
%b

%Wcheck = y*log(x)
%w

%Pcheck = exp(y*log(x))
%p = exp(w)
    