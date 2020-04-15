clc
clear
clear all
u = .5;

% ln(u) = (u-1) - (u-1)^2/2 + (u-1)^3/3 - (u-1)^4/4 + (u-1)^5/5 + ....

%%check if u=0
if u == 0 
    b = q;
else 
    counter = 1;
    c = u-1;
end
    
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
    
    
    
    q =5;
    b = q + c;
    
