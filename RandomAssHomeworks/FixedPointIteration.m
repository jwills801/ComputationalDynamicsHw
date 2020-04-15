function root = FixedPointIteration(g,x0)


tol = 1e-8;
max_iter = 300;
flag = 0;
i = 0;

xprev = g(x0);

while flag == 0
    i = i +1;
    
    x = g(xprev);
    err = abs(x-xprev);
    xprev = x;
    
    if err < tol
        flag = 1;
    end
    
    if i > max_iter
        flag = -1;
    end
end
end
    