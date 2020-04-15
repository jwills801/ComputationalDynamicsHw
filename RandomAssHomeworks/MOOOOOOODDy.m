

Re = 4000:10000:1e8;
RR = [.05 .04 .03 .02 .015 .01 .005 .002 .001 5e-4 2e-4 1e-4 5e-5 1e-5 5e-6 1e-6]; 

hold on
max_iter1 = length(RR);
flag1 = 0;
j = 0;
while flag1 == 0;
    j = j + 1;
  
    

max_iter = length(Re);
flag = 0;
i = 0;
ff = nan(size(Re));

while flag == 0
    i = i +1 ;
    g = @(x) (-2*log10(RR(j)/3.7+2.7/(Re(i)*sqrt(x)))).^-2;
    ff(i) = MyOwnFixedPointIteration(g,.05);
    if i >= max_iter
        flag = -1;
    end
end



plot(Re,ff)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Reynolds Number')
ylabel('friction factor')
legend('.05','.04','.03','.02','.015','.01','.005','.002','.001','5e-4','2e-4','1e-4','5e-5','1e-5','5e-6','1e-6')

if j >= max_iter1
        flag1 = -1;
end
    
end
hold off
 