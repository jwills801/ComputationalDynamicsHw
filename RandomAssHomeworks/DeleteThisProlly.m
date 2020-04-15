tic
HowBig = 1000;
step = 100;

n = 1:step:HowBig;
m = 1:step:HowBig;
p = 1:step:HowBig;


Myt = zeros(1,HowBig);
Matlabt = zeros(1,HowBig);
Attempt = zeros(1,HowBig);
indexer = 1;

for i = 1:step:HowBig
    
    [Myt(indexer),Matlabt(indexer)] = MyOwnMatrixMultiplication(n(indexer),m(indexer),p(indexer));
    Attempt(indexer) = n(indexer);
    indexer = indexer + 1;
end

plot(Attempt,Myt,Attempt,Matlabt)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
legend('My naive function speed','Matlabs function speed')
xlabel('size of matrices')
ylabel('time (sec)')
title('My multiplication versus Matlabs')
    toc
