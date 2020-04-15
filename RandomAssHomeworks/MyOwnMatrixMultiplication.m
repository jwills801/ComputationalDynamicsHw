function [Myt,Matlabt] = MyOwnMatrixMultiplication(n,m,p)

A = rand(n,m);
B = rand(m,p);
C = zeros(n,p);

[Arows,Acolumns] = size(A);
[Brows,Bcolumns] = size(B);

tic
for i = 1:Arows
    for j = 1:Bcolumns
        for k = 1:Acolumns
            C(i,j) = C(i,j)+A(i,k)*B(k,j);
        end
        
        
        
    end
end
Myt = toc;

tic
A*B;
Matlabt = toc;