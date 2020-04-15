
n = 100000;
detvals = nan(1,n);

for k = 1:n

A = nan(6,6);
R = rand(6,6);


for i = 1:6
    for j = 1:6
if R(i,j) > 0.5
    A(i,j) = 1;
else
    A(i,j) = -1;
end
    end
end


detvals(k) = det(A);
end

domain = 1:1:n;
plot(domain,detvals)
xlabel('iteration number')
ylabel('determinant value')
title('6X6 max determinant value')
