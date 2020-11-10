function []= Convergit(mat)

if nargin == 0
   S = load('stock.mat','Bmat');
   mat = getfield(S,'Bmat');
end

[q,k] = size(mat);
n = sqrt(q);
normit = zeros(1,k);

origin = reshape(mat(:,1),n,n);
for i = 1 : k

    temp = reshape(mat(:,i),n,n);
    normit(i) = norm(temp-origin,'fro');
end

err = normit(2:k-1)-ones(1,k-2)*normit(k);

figure;
semilogy(1:k-2,err(1:k-2),'-o')