
x = [0:0.1:4.5];
y = zeros(size(x));
for i = 1 : length(x)
y(i) = EvalFmu(x(i),Betas,C1,alpha);
end

plot(x,y)