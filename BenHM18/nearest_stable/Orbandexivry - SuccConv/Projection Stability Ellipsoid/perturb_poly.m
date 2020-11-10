a = poly([-1+j -1+j -1+j -1+j -1+j])

m = length(a);
n =100;
r = zeros(m-1,n);
figure;
hold on;
b = a(m);
for i = 1: n
    a(m) = b*2*(i-1)/n;
    r(:,i) = roots(a);
    plot(real(r(:,i)),imag(r(:,i)),'o','markersize',4)
end

