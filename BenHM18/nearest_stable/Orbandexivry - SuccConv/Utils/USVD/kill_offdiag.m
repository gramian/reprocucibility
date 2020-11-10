function [R] = kill_offdiag(B)

i = 1;
j = 2;
sign = -1;

if B(j,j) == 0
    c = 0;
    s = 1;
else
    t = sign * B(i,j)/B(j,j);
    c = 1/sqrt(1+t^2);
    s = t*c;
end

R = [c s; -s c];