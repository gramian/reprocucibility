function [] = plot_ellipse(A)
N = 20;


[U D V] = svd(A);

nloop = size(D,1)-2;

for i = 1:nloop
    % generate the ellipsoid at (0,0,0)
    %----------------------------------
    a = 1/sqrt(D(i,i));
    b = 1/sqrt(D(i+1,i+1));
    c = 1/sqrt(D(i+2,i+2));
    [X,Y,Z] = ellipsoid(0,0,0,a,b,c,N);

    %  rotate and center the ellipsoid to the actual center point
    %------------------------------------------------------------
    XX = zeros(N+1,N+1);
    YY = zeros(N+1,N+1);
    ZZ = zeros(N+1,N+1);
    for k = 1:length(X),
        for j = 1:length(Y),
            point = [X(k,j) Y(k,j) Z(k,j)]';
            P = V(i:i+2,i:i+2) * point;
            XX(k,j) = P(1);
            YY(k,j) = P(2);
            ZZ(k,j) = P(3);
        end
    end
    figure;
    mesh(XX,YY,ZZ);
    axis equal
    hidden off
end