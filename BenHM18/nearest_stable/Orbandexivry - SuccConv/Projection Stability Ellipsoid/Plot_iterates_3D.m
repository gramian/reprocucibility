function [] = Plot_iterates_3D()
%Plot_iterates_3D enables the view of the iterates Bi obtained using
%StableMain.m.  It plots the path followed by the iterates in the space of
%the coefficients in 3D. As the space of coefficients has as many
%dimensions as it has coefficients, the three dimensions are chosen as
%follow :
%        - The origin is given by the solution to the problem.
%        - The vector going from the origin to the original (unstable)
%          problem gives the X axis.
%        - The vector going from the origin to the starting point of the
%          algorithm gives the Y axis
%        - The Z axis is obtained by doing a Principal Component Analysis
%          on the iterates. The PCA identifies the dimensions were the
%          largest changes have occured.
% 
%As the axis aren't orthonormal most of the time, an orthonormalisation
%process is achieved.
%
%The iterates that are used were automatically saved (if the option
%compute was activated) as columns in the variables Bmat in stock.mat.
%Copyright : Francois-Xavier Orban de Xivry 
%        for Universite Catholique de Louvain
%
%Creation : 07/04/2010
%Modifications : 01/09/2010 : No modif., the file was checked for mistakes,
%                             but none could be found.
%                05/01/2011 : Making of the header.
%Last modified : 05/01/2011

%% preprocessing
%loading the variable containing the iterates
load stock.mat Bmat

%q is the square of the dimension of each iterate, k is the number of
%iterates (including original instable matrix and starting point)
[q,k] = size(Bmat);

%n is the size of each matrix iterate
n = sqrt(q);

%coord will contain the projection of the position in the n²-dimensional
%space in the 3-D space
coord = zeros(k,3);

%in order to take the best 3-D space you need somme reference points
A0 = Bmat(:,1);%the original problem
A1 = Bmat(:,2);%the starting point 
Ak = Bmat(:,k);%the solution 

%the principal direction out of the original, the starting and the solution
%point will be the 4th reference point (=Z axis)
Bmat(:,[1 2 k]) = [];
Ai = FindPCA(Bmat);
% Ai = Bmat(:,13) % if you want to try a special column as the Zaxis

%setting the reference points in their original n*n form.
A0 = reshape(A0,n,n);
A1 = reshape(A1,n,n);
Ak = reshape(Ak,n,n);
Ai = reshape(Ai,n,n);

%% processing

%setting the origin at the solution
O = Ak;
%the direction of the X, Y and Z axis
OX = (A0-O)/norm(A0-O,'fro');
OY = (A1-O)/norm(A1-O,'fro');
OZ = (Ai-O)/norm(Ai-O,'fro');

%the angles between them (it is needed since the coordinate system is not
%orthogonal)
wxy = acos(real(trace(OX'*OY)));
wyz = acos(real(trace(OY'*OZ)));
wzx = acos(real(trace(OZ'*OX)));

%construction of the symmetric matrix U, it is used to solve the projection
%on each axis
U = [1 cos(wxy) cos(wzx); cos(wxy) 1 cos(wyz); cos(wzx) cos(wyz) 1];

%iterating on each of the iterates except the original, the starting and
%the solution point (hence the k-3)
for j = 1 : k-3
    
    C = reshape(Bmat(:,j),n,n);%current matrix iterate to be projected
    C = C-O;%shifting it at the origin of the coordinate system
    
    %construction of a small 3*3 system to find the projection of the
    %matrix C on the 3D space <OX,OY,OZ>
    v = [real(trace(C'*OX)); real(trace(C'*OY)); real(trace(C'*OZ))];
    
    %the coordinates in <X,Y,Z> basis is given by (index in coord is 2+j
    %since your already now the two first one which are the original and
    %the starting point
    coord(2+j,:) = (U\v).';
end

%% orthonormalisation 

%coordinates in the original coordinate system
coord(1,:) = [norm(A0-O,'fro') 0 0];%corresponding to A0 
coord(2,:) = [0 norm(A1-O,'fro') 0];%corresponding to A1
coord(k,:) = [0 0 norm(Ak-O,'fro')];%corresponding to Ak(it is always [0 0 0] if the origin is Ak)

%coordinates in an orthonormal coordinate system
coordp = [coord(:,1) zeros(k,2)];%X axis doesn't change
coordp = coordp + [coord(:,2)*cos(wxy) coord(:,2)*sin(wxy) zeros(k,1)];%Y axis stays in the YX plane but is made orthogonal to X

%determining the orthogonal projection of OZ in the OXY plane:
%This is tricky since you have to use OX, OYo (o for orthogonal) as a basis
OYo = OY - cos(wxy)*OX;
OYo = OYo/norm(OYo,'fro');
wyoz = acos(real(trace(OYo'*OZ)));
ZPOXY = cos(wzx)*OX+ OYo*cos(wyoz);%projection of OZ in the OXY plane
wzPxy_z = acos(real(trace(ZPOXY'*OZ)));%angle between OZ and ZPOXY

coordp = coordp + [coord(:,3)*cos(wzx) coord(:,3)*cos(wyoz) coord(:,3)*sin(wzPxy_z)];


%% postprocessing

%Plots 
close all
    %first plot : showing the iterates in 3D no correction of the angles
    %(the base is NOT orthonormal)
    figure
    h = plot3(coord(3:(k-1),1),coord(3:(k-1),2),coord(3:(k-1),3),'.')
    hold on;
    h = plot3(coord([1],1),coord([1],2),coord([1],3),'.r');%original problem
    plot3(coord([2],1),coord([2],2),coord([2],3),'.m');%Starting point
    plot3(coord([k],1),coord([k],2),coord([k],3),'.c');%solution
    h = plot3(0,0,0,'o','markersize',10);
    hold off;

    %second plot : showing the iterates in a 3D orthonormal basis
    figure
    hold on;
    h = plot3(coordp([1],1),coordp([1],2),coordp([1],3),'.r');
    plot3(coordp([2],1),coordp([2],2),coordp([2],3),'.m');
    plot3(coordp([k],1),coordp([k],2),coordp([k],3),'.c');
    h = plot3(0,0,0,'o','markersize',10);
    h = plot3(coordp(3:(k-1),1),coordp(3:(k-1),2),coordp(3:(k-1),3),'.')
    hold off;
    
end

%% Principal component anylysis
function [Ai] = FindPCA(Bmat)
%[Ai] = FindPCA(Bmat) is used to compute the Principal Component Analysis
%of the iterates in Plot_iterates_3D.

%computing the size of the original matrices
[q,k] = size(Bmat);
n =sqrt(q);


[U,S,V] = svd(Bmat);

%the matrix returned corresponds to the principal direction linked to the
%max singular value (the singular values are ordered in decreasing order by
%the function SVD(X))
Ai = reshape(U(1,:),n,n);
end
    
