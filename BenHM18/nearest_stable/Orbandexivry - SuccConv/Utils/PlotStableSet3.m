function [varargout] = PlotStableSet3(A,B,B0,C,isdiscrete,varargin)
%[BestSol] = PlotStableSet(A,B,B0,isdiscrete,'nameA','nameB','nameB0')
%plots the limits of the stable in a plane determined by the matrices A, B
%and B0. In the multidimensional space of coefficients, matrices can be
%seen as points and thus define a plane.
% Called as : PlotStableSet()
%             PlotStableSet(A,B,B0)
%             PlotStableSet(A,B,B0,isdiscrete)
%             PlotStableSet(A,B,B0,isdiscrete,'Addname',namevec)
%             PlotStableSet(A,B,B0,isdiscrete,'Addpoints',scoord)
%             PlotStableSet(A,B,B0,isdiscrete,'Addname',namevec,'Addpoints',scoord)
%
%The stable region is drawn by evaluating a succession of matrices given as
%linear combination of (B0-B) and (B-A) (which are normed). The number of
%evaluation is given by the parameters npta and nptb. The points to be
%evaluated are comprised in a interval given by [amin,amax] and
%[bmin,bmax].
%
% Input : A          : The unstable matrix
%         B          : The solution found by StableMain.m
%         B0         : The initial point for the algorithm
%         C          : An additional point which does not lie in the plane
%                      A-B-B0
%         isdiscrete : (optional) 1 if the problem is discrete, 0 if
%                      continuous. Default : 1.
%         namevec    : if the option 'Addname' is used, then namevec is a 3
%                      by 1 string vector which contains the names to write
%                      for A, B an B0 in the rigth order, default names are
%                      added otherwise.
%         scoord     : if the option 'Addpoints' is used, then scoord is a
%                      structure containing all the points to be added in a
%                      field 'value', if a field 'name' is added, then the
%                      strings contained in it are added on the graph
%
% Note : If using the option 'Addname', you must specify either all or none
% of the names of the points.
%
% Output : BestSol : Is an optional output which contains the best solution
%                    found for the problem in the analyzed region.
%
% History  : 14/10/2012 : Creation for PlotStableSet
% Last modified : 14/10/2012

%% Parameters of the method
nptx = 150; %number of samples used for the axis B0-B
npty = 150; %number of samples used for the axis B-A
nptz = 150; %number of samples used for the axis C-B

%coefficient range of the search (B is at center, B0 on the right, A at the
%bottom)
Xnegrange   = 1; % Xneglimit = Xnegrange *  -(B0-B)/norm(B0-B);
Xposrange  = 1; %Xposlimit =  Xposrange * (B0-B)/norm(B0-B);
Ynegrange    = 1; % Yposlimit = Ynegrange* -(A-B)/norm(A-B);
Yposrange = 1; %Yneglimite = Yposrange * (A-B)/norm(A-B);
Znegrange = 1; %Yneglimite = Znegrange * -(C-B)/norm(C-B);
Zposrange = 1; % Yposlimit = Zposrange* -(C-B)/norm(C-B);


%% varargin check
if nargin == 0;
    display('Using the variables from the "base" workspace')
    try
        A = evalin('base', 'A');
        B = evalin('base', 'B');
        B0= evalin('base', 'B0');
        C = evalin('base', 'C');
    catch
        error('Inputs variables A, B, B0 or C does not exist in the workspace')
    end
end
if nargin < 5
    isdiscrete = 0;
end

[nl,nc] = size(A);
if (any([nl,nc] - size(B)) || any([nl,nc] - size(B0)))|| any([nl,nc] - size(C))
    error('Some inputs have different sizes')
end

%% Initialization

%default values
namevec = {'A','B','B0','C'};
addpoints = 0;
lvar = length(varargin);
if mod(lvar,2)~=0
    error('Wrong number of optional arguments, they must come by pair')
end
numoption = lvar/2;
for i = 1 : numoption
    option = lower(varargin{2*i-1});
    switch option
        case 'addname'
            namevec = varargin{2*i};
        case 'addpoints'
            addpoints = 1;
            scoord = varargin{2*i};
        otherwise
            error('The option name is not valid');
    end
end

nameA = namevec{1};
nameB = namevec{2};
nameB0 = namevec{3};
nameC = namevec{4};

%preprocessing for polynomials
if nl < nc
    A = A.';
    B = B.';
    B0= B0.';
    C = C.';
    temp = nl;
    nl = nc;
    nc = temp;
end

if nl == nc
    %does not do anything we deal with a matrix
    normtype = 'fro';
    ispoly = 0;
elseif nc == 1 %We have A,B,B0 as vectors
    %the first index must be one
    A = A/A(1);
    B = B/B(1);
    B0 = B0/B0(1);
    C = C/C(1);
    normtype = 2; %option so we now we have polynomials
    ispoly = 1;
else
    error('-> Input "a" is not a vector or a square matrix');
end

%warning for the type of case
fprintf('\n')
if isdiscrete
    disp('PROCESSING A DISCRETE CASE')
else
    disp('PROCESSING A CONTINUOUS CASE')
end


%epsilon for checking the epsilon-stability of a point (relative to the
%usual boundary)
e_lon = 1e-12;
fprintf('-> Value of the relative additive perturbation\n   of the stability boundary epsilon : %d\n',e_lon)

%% Process : creation of the Basis with 3 points
X = [0 0 0 1]';%A-B;
Y = [0 0 1 0]';%B0-B;
Z = [0 1 0 0]';%C-B;
D1 = X/norm(X,normtype);
D2 = Y/norm(Y,normtype);
D3 = Z/norm(Z,normtype);

%Orthonormalization
D2ortog = D2 - real(trace(D2'*D1))*D1;
D2orth = D2ortog/norm(D2ortog,normtype);
D3ortog = D3 - real(trace(D3'*D1))*D1 -real(trace(D3'*D2orth))*D2orth;
D3orth = D3ortog/norm(D3ortog,normtype);


%angle between A,B0,C 
if ispoly
    theta = acos(real(sum(Y'*D1)));%angle between A and D1
    phi   = pi/2-acos(real(sum(Z'*D3orth))); % angle between C and the plane XY
    D3p = D3 - sin(phi)*D3orth;
    theta2 = acos(real(sum(D3p'*D1)));
else
    theta = acos(real(trace(D2'*D1)));%angle between D1 and D2
    phi   = pi/2-acos(real(trace(Z'*D3orth))); % angle between C and the plane XY
    D3p = D3 - sin(phi)*D3orth;
    theta2 = acos(real(trace(D3p'*D1)));%"theta" angle for C
end
fprintf('-> Relative angle between Y and D1 : %d\n',theta)
fprintf('-> Relative angle between Z and the plane D1D2 : %d\n',phi)
fprintf('-> Relative angle between the projection of Z and D1 : %d\n',theta2)

%limits of the search
Xneglimit =  -norm(X,normtype)*Xnegrange;%min(-norm(X,normtype), norm(Y,normtype)*cos(theta))*Xnegrange;
Xposlimit =  norm(X,normtype)*Xposrange;%max( norm(X,normtype), norm(Y,normtype)*cos(theta))*Xposrange;
Yneglimit = -norm(Y,normtype)*sin(theta)*Ynegrange;
Yposlimit =  norm(Y,normtype)*sin(theta)*Yposrange;
Zneglimit = -norm(Z,normtype)*sin(phi)*Znegrange;
Zposlimit =  norm(Z,normtype)*sin(phi)*Zposrange;


%creating variables
xmap = linspace(Xneglimit,Xposlimit,nptx);
ymap = linspace(Yneglimit,Yposlimit,npty);
zmap = linspace(Zneglimit,Zposlimit,nptz);

coord = int8(zeros(nptx,npty,nptz));
border= int8(zeros(nptx,npty,nptz));
objvalD = zeros(nptx,npty,nptz);
% coord = [];

%Brute force search of stable points within the limits
condsatisf = 0;
minnorm = norm(A-B,'fro');
bestSol = [];
for i = 1 : nptx
    if mod(i,10)== 0
        %         i
    end
    sol1 = B + xmap(i)*D1;
    for j = 1 : npty
        sol2 = sol1 + ymap(j)*D2orth;
        
        for k = 1:nptz
            sol = sol2 + zmap(k)*D3orth;
            %computing the roots of the solution
            if ispoly
                vp = roots(sol);
            else
                vp = eig(sol);
            end
            %verifying the the matrix satisfies its stability criterion
            if isdiscrete && all(abs(vp)<=1+e_lon)
                condsatisf = 1;
            elseif ~isdiscrete && all(vp <= 0+e_lon)
                condsatisf = 1;
            end
            
            %if the stability criterion is satisfy coord(i,j) is 1, 0 otherwise
            if  condsatisf
                coord(i,j,k) = 1;
                normAX = norm(sol-A,normtype);
                objvalD(i,j,k) = normAX;
                %saving the indices of the best solution
                if normAX < minnorm
                    minnorm = normAX;
                    bestSol = [i j k];
                end
%                 %special values for the borders in coord
%                 if xmap(i) == Xposlimit || xmap(i) == Xneglimit
%                     coord(i,j,k) = 2;
%                 end
%                 if ymap(j) == Yposlimit || ymap(j) == Yneglimit
%                     coord(i,j,k) = 3;
%                 end
%                 if zmap(k) == Zposlimit || zmap(k) == Zneglimit
%                     coord(i,j,k) = 4;
%                 end
                %             coord = [coord; xmap(i) ymap(j)];
                
                condsatisf = 0;%the variable has to be put back to 0 for the next loop
            else
                %do nothing right now
            end
        end
    end
end

%Reconstructing the best solution found
if ~isempty(bestSol)
    BestSol = B + xmap(bestSol(1))*D1 + ymap(bestSol(2))*D2orth + zmap(bestSol(3))*D3orth;
end
if nargout == 1
    varargout = {BestSol};
end

%processing the borders
BorderCoeff = 0;
for i =1 : nptx
    for j = 1 : npty
        for k = 1 : nptz
            if coord(i,j,k)~= 0
                if (i>2) && (i<(nptx-1))
                    %regarde à gauche et à droite
                    if (coord(i-1,j,k)==0)||(coord(i+1,j,k)==0)
                        %the points just outside the domain must be set to a
                        %small value for the interpolation to be well performed
                        %by griddata.
                        border(i,j,k) = -1;
                    end
                end
                if (j>2) && (j<(npty-1))
                    %regarde en bas et en haut
                    if (coord(i,j-1,k)==0)||(coord(i,j-1,k)==0)
                        %the points just outside the domain must be set to a
                        %small value for the interpolation to be well performed
                        %by griddata.
                        border(i,j,k) = -1;
                    end
                end
                if (k>2) && (k<(nptz-1))
                    %regarde en bas et en haut
                    if (coord(i,j,k-1)==0)||(coord(i,j,k+1)==0)
                        %the points just outside the domain must be set to a
                        %small value for the interpolation to be well performed
                        %by griddata.
                        border(i,j,k) = -1;
                    end
                end
            end
        end
    end
end

if ~isempty(coord)
    %% Computing the orthonormalized coordinates
    %Computing the coordinates of the points that are stable or that forms
    %the boundary (inside the limits) and will be plotted
    [Ind] = find(border == -1);%warning J is a linear index over the 2nd and 3rd dimensions
    [I,J,K] = ind2sub(size(coord),Ind);
    x = ((I-1)*(Xposlimit-Xneglimit)/(nptx-1))+ Xneglimit;
    y = ((J-1)*(Yposlimit-Yneglimit)/(npty-1))+ Yneglimit;
    z = ((K-1)*(Zposlimit-Zneglimit)/(nptz-1))+ Zneglimit;
    
%     %computing the coordinates of the points that are on the borders either
%     %lateral or above and below
%     [Ind] = find(border>0);
%     [L,M,N] = ind2sub(size(coord),Ind);
%     l = ((L-1)*(Xposlimit-Xneglimit)/(nptx-1))+ Xneglimit;
%     m = ((M-1)*(Yposlimit-Yneglimit)/(npty-1))+ Yneglimit;
%     n = ((N-1)*(Zposlimit-Zneglimit)/(nptz-1))+ Zneglimit;
    
%     % Doing the same for the border inside the domain
%     [Ind]=find(border<0);
%     [O,P,Q] = ind2sub(size(coord),Ind);
%     o = ((O-1)*(Xposlimit-Xneglimit)/(nptx-1))+ Xneglimit;
%     p = ((P-1)*(Yposlimit-Yneglimit)/(npty-1))+ Yneglimit;
%     q = ((Q-1)*(Zposlimit-Zneglimit)/(nptz-1))+ Zneglimit;
    
    
    %Coordinates of the original points
    coordB  = [0 0 0];
    coordA = [real(trace(D1'*AB)) real(trace(D2orth'*AB)) real(trace(D3orth'*AB))];
    coordB0  = [real(trace(D1'*B0B)) real(trace(D2orth'*B0B)) real(trace(D3orth'*B0B))];
    coordC  = [real(trace(D1'*CB)) real(trace(D2orth'*CB)) real(trace(D3orth'*CB))];
    
    
    %% Compute additional points to the graphs if the option 'Addpoints' is used
    
    %At this stage we only compute the coordinates of the projections in
    %the plane defined by D1 and D2
    if addpoints == 1
        
        %collecting data and instaciating appropriate sized vector of
        %coordinates
        numAddpoints = size(scoord, 2 );
        Addcoord = zeros(numAddpoints,3);
        Addname = cell(numAddpoints,1);
        
        %computing the coordinates
        for i = 1 :numAddpoints
            
            %current point
            curpoint = scoord(1,i);
            curval = curpoint.value;
            if isfield(scoord,'name')
                Addname{i} = curpoint.name;
            else
                Addname{i} = sprintf('P%i',i);
            end
            
            %coordinates
            curval = curval-B;%So that we set the origin at B
            proj1 = real(trace(D1'*curval));
            proj2 = real(trace(D2orth'*curval));
            proj3 = real(trace(D3orth'*curval));
            Addcoord(i,:) = [proj1 proj2 proj3];
            
            %computing the distance to the plane
            pcoord = proj1*D1+proj2*D2orth + proj3*D3orth;
            vnormal = curval-pcoord;
            vndist = norm(vnormal,'fro');
            fprintf('-> Distance from the point %s to the plane : %4.2d\n',Addname{i},vndist);
        end
    end
    %% Check before the figure
    
    if isempty(x)
        warning('Figures cannot be drawn, no additional stable points detected around the 3 points')
    else
        %% First fig : Stable set (orthonormalized base)
%         figure;
%         hold on;
        %Adding the stable points, and the limit of the search
        colored_points
        %plot3(x,y,z,'.','MarkerSize',4)
%         plot3(l,m,n,'.','MarkerSize',6,'MarkerEdgecolor',[0.5 0.2 0.2])
        
        %Adding the original points
        plot3(coordB(1),coordB(2),coordB(3),'.m','MarkerSize',6)
        plot3(coordB0(1),coordB0(2),coordB0(3),'.m','MarkerSize',6)
        plot3(coordA(1),coordA(2),coordA(3),'.m','MarkerSize',6)
        plot3(coordC(1),coordC(2),coordC(3),'.m','MarkerSize',6)
        
        %computing the offset before putting the label of the points(not
        %used anymore)
        [xlim] = get(gca,'Xlim');
        [ylim] = get(gca,'Ylim');
        [zlim] = get(gca,'Zlim');
        
        %adding the text associated to them
        text(coordB(1),coordB(2),coordB(3),nameB,'color','m','Verticalalignment','bottom','horizontalalignment','right')
        text(coordB0(1),coordB0(2),coordB0(3),nameB0,'color','m','Verticalalignment','bottom','horizontalalignment','right')
        text(coordA(1),coordA(2),coordA(3),nameA,'color','m','Verticalalignment','bottom','horizontalalignment','right')
        text(coordC(1),coordC(2),coordC(3),nameC,'color','m','Verticalalignment','bottom','horizontalalignment','right')
        
%         %adding the additional points if the option is activated
%         if addpoints == 1
%             for i = 1 : numAddpoints
%                 cox = Addcoord(i,1);
%                 coy = Addcoord(i,2);
%                 coz = Addcoord(i,3);
%                 namep = Addname{i};
%                 plot3(cox,coy,coz,'.m','Markersize',6);
%                 text(cox,coy,coz,namep,'color','m','Verticalalignment','bottom','horizontalalignment','right');
%             end
%         end
        
%         axis equal
        
        hold off;
%         
%         %% Third figure : level curves
%         %     clf
%         figure;
%         hold on;
%         
%         %keeping the distances of the stable points only
%         val  = objvalD(sub2ind(size(objvalD),I,J,K));
%         
%         %creating the uniform grid
%         ax = linspace(xlim(1),xlim(2),coeffx*nptx);
%         ay = linspace(ylim(1),ylim(2),coeffy*npty);
%         [axm,aym]= meshgrid(ax,ay);
%         
%         try
%             %interpolating the values of the stable distances on a uniform grid
%             Z = griddata(x,y,val,axm,aym,'linear');
%             
%             
%             %---computing and adding the contour---
%             %noticable level that we want to have vmin,B,B0 and vmax
%             vmin = minnorm;
%             vmax = max(val);
%             vb0 = norm(A-B0,normtype);
%             vb  = norm(A-B,normtype);
%             vl = min(vb,vb0);
%             vu = max(vb,vb0);
%             
%             %a ratio is computed so that the number of curves is spread evenly
%             %outside of the 4 noticable points.
%             rap1 = ceil(((vl-vmin)/(vmax-vmin))*10);
%             rap2 = ceil(((vu-vl)/(vmax-vmin))*10);
%             rap3 = ceil(((vmax-vu)/(vmax-vmin))*10);
%             if rap1 < eps, rap1 =1; end;
%             if rap2 < eps, rap2 =1; end;
%             if rap3 < eps, rap3 =1; end;
%             
%             %tracing the contours
%             v = [[vmin:(vl-vmin)/rap1:vl],[vl+(vu-vl)/rap2:(vu-vl)/rap2:vu-(vu-vl)/rap2],[vu:(vmax-vu)/rap3:vmax]];
%             [C,h]=contour(axm,aym,Z,v);
%             
%             %         %If you want to see where are the boundaries of the stable set
%             %         bdout = plot(k,l,'.','markersize',4);
%             %         bdin = plot(o,p,'.','markersize',4);
%             
%             %adding the original points
%             plot(coordB(1),coordB(2),'.m','MarkerSize',6)
%             plot(coordB0(1),coordB0(2),'.m','MarkerSize',6)
%             plot(coordA(1),coordA(2),'.m','MarkerSize',6)
%             
%             %computing the offset before putting the label of the points
%             [xlim] = get(gca,'Xlim');
%             xoffs2 = (xlim(2)-xlim(1))/50;
%             [ylim] = get(gca,'Ylim');
%             yoffs2 = (ylim(2)-ylim(1))/30;
%             
%             %adding the text associated to them
%             text(coordB(1),coordB(2),nameB,'color','m','Verticalalignment','bottom','horizontalalignment','right')
%             text(coordB0(1),coordB0(2),nameB0,'color','m','Verticalalignment','bottom','horizontalalignment','right')
%             text(coordA(1),coordA(2),nameA,'color','m','Verticalalignment','bottom','horizontalalignment','right')
%             
%             %adding the additional points if the option is activated
%             if addpoints == 1
%                 for i = 1 : numAddpoints
%                     cox = Addcoord(i,1);
%                     coy = Addcoord(i,2);
%                     namep = Addname{i};
%                     plot(cox,coy,'.m','Markersize',6);
%                     text(cox,coy,namep,'color','m','Verticalalignment','bottom','horizontalalignment','right');
%                 end
%             end
%             
%             axis equal
%             
%             hold off;
%             
%         catch
%             %impossible to plot the contours, close the figure, issue a warning
%             close(gcf);
%             warning('PSS:NotEnoughData','The contour could not be drawn\n This is probably due to a lack of stable points in our search area.')
%             rethrow(lasterror)
%         end
%         
%         %Drawing an additional circle on the last plot, the utility can be to
%         %check that the angle is well computed.
%         additionalcircle = 0;
%         if additionalcircle
%             display('Displaying additional circle. \n Deactivation inside the file is possible : look for ''additionalcircle''\n')
%             
%             %the radius is userdefined as well as the spread of the angle
%             center = coordA;%coordinate (X,Y) in a line vector [X Y]
%             radius = -0.05-center(2);%scalar
%             angleinterval = [pi/4 3*pi/4];%in radian
%             
%             %process
%             hold on;
%             angle = linspace(angleinterval(1),angleinterval(2),500);
%             angle = reshape(angle,[],1);
%             circlepoints = ones(size(angle))*center + radius*[cos(angle) sin(angle)];
%             plot(circlepoints(:,1),circlepoints(:,2))
%             hold off;
%         end
    end
    
else
    display('-> No stable value for alpha and beta')
end



