function [varargout] = PlotStableSet_colored(A,B,B0,isdiscrete,varargin)
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
% History  : 05/01/2011 : Making of the header
%            30/03/2011 : Enabling the handling of polynomials
%            04/11/2011 : Enabling the option 'isdiscrete' + return of a
%                         optional output which contains the best solution
%                         found.
%            08/11/2011 : Adding the axis equal command to the graphs +
%                         enabling the use of variables from the workspace
%                         if the arguments are not defined.
%            12/07/2012 : Added the possibility of drawing circles at the
%                         end of the file.
%            13/07/2012 : Cleaning of the file, change of names, paramters
%                         at the beginning of the file, interpolation
%                         finally done properly (not addiotonal lines
%                         appears once the contour is drawn.
%            24/08/2012 : Added e_lon as a harcoded parameter
%            23/09/2012 : Enables the possibility of representing the
%                         projection of additional points in the plane.
% Last modified : 23/09/2012

%% Parameters of the method
nptx = 30; %number of samples used for the axis B0-B
npty = 30; %number of samples used for the axis B-A

%coefficient range of the search (B is at center, B0 on the right, A at the
%bottom)
leftrange   = 1.5; % leftlimit = leftrange *  -(B0-B)/norm(B0-B);
rightrange  = 1.5; %rightlimit =  rightrange * (B0-B)/norm(B0-B);
toprange    = 1; % toplimit = toprange* -(A-B)/norm(A-B);
bottomrange = 1; %bottomlimite = bottomrange * (A-B)/norm(A-B);

%coefficient for the numberof points on which we carry the interpolation
%number of integration points = [coeffx*nptx coeffy*npty]
coeffx = 3;
coeffy = 3;

%% varargin check
if nargin == 0;
    display('Using the variables from the "base" workspace')
    try
        [A, B, B0] = evalin('base', ['A', 'B', 'B0']);
    catch
        error('Inputs variables A, B or B0 does not exist in the workspace')
    end
end
if nargin < 4
    isdiscrete = 0;
end

[nl,nc] = size(A);
if any([nl,nc] - size(B)) || any([nl,nc] - size(B0))
    error('Some inputs have different sizes')
end

%% Initialization 

%default values
namevec = {'A','B','B0'};
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

%preprocessing for polynomials
if nl < nc
    A = A.';
    B = B.';
    B0= B0.';
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
e_lon = 0;
fprintf('-> Value of the relative additive perturbation\n   of the stability boundary epsilon : %d\n',e_lon)

%% Process : creation of the Basis with 3 points
X = B0-B;
Y = B-A;
D1 = X/norm(X,normtype);
D2 = Y/norm(Y,normtype);

%angle between D1 and D2
if ispoly
    theta = acos(real(sum(D2'*D1)));%angle between D1 and D2
else
    theta = acos(real(trace(D2'*D1)));%angle between D1 and D2
end
fprintf('-> Relative angle between the axis : %d\n',theta)

%Orthonormalization
D2 = D2- norm(D2,normtype)*cos(theta)*D1;
D2 = D2/norm(D2,normtype);

%limits of the search
leftlimit   = min(-norm(X,normtype), norm(Y,normtype)*cos(theta))*leftrange;
rightlimit  = max( norm(X,normtype), norm(Y,normtype)*cos(theta))*rightrange;
bottomlimit = -norm(Y,normtype)*sin(theta)*bottomrange;
toplimit    =  norm(Y,normtype)*sin(theta)*toprange;


%creating variables 
xmap = linspace(leftlimit,rightlimit,nptx);
ymap = linspace(bottomlimit,toplimit,npty);
coord = int8(zeros(nptx,npty));
border= int8(zeros(nptx,npty));
objvalD = zeros(nptx,npty);
% coord = [];

%Brute force search of stable points within the limits 
condsatisf = 0;
minnorm = 1e100;
bestSol = [];
for i = 1 : nptx
    if mod(i,10)== 0
        %         i
    end
    sol1 = B + xmap(i)*D1;
    for j = 1 : npty
        sol = sol1 + ymap(j)*D2;

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
            coord(i,j) = 1;
            normAX = norm(sol-A,normtype);
            objvalD(i,j) = normAX;
            %saving the indices of the best solution
            if normAX < minnorm
                minnorm = normAX;
                bestSol = [i j];
            end
            %special values for the borders in coord
            if xmap(i) == rightlimit || xmap(i) == leftlimit
                coord(i,j) = 2;
            end
            if ymap(j) == toplimit || ymap(j) == bottomlimit
                coord(i,j) = 3;
            end
            %             coord = [coord; xmap(i) ymap(j)];

            condsatisf = 0;%the variable has to be put back to 0 for the next loop
        else
            %do nothing right now
        end
    end
end

%Reconstructing the best solution found
if ~isempty(bestSol)
    BestSol = B + xmap(bestSol(1))*D1 + ymap(bestSol(2))*D2;
end
if nargout == 1
    varargout = {BestSol};
end

%processing the borders
BorderCoeff = 0;
for i =1 : nptx
    for j = 1 : npty
        if coord(i,j)~= 0
            %Parts of the limits that are stable
            if (i == nptx || i ==1) ||...
                    (j==npty || j==1)
                border(i,j) = 1;
                objvalD(i,j) = BorderCoeff*minnorm;
            end
        else
            %boundary of the stable set within the limits.
            %Here we set the "boundary" just outside the stable set so that
            %we do not destroy interesting data (this boundary is rightly
            %interesting)
            if (i>2) && (i<(nptx-1))
                %regarde à gauche et à droite
                if (coord(i-1,j)==1)||(coord(i+1,j)==1)
                    %the points just outside the domain must be set to a
                    %small value for the interpolation to be well performed
                    %by griddata.
                    objvalD(i,j) = BorderCoeff*minnorm; 
                    border(i,j) = -1;
                end
            end
            if (j>2) && (j<(npty-1))
                %regarde en bas et en haut
                if (coord(i,j-1)==1)||(coord(i,j-1)==1)
                    %the points just outside the domain must be set to a
                    %small value for the interpolation to be well performed
                    %by griddata.
                    objvalD(i,j) = BorderCoeff*minnorm;
                    border(i,j) = -1;
                end
            end
        end
    end
end

if ~isempty(coord)
%% Computing the orthonormalized coordinates
    %Computing the coordinates of the points that are stable or that forms
    %the boundary (inside the limits) and will be plotted 
    [I,J] = find((coord(:,:) > 0) |(border(:,:) == -1));
    u = ((I-1)*(rightlimit-leftlimit)/(nptx-1))+ leftlimit;
    v = ((J-1)*(toplimit-bottomlimit)/(npty-1))+ bottomlimit;
    x = u;
    y = v;

    %computing the coordinates of the points that are on the borders either
    %lateral or above and below
    [K,L] = find(border>0);
    k = ((K-1)*(rightlimit-leftlimit)/(nptx-1))+ leftlimit;
    l = ((L-1)*(toplimit-bottomlimit)/(npty-1))+ bottomlimit;
    
    % Doing the same for the border inside the domain
    [O,P]=find(border<0);
    o = ((O-1)*(rightlimit-leftlimit)/(nptx-1))+ leftlimit;
    p = ((P-1)*(toplimit-bottomlimit)/(npty-1))+ bottomlimit;

    
    %Coordinates of the original points
    coordB  = [0 0];
    coordB0 = [norm(X,normtype) 0];
    coordA  = [-norm(Y,normtype)*cos(theta) -norm(Y,normtype)*sin(theta)] ;
    
    
%% Compute additional points to the graphs if the option 'Addpoints' is used

    %At this stage we only compute the coordinates of the projections in
    %the plane defined by D1 and D2
    if addpoints == 1
        %creating an orthonormalized base from D1(or x) to simplify the
        %projection 
        D2ortog = D2 - real(trace(D2'*D1))*D1;
        D2orth = D2ortog/norm(D2ortog,'fro');
        
        %collecting data and instaciating appropriate sized vector of
        %coordinates
        numAddpoints = size(scoord, 2 );
        Addcoord = zeros(numAddpoints,2);
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
            % WARNING : it is different that the computation of x and y
            % above because we have to work on an orhtornomalized base to
            % get the correct coordinates of the projection of the point in
            % the plane.
            curval = curval-B;%So that we set the origin at B
            proj1 = real(trace(D1'*curval));
            proj2 = real(trace(D2orth'*curval));           
            Addcoord(i,:) = [proj1 proj2]; 
            
            %computing the distance to the plane
            pcoord = proj1*D1+proj2*D2orth;
            vnormal = curval-pcoord;
            vndist = norm(vnormal,'fro');
            fprintf('-> Distance from the point %s to the plane : %4.2d\n',Addname{i},vndist);
        end
    end
%% Check before the figure

    if isempty(x)&& isempty(k)
        warning('Figures cannot be drawn, no additional stable points detected around the 3 points')
    else
%% First figure : Stable (non-orthonormalized base)

        figure;
        hold on;
        plot(u,v,'.','MarkerSize',4)
        plot(0,0,'.m','MarkerSize',6)
        plot(norm(X,normtype),0,'.m','MarkerSize',6)
        plot(norm(Y,normtype)*cos(theta),-norm(Y,normtype)*sin(theta),'.m','MarkerSize',6)
        hold off;

%% second fig : Stable set (orthonormalized base)
    figure;
    hold on;
    %Adding the stable points, and the limit of the search
    plot(x,y,'.','MarkerSize',4)
    plot(k,l,'.','MarkerSize',6,'MarkerEdgecolor',[0.5 0.2 0.2])
    
    %Adding the original points
    plot(coordB(1),coordB(2),'.m','MarkerSize',6)
    plot(coordB0(1),coordB0(2),'.m','MarkerSize',6)
    plot(coordA(1),coordA(2),'.m','MarkerSize',6)

    %computing the offset before putting the label of the points
    [xlim] = get(gca,'Xlim');
    xoffs = (xlim(2)-xlim(1))/50;
    [ylim] = get(gca,'Ylim');
    yoffs = (ylim(2)-ylim(1))/30;
    
    %adding the text associated to them
    text(coordB(1),coordB(2),nameB,'color','m','Verticalalignment','bottom','horizontalalignment','right')
    text(coordB0(1),coordB0(2),nameB0,'color','m','Verticalalignment','bottom','horizontalalignment','right')
    text(coordA(1),coordA(2),nameA,'color','m','Verticalalignment','bottom','horizontalalignment','right')
    
    %adding the additional points if the option is activated
    if addpoints == 1
        for i = 1 : numAddpoints
            cox = Addcoord(i,1);
            coy = Addcoord(i,2);
            namep = Addname{i};
            plot(cox,coy,'.m','Markersize',6);
            text(cox,coy,namep,'color','m','Verticalalignment','bottom','horizontalalignment','right');
        end
    end

    axis equal

    hold off;

%% Third figure : level curves
    %     clf
    figure;
    hold on;

    %keeping the distances of the stable points only
    val  = objvalD(sub2ind(size(objvalD),I,J));

    %creating the uniform grid
    ax = linspace(xlim(1),xlim(2),coeffx*nptx);
    ay = linspace(ylim(1),ylim(2),coeffy*npty);
    [axm,aym]= meshgrid(ax,ay);

    try
        %interpolating the values of the stable distances on a uniform grid
        Z = griddata(x,y,val,axm,aym,'linear');       


        %---computing and adding the contour---
        %noticable level that we want to have vmin,B,B0 and vmax
        vmin = minnorm;
        vmax = max(val);
        vb0 = norm(A-B0,normtype);
        vb  = norm(A-B,normtype);
        vl = min(vb,vb0);
        vu = max(vb,vb0);

        %a ratio is computed so that the number of curves is spread evenly
        %outside of the 4 noticable points.
        rap1 = ceil(((vl-vmin)/(vmax-vmin))*10);
        rap2 = ceil(((vu-vl)/(vmax-vmin))*10);
        rap3 = ceil(((vmax-vu)/(vmax-vmin))*10);
        if rap1 < eps, rap1 =1; end;
        if rap2 < eps, rap2 =1; end;
        if rap3 < eps, rap3 =1; end;

        %tracing the contours
        v = [[vmin:(vl-vmin)/rap1:vl],[vl+(vu-vl)/rap2:(vu-vl)/rap2:vu-(vu-vl)/rap2],[vu:(vmax-vu)/rap3:vmax]];
        [C,h]=contour(axm,aym,Z,v);
        
%         %If you want to see where are the boundaries of the stable set
%         bdout = plot(k,l,'.','markersize',4);
%         bdin = plot(o,p,'.','markersize',4);

        %adding the original points
        plot(coordB(1),coordB(2),'.m','MarkerSize',6)
        plot(coordB0(1),coordB0(2),'.m','MarkerSize',6)
        plot(coordA(1),coordA(2),'.m','MarkerSize',6)

        %computing the offset before putting the label of the points
        [xlim] = get(gca,'Xlim');
        xoffs2 = (xlim(2)-xlim(1))/50;
        [ylim] = get(gca,'Ylim');
        yoffs2 = (ylim(2)-ylim(1))/30;
    
        %adding the text associated to them
        text(coordB(1),coordB(2),nameB,'color','m','Verticalalignment','bottom','horizontalalignment','right')
        text(coordB0(1),coordB0(2),nameB0,'color','m','Verticalalignment','bottom','horizontalalignment','right')
        text(coordA(1),coordA(2),nameA,'color','m','Verticalalignment','bottom','horizontalalignment','right')

        %adding the additional points if the option is activated
        if addpoints == 1
            for i = 1 : numAddpoints
                cox = Addcoord(i,1);
                coy = Addcoord(i,2);
                namep = Addname{i};
                plot(cox,coy,'.m','Markersize',6);
                text(cox,coy,namep,'color','m','Verticalalignment','bottom','horizontalalignment','right');
            end
        end
    
        axis equal

        hold off;

    catch
        %impossible to plot the contours, close the figure, issue a warning
        close(gcf);
        warning('PSS:NotEnoughData','The contour could not be drawn\n This is probably due to a lack of stable points in our search area.')
        rethrow(lasterror)
    end
    
    %Drawing an additional circle on the last plot, the utility can be to
    %check that the angle is well computed.
    additionalcircle = 0;
    if additionalcircle
        display('Displaying additional circle. \n Deactivation inside the file is possible : look for ''additionalcircle''\n')

        %the radius is userdefined as well as the spread of the angle
        center = coordA;%coordinate (X,Y) in a line vector [X Y]
        radius = -0.05-center(2);%scalar
        angleinterval = [pi/4 3*pi/4];%in radian

        %process
        hold on;
        angle = linspace(angleinterval(1),angleinterval(2),500);
        angle = reshape(angle,[],1);
        circlepoints = ones(size(angle))*center + radius*[cos(angle) sin(angle)];
        plot(circlepoints(:,1),circlepoints(:,2))
        hold off;
    end
    end

else
    display('-> No stable value for alpha and beta')
end



