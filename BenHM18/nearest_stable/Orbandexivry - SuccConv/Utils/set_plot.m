function [] = set_plot(width,height,centered_flag,show_axe,linewidth,fontsize,markersize,faxe)
%set_plot(width,height,centered_flag,show_axe,linewidth,fontsize,markersize,faxe)
%modify the dimension of the axes so that the graph has the right
%dimension,the right graph properties and is well positioned inside its
%frame.  
%
%All inputs are optional.
%   Input : width         : the width of the graph(white area) if 0 the width
%                           will not be changed.
%           height        : the height of the graph(white area) if 0 the
%                           height will not be changed.
%           centered_flag : if 1 margin will be added on the right so that
%                           the graph lies in the center of the window.
%           show_axe      : if 0 won't display the ticks and labels of the
%                           graph and only the content will be displayed
%                           (no space around). Default value :1. If it is
%                           activated the size of the figure, not the axes,
%                           will be of width*height. It is automatically
%                           centered
%           linewidth     : the linewidth of the lines in the graph (if
%                           any). Will not be changed if 0.
%           fontsize      : the fontsize to be applied on the text inside the
%                           graph and on the axes as well. Will not be
%                           changed if 0.
%           markersize    : the size of the markers inside the graph (if
%                           any).Will not be changed if 0.
%           faxe          : a valid handle on which we apply set_plot
%
% History : ../../2009 : Creation
%           28/08/2012 : Adding the option show_axe
% Last modified : 28/08/2012
%% Preprocessing
if nargin < 2
    width = 0;
    height = 0;
end
if nargin < 3
    centered_flag = 0;
end
if nargin < 5
    linewidth = 0;
end
if nargin < 6
    fontsize = 0;
end
if nargin < 7
    markersize = 0;
end
if nargin < 8
    faxe = gca;
end
if nargin < 4
    visibility = get(gca,'Visible')
    if strcmp(visibility,'on')
        show_axe = 1;
    else
        show_axe = 0;
    end
end
if nargin > 8
    error('Wrong number of inputs argument')
end


%setting the current axes
if nargin >=8
    if ~ishandle(faxe)
        error('The given handle is not a valid handle')
    else
        axes(faxe)%gca is now equal to faxe
    end
end 

%check center_flag is valid
if isnumeric(centered_flag) == 0
    centered_flag = lower(centered_flag);
    if strcmp(centered_flag,'centered') == 1
        centered_flag = 1;
    else
        centered_flag = 0;
    end
end

doLine = 1;
doFont = 1;
doMarker = 1;

if linewidth ==0
    doLine = 0;
end
if fontsize == 0
    doFont = 0;
end
if markersize == 0
    doMarker = 0;
end

%-end preprocessing
%% measurements
%%%
set(gcf,'Units','centimeters');
set(gca,'Units','centimeters');

%% linewidth, fontsize and marker size
if doFont && show_axe
    set(gca,'Fontsize',fontsize);
    set(get(gca,'Xlabel'),'Fontsize',fontsize);
    set(get(gca,'Ylabel'),'Fontsize',fontsize);
    set(get(gca,'Zlabel'),'Fontsize',fontsize);
end

% axis normal;

%%%
%Changing the linewidth and the fontsize of the elements inside the axes.

hh= get(gca,'children');%éléments dans les axes

if length(hh) ==1

    %plot de la convergence
    if doLine
        set(hh(1),'LineWidth',linewidth);
    end
    if doMarker
        set(hh(1),'Markersize',markersize);
    end

else
    %this only works for some graphs
    type = get(hh,'type');
    m = length(hh);
    for i =1:m
        if strcmp(type(i),'text')
            if doFont
                set(hh(i),'Fontsize',fontsize)
            end
        else
            if doLine
                set(hh(i),'linewidth',linewidth)
            end
            if doMarker
                set(hh(i),'Markersize',markersize)
            end
        end
    end
end

%% Dimensions
add_space = 0;%dimensions in centimeter to add so that nothing is cut when 
%the axes are not displayed
if ~show_axe
    add_space = 0.1;
end
    
%Setting the dimensions at the required width and height.
%pos = [start_hor_pos start_vert_pos width height]
pos = get(gca,'position');
if width == 0;
    width = pos(3);%the size of the axes should not change if no size was specified
else
    width = width-add_space;
end
if height == 0;
    height = pos(4);
else
    height = height-add_space;
end
%tightinset contains dimensions to add to pos to contain all things inside
%the figure, its configuration is also pos = [ideal_start_hor_pos ideal_start_vert_pos
%add_width add_height] where add_witdh and add_height are increments to add
%to width and height while ideal_start_hor_pos and ..._vert_pos replace the
%ones obtained using pos (its really not a good conception).
%initialization that will keep only the content of the axes
if show_axe
    inset = get(gca,'tightinset');
else
    inset = add_space*ones(1,4) ;%element 1 and 2 are used in 'position' of gca
                                 %element 3 and 4 are used in 'position' of
                                 %gcf
end

%I do not increment right away the width and height because it might need
%some correction once I modify the two first arguement of pos
newpos = [inset(1) inset(2) width height];
set(gca,'position',newpos);

%once your change your position the ticks might have change as well as the
%space they need, thus computing a correction
inset_corr = inset;
if show_axe
    inset_corr = get(gca,'tightinset');
end
if any(inset_corr ~= inset)%does not enter if show_axe = 0
    %so we correct the space needed
    newpos = [inset_corr(1) inset_corr(2) width height];
    inset = inset_corr;
    set(gca,'position',newpos);
end

%now we look at the space needed on the right and on the top of the graph
poso = get(gca,'outerposition');

newposo  = [poso(1:2) newpos(1:2)+newpos(3:4)+inset(3:4)];

% set(gca,'outerposition',newposo) doesn't do it properly instead we use
% the position property of the FIGURE not the axes;
posf(3:4) = newposo(3:4);
set(gcf,'position',posf);


%% Appearance
%display the axes
if show_axe
    axis on
else
    axis off;
end

%this part set the dimension of the box so that the graph part is centered
%(this optioin is usefull if you have a lot in your margin so that the
%graph look quite out of the frame whereas it should be at its center)

if centered_flag == 1
    %puts the same distance on the right than on the left.
    add_width = inset(1) - inset(3);
    posf(3) = posf(3)+add_width;
    set(gcf,'position',posf);
end

%% Print preprocessing
%puts back the normal properties so everything in matlab works good
%and center the figures on the screen

set(gcf,'units','normalized')
posf = get(gcf,'position');
posf(1) = (1-posf(3))/2;% coordinates for centering the graph
posf(2) = (1-posf(4))/2;

set(gcf,'position',posf)%center it

%default values
set(gcf,'paperpositionmode','auto')
set(gcf,'units','default')
set(gca,'units','default')

%%%

%preprocessing for the printing
%The next three lines set the papersize the exact dimension of the graph so
%it does not appear on an A4 paper when printed as a .pdf
set(gcf,'paperunits','centimeters')
papersize = [newposo(3)+0.2 newposo(4)+0.2];%the 0.2 is an experimental value so that the graph is not cut even slightly

if centered_flag
    papersize(1) = papersize(1) + add_width;
end

set(gcf,'papersize',papersize)
set(gcf,'paperunits','default')
