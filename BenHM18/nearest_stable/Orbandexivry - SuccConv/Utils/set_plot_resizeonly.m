function [] = set_plot_resizeonly(varargin)
%set_plot_label(width,height,object,type,centered)modifies the size of the
%figure so that it has the required size. The difference with set_plot.m is
%that it resize the figure so that it includes(or not)the labels and title
%of the figure. It does not change the fontsize or linewidth.
%The size can either be an overall size including the labels and title or
%just the size of the graph. 
%The paper property of figure is then set for an easy .pdf impression
%All inputs are optional
%     width : the width of the object (cm) if omitted width and height will
%             be kept at their current value.
%     height: the height of the object (cm) if omitted width and height will
%             be kept at their current value.
%     object: the object on which to apply the dimensions:'figure' or
%             'axes'.
%     type  : 'tight' if you dont want labels and title to be included (in
%              which case object will systematically be 'axes')
%             'extended' if you want the labels and title to be included
%     centered : 'centered' if you want to have you graph centered inside
%                the window. Having a centered graph when 'figure' is
%                selected can lead to some weird dimension of the axes

%default
height = 0;
width = 0;
type = 'extended';
object = 'figure';
centeredflag = 0;
add_width = 0;

%input processing
for i = 1 : nargin-1
    if isnumeric(varargin{i})&& isnumeric(varargin{i+1})
        width = varargin{1};
        height  = varargin{2};
    end
end
for i = 1 : nargin
    if ischar(varargin{i})
        switch lower(varargin{i})
            case {'tight','extended'}
                type = lower(varargin{i});
            case {'figure','axes'}
                object = lower(varargin{i});
            case 'centered'
                centeredflag = 1;
            otherwise
                error('The input string "%s" does not exist',varargin(i))
        end
    end
end

if strcmp(type,'tight')
    object = 'axes';
end

%permet de modifier les unit�s avec lesquelles on mesure les tailles

set(gca,'units','centimeters');
set(gcf,'units','centimeters');

posf = get(gcf,'position');
posa = get(gca,'position');
pos = posa;
if width == 0
    switch object
        case 'axes'
            width = posa(3);
            height = posa(4);
        case 'figure'
            width = posf(3);
            height = posf(4);
    end
end

switch type
    case 'tight'
        add_space = 0.2;
        in = zeros(1,4);
        correction = [0 0];%because tightinset is not always satisfactory
    case 'extended'
        in = get(gca,'tightinset');
        add_space = 0;
        correction = [0 0.2];%because tightinset is not always satisfactory
        if centeredflag
            add_width = in(1)-in(3);
        end
end

switch object
    case 'axes'
        pos(3:4) = [width height];
        posf(3:4) = in(1:2)+correction +pos(3:4)+ in(3:4)+add_space*ones(1,2)+[add_width 0]; %[0 0.1] added because the thight inset doesn't let enough space for the labels below       
    case 'figure'
        posf(3:4) = [width height];
        pos(3:4) = [width height]-in(1:2)-in(3:4)-correction-[add_width 0];
end

pos(1:2) = in(1:2)+correction + add_space/2*ones(1,2);%correction added because the thight inset doesn't let enough space for the labels below

%setting the new dimensions in effect
set(gcf,'position',posf);
set(gca,'position',pos);

%preprocessing for the printing
%The next three lines set the papersize the exact dimension of the graph so
%it does not appear on an A4 paper when printed as a .pdf
set(gcf,'paperunits','centimeters')

newposf = get(gcf,'position');
papersize = [newposf(3) newposf(4)];

% add_width = 1.1*fontsize;
% if centered_flag
%     papersize(1) = papersize(1) + add_width;
% end

set(gcf,'papersize',papersize)
set(gcf,'paperunits','default')

%remettre tout en place
set(gcf,'units','normalized')
set(gca,'units','normalized')
set(gcf,'paperpositionmode','auto')%->this is actually very important it adapts automaticaly the position and size of the paper to the size on the screen
set(gca,'activepositionproperty','outerposition')