%Stable_param.m defines parameters are accessible to the user in order to
%tune the software at his convienence. Most regards post-processing but
%precision can be tuned using this script.

%----------------
%POST PROCESSING
%----------------

%noms de fichier à donner
Sol_filename = 'Sol';
Conv_filename = 'Convergence_rate';

%'Article' style
w_artcl = 5; %width [cm]
h_artcl = 5; % height [cm]
cent_artcl = ''; %figure is centered in the output file
showax_artcl = 0;%if 1 display the axes otherwise won't
ms_artcl = 5; %Marker size
fs_artcl = 12; %Font size
lw_artcl = 1.7; % Line width
ext_artcl = 'pdf';

%'Poster' style
w_poster = 12;
h_poster = 12;
cent_poster = '';
showax_poster = 0;
ms_poster = 12;
fs_poster = 20;
lw_poster = 3;
ext_poster = 'svg';

%'Personnalized' style (currently talk)
w_perso = 5;
h_perso = 5;
cent_perso = '';
showax_perso = 1;
ms_perso = 5;
fs_perso = 12;
lw_perso = 1.7;
ext_perso = 'pdf';


bleufonce = [0 38 116]/255;
lightblue = [74 129 212]/255;
degradebleu = [213 229 255; 170 204 255; 128 179 255; 103 164 255; ...
            85 153 255; 51 134 255; 42 127 255; 0 102 255;...
            0 85 212; 0 68 170]/255;
orangeflash = [254 143 25]/255;
black = [0 0 0];
%couleur des graphes !!Valid for the poster and personnalized mode only !!
% put [] for default colors.
axes_color = bleufonce;%black
iterate_color = [bleufonce;lightblue];%[]
startingpoint_color = bleufonce;
originalpoint_color = bleufonce;
converg_color = [];