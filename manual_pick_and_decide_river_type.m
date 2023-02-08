%  this script allows users to manually pick rivers from a DEM, and
%  visualize river profiles to manually decide river types.
% Yanyan Wang (wangyanyan0607@hotmail.com)
% last update on Feb.8, 2023


%% Inputs
demtxt = 'yourDEMtiff'; % tif file name of your DEM, include the extension, e.g. 'MDG.tif'
preciptxt = 'yourPRECIPITATIONtiff';% tif file name of your precipitation include the extension, e.g. 'MDG_precip.tif'
crita = 1e6; % minimum drainage area of channel initiation, [m^2]

mn = 0.45; % the concavity used for Chi calculation
A0 = 1e6; % the reference area used for Chi calculation, [m^2]
smoWin = 2e3; % the smooth window size used for Chi calculation, [m]
P0 = 1; % the reference precipitation for Chi calculation, [mm/yr]

%% Extract streams
DEM = GRIDobj('demtxt'); % DEM
FD = FLOWobj(DEM,'preprocess','carve'); % flow direction
A = flowacc(FD)*DEM.cellsize^2; % drainage area, [m^2]
S0 = STREAMobj(FD,flowacc(FD)>crita/DEM.cellsize^2); % streams

% calculate chi and ksn
cs = DEM.cellsize;

% option 1, not correct precipitation for Chi
PRECIP = DEM;
PRECIP.Z = PRECIP.Z./PRECIP.Z;

% % option 2, correct precipitation for Chi
% PRECIP = GRIDobj(preciptxt);

% calculate Chi in GRIDobj
[ChiGrid,ksnG] = chi_FD_AP_modified(DEM, FD, 'crita',crita,'mn',mn,'Ao',Ao, ...
    'smoWin', smoWin,'precipitation',PRECIP,'Po',P0);

%% Manually pick rivers from the mapview DEM and streams
flowpathapp(FD,DEM,S0);
cont_opt = input('When you are finished type in the name of your STREAMobj variable here: ', 's');
disp(' ');
varTest = 0;
while varTest == 0
    varTest = exist(cont_opt,'var');
    if varTest ~= 1
        txt = sprintf(['\nYour STREAMobj variable name and the name that you\n'...
            'just input do not match.\n']);
        disp(txt);
        cont_opt = input('Retype your STREAMobj variable name: ', 's');
    end
end
% close flowpathapp figures
close('Main');
close('Profiles');
S = eval(cont_opt); % S are the picked rivers

%% Manually decide river types
% the type stores the river types of all the pick rivers
type = escarpment_river_type(DEM,FD,S,ChiGrid);

%% Can proceed to analyze the statistics, or location of different river types



























