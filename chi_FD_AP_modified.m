function [ChiGrid, ksnG] = chi_FD_AP_modified(DEM, FD,  varargin)
% chi_profiler.m allows for interative river profile analysis via
% the integral or chi method (e.g. Perron and Royden, 2013). Users can
% regress through channel segments to get ksn and select knickpoints along
% river profiles
%
% Inputs:
%       1) DEM: TopoToolBox DEM GRIDobj (required)
%       2) S: TopoToolBox STREAMobj (required)
%       4) crita: Critical drainage area for channel head initiation in map
%          units. (optional) {default --> 1e6}
%       5) mn: reference m/n (theta) value (optional) {default --> 0.45}
%       6) Ao: reference drainage area for chi integration (optional)
%          {default --> 1}. Note: it is recommended to always use an 'Ao'
%          of 1 as this results in chi versus elevation plots with a slope
%          that is equal to the normalized steepness index.
%       7) smoWin: size of window used to smooth elevation data (set this
%          to the cell size of the DEM if you don't want the data to be
%          smoothed) (optional) {default --> 250}
%       8) flowOption: string of 'fill' or 'carve' for flow routing.
%          (optional) {default --> empty assumes that DEM is already 
%          filled or carved}
%       
% Outputs:
%       No variable outputs to Matlab, but chi_profiler.m will produce a
%       series of files that can be imported into Matlab and ArcGIS.
%
%       - Data tables of each of the stream channels analyzed
%       - Shapefiles of chi, ksn, regressed ksn segments and knickpoints
%       - Tables of regressed ksn and knickpoint statistics
%       Optional:
%       - river profiler figures
%       - chi and ksn map of entire river network in DEM as a shapefile
%
%       See 'User_Guide.docx' for more information on the how to run
%       chi_profiler.m and the output file organization
%
% Example:
%  chi_profiler(DEM,S,'my_proj','crita',1e6,'mn',0.45,'Ao',1,'smoWin',250);
%
%
% Author: Sean F. Gallen
% Date Modified: 02/20/2017
% email: sean.gallen[at]erdw.ethz.ch

% this script is modified to include option of precipitation for CHI calculation and
% use flowdirection (e.g. hydroSHED corrected ) other than derived from
% DEM. Add on March. 2020, by Yanyan Wang, wangyanyan0607@hotmail.com

% Parse Inputs
p = inputParser;         
p.FunctionName = 'chi_FD_AP';

% required inputs
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
addRequired(p,'FD',@(x) isa(x,'FLOWobj'));

% optional inputs
addOptional(p,'crita', 1e6, @(x) isscalar(x));
addOptional(p,'mn', 0.45, @(x) isscalar(x));
addOptional(p,'Ao', 1, @(x) isscalar(x));
addOptional(p,'smoWin', 250, @(x) isscalar(x));
addOptional(p,'precipitation',  @(x) isa(x,'GRIDobj'));
addOptional(p,'Po', 1000, @(x) isscalar(x));

parse(p,DEM, FD,  varargin{:});
DEM   = p.Results.DEM;
FD = p.Results.FD;

crita    = p.Results.crita;
mn = p.Results.mn;
Ao = p.Results.Ao;
smoWin = p.Results.smoWin;
P = p.Results.precipitation;
Po =  p.Results.Po;


%% resize grids to make DEM and PPRECIP align
% make sure grids are the same size
[xMin(1,1),xMax(1,1),yMin(1,1),yMax(1,1)] = findCorners(DEM);
[xMin(2,1),xMax(2,1),yMin(2,1),yMax(2,1)] = findCorners(P);

% resizing grids to largest overlapping area
xMinP = max(xMin);
xMaxP = min(xMax);
yMinP = max(yMin);
yMaxP = min(yMax);

DEM = gridReSize(DEM,xMinP,xMaxP,yMinP,yMaxP);
P = gridReSize(P,xMinP,xMaxP,yMinP,yMaxP);

P.refmat = DEM.refmat;

% @YW, fix bug: precipitation is empty at some non-nan DEM grid 
nanDEMgrid = isnan(DEM.Z(:));
nanPgrid = isnan(P.Z(:));

pemtpy = nanPgrid.*~nanDEMgrid ;

IND = find(pemtpy);

for i = 1:length(IND)
    
    [Ip, Jp] = ind2sub(P.size,IND(i));
    
    % find neighbour nodes    
    nbI= [Ip, Ip+1, Ip, Ip+1, Ip-1, Ip-1,Ip-1, Ip+1]; 
    nbJ = [Jp-1, Jp+1, Jp+1,  Jp, Jp, Jp-1, Jp+1, Jp-1];
    
    boundary = [0 0  P.size(2)+1 P.size(1)+1];
    
    [lia1, ~] = ismember(nbI, boundary);
    [lia2, ~] = ismember(nbJ, boundary);
    
    lia = lia1|lia2;        
    neighbour = sub2ind(P.size, nbI(~lia),nbJ(~lia));    
    % empty node value is extrapolated from average values of its direct 8 neighbours
    P.Z(IND(i)) = mean( P.Z(neighbour),'omitnan');
    
end

% FD = FLOWobj(DEM,'preprocess','carve');

%% create varables with topotoolbox functions
% declare cellsize
cs = DEM.cellsize;

% Calculate flow accumulation (A)
A   = flowacc(FD).*(cs^2);

S1 = STREAMobj(FD,'minarea',crita/(DEM.cellsize^2));
S1.refmat = DEM.refmat;

% Declare STREAMobj variables for faster processing through forloop
ordList = S1.orderednanlist;
strmBreaks = find(isnan(ordList));

GridID = S1.IXgrid;

Sz = double(DEM.Z(GridID));                 % elevation
SmoZ = Sz;               % dumby vector to get smoothed data



% get variables ready for chi integration
chis = zeros(size(S1.distance));
Six = S1.ix;
Sixc = S1.ixc;
Sx = S1.distance;
AP = A.Z(S1.IXgrid).*P.Z(S1.IXgrid);
Sa = (Ao*Po./AP).^mn;

h = waitbar(0,'calculating \chi for full stream network...');
% calculating chi for the entire river network
for lp = numel(Six):-1:1
    chis(Six(lp)) = chis(Sixc(lp)) + (Sa(Sixc(lp))+(Sa(Six(lp))-Sa(Sixc(lp)))/2) *(abs(Sx(Sixc(lp))-Sx(Six(lp))));
    f = (numel(Six)+1 - lp)/numel(Six);
    waitbar(f,h);
end
close(h);

ChiGrid = DEM;
ChiGrid.Z = nan(size(DEM.Z));
ChiGrid.Z(S1.IXgrid) = chis;

% plot all of the river river profile data as thin gray lines
h = waitbar(0,'Smoothing elevation data for full stream network...');
id1 = 0;
for i = 1:length(strmBreaks)
    strmInds = ordList(id1+1:strmBreaks(i)-1);
    SmoZ(strmInds) = smoothChannelZ(Sz(strmInds),smoWin,cs);
    id1 = strmBreaks(i);
    f = i/length(strmBreaks);
    waitbar(f,h);
end
close(h)

%!!!!! add piece to examine Sz, Y.WANG in Feb. 2021. Problem arises if use
%HydroSHED flow direction for stream extraction, streams might flow on DEM
%pixels with NaN value, expecially at coastal areas. So need to assign
%values for these elevation nan pixels at stream point. Found this problem
%for Madagascar, and it is the coastal empty pixels, so we assign elevation
%of these pixels to be 0.
SmoZ(isnan(SmoZ)) = 0;
%!!!!! add piece to examine Sz, Y.WANG in Feb. 2021.


%%Comment out if don't want to calcualte ksn from chi, the binnedKsn deal with imcomplemte situations, could end in error
riversize = ceil(smoWin/DEM.cellsize); % the river size is in unit of grid size
ksnSTREAMs = binnedKsn_modified(S1,chis,SmoZ,smoWin,cs,Ao,mn,riversize);
ksnG = DEM;
ksnG.Z = nan(size(DEM.Z));
ksnG.Z(GridID) = ksnSTREAMs;
% 
% %% finally ask the user if they would like to make a chi map and ksn map
% % written as a shapefile for the entire drainage network
% 
% S = STREAMobj(FD,'minarea',50e6/(DEM.cellsize^2));
% S.refmat = DEM.refmat;
% 
% MS = STREAMobj2mapstruct(S1,'seglength',seglenth,'attributes',...
%     {'chi' ChiGrid @mean...
%     'ksn' ksnG @mean});
% fileName = [chanDir,'/','MDGp1418_1e6', '_chi_ksn_map.shp'];
% shapewrite(MS,fileName);
% 

end


