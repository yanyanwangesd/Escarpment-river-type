function type = escarpment_river_type(DEM,FD,S,Chi)
% Function name: escarpment_river_type.m
% Author: Yanyan Wang
% Date modified: 08/02/2023
% Purpose: quanlitatively determine escarpment river types with normal
% river profile, Chi-elevation profile, and Chi-steepness profile. The
% potential application includes capture analysis, divide asymmetries

%
% Inputs: 
%       1) DEM, GRIDobj
%       2) FD, flow direction, FLOWobj
%       3) S, streams to be analyzed, STREAMobj
%       4) Chi, GRIDobj

%
% Output:
%       1) type, array of river types determined by user
%
% Related publication: Wang et al 2023, The role of weathering on morphology and rates of escarpment retreat of the rift margin of Madagascar
%                      on JGR-Earth Surface
% Please cite this paper if you use this function.


p = inputParser;         
p.FunctionName = 'escarpment_river_type';

% required inputs
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
addRequired(p,'FD', @(x) isa(x,'FLOWobj'));
addRequired(p,'S', @(x) isa(x,'SREAMobj'));
addRequired(p,'Chi', @(x) isa(x,'GRIDobj'));

parse(p,DEM,FD,S,Chi);
S   = p.Results.S;
DEMesp = p.Results.DEM;
CHIesp = p.Results.Chi;
FDesp =  p.Results.FD;

%% 1) stream nodal metrics
IXgrid = S.IXgrid;
ordList = S.orderednanlist; 
strmBreaks = find(isnan(ordList));

%% 2) find channel head GRID index
headsGRID = nan(length(strmBreaks),1);
id=0;
for i = 1:length(strmBreaks)    
    headsGRID(i) = ordList(id+1);
    id = strmBreaks(i);    
end
headsGRID = IXgrid(headsGRID);

%% 3) sort rivers according to the channel head latitude, is optional here
% [x,y] =  getcoordinates(DEMesp);
% Y2 = repmat(y,1,length(x));
% [~, neworder]= sort(Y2(headsGRID));
% headsGRIDsort = headsGRID(neworder);
headsGRIDsort = headsGRID; % do not re-order the river sequences

%% 4) plot profiles and determine river type
type = nan(length(strmBreaks),1);
for i = 1:length(strmBreaks)

    Splot = STREAMobj(FDesp,'channelheads',headsGRIDsort(i));           
    dist2 = Splot.distance;
    IXgrid2 = Splot.IXgrid;
    chip2 = CHIesp.Z(IXgrid2);
    zp2 =double(DEMesp.Z(IXgrid2));    
    zp2(isnan(zp2)) = 0; % this is for basins covering coastal plain, 
    ordList2 = Splot.orderednanlist; 
    strmBreaks2 = find(isnan(ordList2));
    id2=0;
    ii=1;
    strmInds2 = ordList2(id2+1:strmBreaks2(ii)-1);   
    zp2(strmInds2) = smoothChannelZ(zp2(strmInds2),2000,DEMesp.cellsize);
    displot = dist2(strmInds2);
    zplot = zp2(strmInds2);
    chiplot = chip2(strmInds2);
    idnan = ~isnan(chiplot);
    chiplot2 = interp1(displot(idnan),chiplot(idnan),displot);
    chiplot2 = smoothChannelZ(chiplot2,2000,DEMesp.cellsize);
    ksn = [diff(zplot); 0]./[diff(chiplot2);0]; 

   
    fh = figure(i+100);
    subplot(3,1,1)
    plot(displot/1000,zplot,'-','lineWidth',4,'color','k'); 
    xlabel('Distance (km)','FontSize',24); ylabel('Elevation (m)'); 
    set(gca,'FontSize',16)
    axis tight
    box on
    grid on
    subplot(3,1,2)
    plot( chiplot2,zplot,'k-','lineWidth',4,'color','r'); hold on
    xlabel('\chi ','FontSize',24); ylabel('Elevation (m)');
    set(gca,'FontSize',16)
    box on
    grid on
    axis tight
    subplot(3,1,3)
    hold on
    plot(chiplot2,ksn,'m-','lineWidth',4,'color','m');   
    xlabel('\chi'); ylabel('K_{sn}');  
    set(gca,'FontSize',16)
    axis tight
    box on
    grid on
    
    pos = get(fh,'position');
    set(fh,'position',pos.*[1 1 1.5 1.5]) % make window size 1.5*1.5 bigger

    % User type in the river type in the Command window
    fprintf('\n ********* User Input Needed *********\n')
    fprintf('\n CLASSIFY the river type with a number, for example, 1=divide-type, 2=knickzone-type, etc... \n')
    prompt = 'Your river type is: ';
    type(i) = input(prompt);

   
    close all
    
          
end



