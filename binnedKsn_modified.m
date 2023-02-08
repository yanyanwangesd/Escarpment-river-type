function [ksnSTREAMs] = binnedKsn_modified(S,chi,Sz,win,cs,Ao,mn,riversize)
% binnedKsn.m calcuates ksn from chi elevation data using a moving window
%
% Inputs:
%       1) S: a TopoToolBox STREAMobj
%       2) chi: topologically ordered vector of chi values
%       3) Sz: topologically ordered vector of elevation values
%       4) win: window size in map units
%       5) cs: cellsize in map units
%       6) Ao: reference drainage area used for chi integration
%       7) mn: reference m/n (theta) value
%
% Outputs:
%       1) topologically ordered vector of ksn values
%
% Author: Sean F. Gallen
% Date Modified: 12/31/2015

% Add note: Y.WANG: for DEM with flat areas, small drainages make the map
% looks not clean. This is a modified version of Sean Gallen's binnedKsn.m. This function
% use the number of stream nodes as a criteria to decide which drainge to
% bin Ksn. The number of river nodes is given by riversize in the input. on
% Feb. 2020

% should consider when river length is smaller than smothwindow, could use
% river length as smooth window size

% should make sure the input Sz doesn't contain NaN values, the input chi
% and Sz should be strictly one to one, each non-nan Sz has its own non-nan
% chi value. Otherwise, the function will fail in the loop, gives warning
% of array size not match (Y.WANG, Feb. 2021)

%!!!!! add piece to examine Sz, Y.WANG in Feb. 2021. Problem arises if use
%HydroSHED flow direction for stream extraction, streams might flow on DEM
%pixels with NaN value, expecially at coastal areas. So need to assign
%values for these elevation nan pixels at stream point. Found this problem
%for Madagascar, and it is the coastal empty pixels, so we assign elevation
%of these pixels to be 0.
Sz(isnan(Sz)) = 0;
%!!!!! add piece to examine Sz, Y.WANG in Feb. 2021.



step = (win/2);
indStep = ceil(step/cs);
% declare necessary variables from stream object
ordList = S.orderednanlist;         % ordered list of streams split by nans
strmBreaks = find(isnan(ordList));  % get position of nans
Six = S.ix;                         % doners
Sixc = S.ixc;                       % recievers
Sd = S.distance;                    % distance from mouth

% cast dumby vector to catch data
ksnSTREAMs = nan(size(Sd));

% from the donor position identifiy receiver nodes in stream object order
receiverID = nan(size(Sd));
donorID = nan(size(Sd));
for i = numel(Six):-1:1
    receiverID(Six(i)) = Sixc(i);
    donorID(Six(i)) = Sixc(i);
end


id1 = 0;
h = waitbar(0,'Calculating ksn from chi for ksn map...');
for i = 1:length(strmBreaks)
    strmInds = ordList(id1+1:strmBreaks(i)-1);
    if  length(strmInds) >= riversize  % !!!!! only smooth for rivers longer than riverzie, added by Y.WANG
        tChi = chi(strmInds);
        tZ = Sz(strmInds);
        tD = Sd(strmInds);
        trecID = receiverID(strmInds);
        tKsn = nan(length(strmInds),1);
        % if the stream is a trunk channel run smoothing window
        if isnan(trecID(end))
            for j = 1:length(tChi)
                % find within the window size.
                chiWin = tChi(tD > tD(j)-win & tD <= tD(j)+win);
                zWin = tZ(tD > tD(j)-win & tD <= tD(j)+win);

                % make sure there are enough points for the regression
                if length(chiWin) > 2
                    chi_segMat = [ones(size(zWin)) chiWin];
                    [b,~,~,~,~] = regress(zWin,chi_segMat,0.05);
                    tKsn(j) = b(2)*Ao^mn;
                else
                    tKsn(j) = nan;
                end
            end
        % if the stream is a tributary allow smoothing window to continue down
        % stream by 1/2 of moving window width
        else
            addChi = nan(1,indStep);
            addZ = nan(1,indStep);
            addD = nan(1,indStep);
            nInd = trecID(end);
            n = 1;
            while n < length(addChi)
                addChi(n) = chi(nInd);
                addZ(n) = Sz(nInd);
                addD(n) = Sd(nInd);
                nInd = receiverID(nInd);
                if isnan(nInd) == 1 || isempty(addChi(isnan(addChi))) == 1
                    n = length(addChi) + 1;
                else
                    n = n+1;
                end
            end
                        
            addChi = addChi(~isnan(addChi));
            addZ = addZ(~isnan(addZ));
            addD = addD(~isnan(addD));
            
            tChi2 = nan(length(tChi)+length(addChi),1);
            tZ2 = nan(length(tZ)+length(addZ),1);
            tD2 = nan(length(tD)+length(addD),1);

            tChi2(1:length(tChi)) = tChi;
            tChi2(length(tChi)+1:end) = addChi;
            tZ2(1:length(tZ)) = tZ;
            tZ2(length(tZ)+1:end) = addZ;
            tD2(1:length(tD)) = tD;
            tD2(length(tD)+1:end) = addD;

            for j = 1:length(tChi)
                % find within the window size.
                chiWin = tChi2(tD2 > tD2(j)-win & tD2 <= tD2(j)+win);
                zWin = tZ2(tD2 > tD2(j)-win & tD2 <= tD2(j)+win);

                % make sure there are enough points for the regression
                if length(chiWin) > 2
                    chi_segMat = [ones(size(zWin)) chiWin];
                    [b,~,~,~,~] = regress(zWin,chi_segMat,0.05);
                    tKsn(j) = b(2)*Ao^mn;
                else
                    tKsn(j) = nan;
                end
            end
        end
        tKsn(tKsn < 0) = 0;
        ksnSTREAMs(strmInds) = tKsn;
    else
         ksnSTREAMs(strmInds) = nan;
        
        
    end
% %%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%  buggy = num2str(i);
%  buggy = strcat('i = ',buggy);
%  disp(buggy)
% 
% %%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    f = i/length(strmBreaks);
    waitbar(f,h);
    id1 = strmBreaks(i);
end


close(h)


end


