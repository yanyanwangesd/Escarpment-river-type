function outGrid = gridReSize(inGrid,xMinP,xMaxP,yMinP,yMaxP)
% gridReSize takes five inputs, an grid and maximum and minimum
% corrdinates you want to strink the grid to. The output is the resized
% grid.
% Created: 11/06/2013
% Author: Sean F. Gallen

ncols = inGrid.size(2);
nrows = inGrid.size(1);
R = inGrid.refmat;
[RSmaxR, RSminC] = UTMlatlon2pix(R,yMinP, xMinP);
[RSminR, RSmaxC] = UTMlatlon2pix(R,yMaxP, xMaxP);
GminR = 1;
GminC = 1;
GmaxC = ncols;
GmaxR = nrows;

%% This block of code resizes the grids according to smallest dimentions 
% and creates headers for the ascii file associated witht the resized grids

if RSminR >= GminR;
    minR = RSminR;
else
    minR = GminR;
end

if RSmaxR <= GmaxR;
    maxR = RSmaxR;
else
    maxR = GmaxR;
end

if RSminC >= GminC;
    minC = RSminC;
else
    minC = GminC;
end

if RSmaxC <= GmaxC;
    maxC = RSmaxC;
else
    maxC = GmaxC;
end
 outGrid = inGrid;
 
 minR = round(minR);
 maxR = round(maxR);
 minC = round(minC);
 maxC = round(maxC);
 
 rows = minR:maxR;
 cols = minC:maxC;
 outGrid.Z = inGrid.Z(rows,cols);
 outGrid.size = size(outGrid.Z);
 [yll,xll] = pix2UTMlatlon(inGrid.refmat,minR,minC);
 outGrid.refmat = makerefmat(xll, yll, inGrid.cellsize, -inGrid.cellsize);
 
end