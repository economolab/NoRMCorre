function [extended_grid,i,j,k] = getExtendedGrid(xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap,sz)
% constructs grid coordinates
% INPUTS:
% 
% overlap:      amount of overlap
% sz:           spatial size of X

% OUTPUT:
% extended_grid:            output grid coordinates

% Written by Eftychios A. Pnevmatikakis, Simons Foundation, 2016
% edited by Will Cunningham, Economo Lab, 2024

nd = length(sz);
if nd == 2; sz(3) = 1; end
[i, j, k] = ndgrid(1:length(xx_s), 1:length(yy_s), 1:length(zz_s));

extended_grid = [
    max(xx_s(i(:)) - overlap(1), 1); min(xx_f(i(:)) + overlap(1), sz(1)); ...
    max(yy_s(j(:)) - overlap(2), 1); min(yy_f(j(:)) + overlap(2), sz(2)); ...
    reshape(max(zz_s(k(:)) - overlap(3), 1),1,[]); reshape(min(zz_f(k(:)) + overlap(3), sz(3)),1,[])
];
extended_grid = reshape(extended_grid,6,length(xx_s),length(yy_s),length(zz_s));