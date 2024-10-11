function X = cell2mat_ov_sum(I,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap,sz,Bs)

% converts a cell array to a matrix when the cell elements overlap
% INPUTS:
% I:            cell array
% grid_size:    true size of each element
% overlap:      amount of overlap in each direction
% d1:           number of rows of matrix
% d2:           number of columns of matrix

% OUTPUT:
% X:            output matrix

% Old version Written by Eftychios A. Pnevmatikakis, Simons Foundation, 2016
% New version (much faster) written by Will Cunningham Economo Lab 2024

if isa(I{1},'gpuArray')
    dType = 'gpuArray';
else
    dType = 'single';
end
patchSize = size(I{round(size(I,1)/2),round(size(I,2)/2)});
sameSize = cellfun(@(X) all(size(X)==patchSize),I);
goodCols = any(sameSize~=0,1);
goodRows = any(sameSize~=0,2);
temp = cat(3,I{goodRows,goodCols});
temp = reshape(temp,patchSize(1),patchSize(2),sum(goodRows),sum(goodCols));
temp = permute(temp,[1,3,2,4]);
temp = reshape(temp,patchSize(1)*sum(goodRows),patchSize(2)*sum(goodCols));

idxBadRows = find(~goodRows);
preBadRows = idxBadRows(idxBadRows<find(goodRows,1));
postBadRows = idxBadRows(idxBadRows>find(goodRows,1));
for i = 1:sum(preBadRows>0)
    temp = cat(1,cat(2,I{preBadRows(i),goodCols}),temp);
end
for i = 1:numel(postBadRows)
    temp = cat(1,temp,cat(2,I{postBadRows(i),goodCols}));
end
idxBadCols = find(~goodCols);
preBadCols = idxBadCols(idxBadCols<find(goodCols,1));
postBadCols = idxBadCols(idxBadCols>find(goodCols==0,1));
for i = 1:numel(preBadCols)
    temp = cat(2,cat(1,I{:,preBadCols(i)}),temp);
end
for i = 1:numel(postBadCols)
    temp = cat(2,temp,cat(1,I{:,postBadCols(i)}));
end
% figure;imagesc(temp);axis image;colormap gray;





% Initialize the row transformation matrix
A = zeros(sz(1),size(temp,1),dType);

for i = 1:ceil(size(temp,1)/patchSize(1))
    idx2 = patchSize(1)*(i-1)+1; idx2 = idx2:idx2+patchSize(1)-1;
    idx = idx2-2*overlap(1)*(i-1);
    toAdd = (eye(patchSize(1),dType).*[1:overlap(1) (overlap(1)*ones(1,patchSize(1)-2*overlap(1),dType)) overlap(1):-1:1]);
    goodDex1 = idx>0 & idx<=sz(1); idx = idx(goodDex1);
    goodDex2 = idx2>0 & idx2<=size(temp,1); idx2 = idx2(goodDex2);
    A(idx,idx2) = A(idx,idx2) + toAdd(goodDex1,goodDex2);
end
A = A./sum(A,2);

% Initialize the row transformation matrix
B = zeros(sz(2),size(temp,2),dType);
for i = 1:ceil(size(temp,2)/patchSize(2))
    idx2 = patchSize(2)*(i-1)+1; idx2 = idx2:idx2+patchSize(2)-1;
    idx = idx2-2*overlap(2)*(i-1);
    toAdd = (eye(patchSize(2),dType).*[1:overlap(2) (overlap(2)*ones(1,patchSize(2)-2*overlap(2),dType)) overlap(2):-1:1]);
    goodDex1 = idx>0 & idx<=sz(2); idx = idx(goodDex1);
    goodDex2 = idx2>0 & idx2<=size(temp,2); idx2 = idx2(goodDex2);
    B(idx,idx2) = B(idx,idx2)  + toAdd(goodDex1,goodDex2);
end
B = B./sum(B,2);


X = (A * temp * B');
% figure;moviesc(vm(gather(cat(3,X,(A * temp * B')./denom))));
% figure;moviesc(vm(gather(cat(3,B,A * ones(size(temp)) * B2'))));axis image;colormap gray;

end

function X = oldVersion(I,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap,sz,Bs)
tic
if isa(I{1},'gpuArray')
    dType = 'gpuArray';
else
    dType = 'single';
end
if nargin < 10 || isempty(Bs)
    Bs = cellfun(@(x) ones(size(x),dType), I,'un',0);
end
X = zeros([sz,size(I{1,1},length(sz)+1)],dType);
B = zeros(size(X),dType);
if length(sz) == 2; sz(3) = 1; end

for i = 1:length(xx_f)
    for j = 1:length(yy_f)
        for k = 1:length(zz_f)
            extended_grid = [max(xx_s(i)-overlap(1),1),min(xx_f(i)+overlap(1),sz(1)),max(yy_s(j)-overlap(2),1),min(yy_f(j)+overlap(2),sz(2)),max(zz_s(k)-overlap(3),1),min(zz_f(k)+overlap(3),sz(3))];
            %W = construct_weights([xx_s(i),xx_f(i),yy_s(j),yy_f(j),zz_s(k),zz_f(k)],extended_grid)';   
            
            Xtemp = zeros(size(I{i,j,k}),dType);
            Btemp = Xtemp;
            ind = ~isnan(I{i,j,k});
            Xtemp(ind) = Bs{i,j,k}(ind).*I{i,j,k}(ind);
            Btemp(ind) = Bs{i,j,k}(ind);
            %X(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) = X(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) + Bs{i,j,k}.*I{i,j,k}; 
            %B(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) = B(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) + Bs{i,j,k}.*(I{i,j,k}~=0);
            X(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) = X(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) + Xtemp; 
            B(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) = B(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) + Btemp;            
        end
    end
end

X = X./B;
toc
end
%