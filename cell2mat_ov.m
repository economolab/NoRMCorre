function X = cell2mat_ov(I,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap,sz)

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
nonBorderCols = any(sameSize~=0,1);
nonBorderRows = any(sameSize~=0,2);
temp = cat(3,I{nonBorderRows,nonBorderCols});
temp = reshape(temp,patchSize(1),patchSize(2),sum(nonBorderRows),sum(nonBorderCols));
temp = permute(temp,[1,3,2,4]);
temp = reshape(temp,patchSize(1)*sum(nonBorderRows),patchSize(2)*sum(nonBorderCols));

idxBorderRows = find(~nonBorderRows);
preBorderRows = idxBorderRows(idxBorderRows<find(nonBorderRows,1));
postBorderRows = idxBorderRows(idxBorderRows>find(nonBorderRows,1));
for i = 1:sum(preBorderRows>0)
    temp = cat(1,cat(2,I{preBorderRows(i),nonBorderCols}),temp);
end
for i = 1:numel(postBorderRows)
    temp = cat(1,temp,cat(2,I{postBorderRows(i),nonBorderCols}));
end
idxBadCols = find(~nonBorderCols);
preBadCols = idxBadCols(idxBadCols<find(nonBorderCols,1));
postBadCols = idxBadCols(idxBadCols>find(nonBorderCols==0,1));
for i = 1:numel(preBadCols)
    temp = cat(2,cat(1,I{:,preBadCols(i)}),temp);
end
for i = 1:numel(postBadCols)
    temp = cat(2,temp,cat(1,I{:,postBadCols(i)}));
end
% figure;imagesc(temp);axis image;colormap gray;

% Initialize the row transformation matrix
A = zeros(sz(1),size(temp,1),dType);

toAdd = (eye(patchSize(1),dType).*[1:overlap(1) (overlap(1)*ones(1,patchSize(1)-2*overlap(1),dType)) overlap(1):-1:1]);
for i = 1:ceil(size(temp,1)/patchSize(1))
    
    idx2 = -overlap(1) + patchSize(1)*(i-1)+1; idx2 = idx2:idx2+patchSize(1)-1;
    idx = idx2-2*overlap(1)*(i-1);
    goodDex1 = idx>0 & idx<=sz(1); idx = idx(goodDex1);
    goodDex2 = idx2>0 & idx2<=size(temp,1); idx2 = idx2(goodDex2);
    A(idx,idx2) = A(idx,idx2) + toAdd(goodDex1,goodDex2);
end
A = A./sum(A,2);

% Initialize the row transformation matrix
B = zeros(sz(2),size(temp,2),dType);
for i = 1:ceil(size(temp,2)/patchSize(2))
    idx2 = -overlap(2) + patchSize(2)*(i-1)+1; idx2 = idx2:idx2+patchSize(2)-1;
    idx = idx2-2*overlap(2)*(i-1);
    goodDex1 = idx>0 & idx<=sz(2); idx = idx(goodDex1);
    goodDex2 = idx2>0 & idx2<=size(temp,2); idx2 = idx2(goodDex2);
    B(idx,idx2) = B(idx,idx2)  + toAdd(goodDex1,goodDex2);
end
B = B./sum(B,2);


X = (A * temp * B');
% figure; imagesc(vm(cat(3,XC,X)))
% figure;moviesc(vm(gather(cat(3,X,(A * temp * B')./denom))));
% figure;moviesc(vm(gather(cat(3,B,A * ones(size(temp)) * B2'))));axis image;colormap gray;

end
%