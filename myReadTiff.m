function out = myReadTiff(fullfn)
    out = tiffreadVolume(fullfn);
    info = imfinfo(fullfn);   
    Nrow = info(1).Height;
    Ncol = info(1).Width;
    [nchans,nplanes] = parseTiffChansPlanesInfo(info(1));
    out = reshape(out, Nrow, Ncol, nchans, nplanes);
end