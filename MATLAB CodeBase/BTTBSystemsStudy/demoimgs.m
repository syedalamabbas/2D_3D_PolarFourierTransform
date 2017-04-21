function varargout=demoimgs
    pth = fileparts(which('cameraman.tif'));
    D = dir(pth);
    C = {'.tif';'.jp';'.png';'.bmp'};
    idx = false(size(D));
    for ii = 1:length(C)
        idx = idx | (arrayfun(@(x) any(strfind(x.name,C{ii})),D));
    end
    D = D(idx);
    for ii = 1:numel(D)
        fprintf('%s\n',D(ii).name)
    end
    if nargout, varargout{1}=pth; 
    end
end