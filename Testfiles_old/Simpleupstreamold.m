function mup = Simpleupstreamold(G,flag, m)
N  = double(G.faces.neighbors);
intInx = all(N ~= 0, 2);
N  = N(intInx, :);
if numel(flag) == 1
    flag = repmat(flag, size(N, 1), 1);
end
assert(numel(flag) == size(N, 1) && islogical(flag), ...
    'One logical upstream flag must'' be supplied per interface.');
upCell       = N(:,2);
upCell(flag) = N(flag,1);
if isnumeric(m)
    % m is a simple matrix, we just extract values using the cells
    mup = m(upCell, :);
else
    display('mobilities are likely ADI framework')
end