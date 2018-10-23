function varargout = applyoffset(xoff,varargin)
%% Apply offset in x-direction
% First input must be offset as single number. Afterwards, any number of
% point struct cells can be input and the offset will be applied to each of
% them in turn

% Output should be same size cell as input
varargout = cell(size(varargin));

for ii=1:length(varargin)
    
    new = varargin{ii};
    
    for i=1:length(new)
        new(i).x = new(i).x+xoff;
    end
    
    varargout{ii} = xyztopoints(new);
end