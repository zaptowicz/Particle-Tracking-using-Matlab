function dx = mkpdf(tr,dt,varargin)

% mkpdf  --  Eric R. Weeks  --  new version 5-20-05
%
% see http://www.physics.emory.edu/~weeks/idl/
%    for more information and software updates
%
%
%
% basically this utilizes very similar code to getdx.pro, which
% I wrote a while ago.

%% Reading and setting parameters
% Set default values for optional parameters
default_dim = 2;

% Create fields for all optionals inputs
p = inputParser;

% Variables
addParameter(p,'dim',default_dim,@isnumeric)

% Keywords
%addOptional(p,'flag', default_flag)

% populate optional parameters from inputs
parse(p,varargin{:})
dim_input = p.Results.dim;

if isempty(dt) 
    dt=1;
end

dx=getdx(tr,dt,dim=dim_input);
dim=dim_input;
w=find(dx(:,dim+1) > -0.5); nw=numel(w);
if (nw > 0)
	dx=dx(w,1:dim);
else
	dx = 0.0
end

end