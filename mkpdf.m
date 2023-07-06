function dx = mkpdf(tr,dt,varargin)
% FUNCTION NAME:
%   mkpdf
%
% DESCRIPTION:
%   Calculates the change in particle positions (dx, dy) for a delay time, dt. 
%
% INPUT (REQUIRED)
%              tr: Track array.
%		       dt: delay time between frames to compare
%
% INPUT (OPTIONAL)
%		       dim: [default: 2] Spatial dimensions
%
% OUTPUT:
%		    output: [dx, dy]
%
% CALLING SEQUENCE:
%   dxdy = mkpdf(t,1)
%
% NOTES :
%   IDL Version - This code was translated from mkpdf.pro code provided on
%           Eric Weeks' website:
%           https://physics.emory.edu/faculty/weeks/idl/kit/mkpdf.pro
%
% REVISION HISTORY:
%   05/20/2005 - Eric Weeks
%       * Write orginal version
%   07/05/2023 - Kevin Aptowicz
%       * Translated to MATLAB
%

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