function res = bpass(image, lnoise, lobject, varargin)
% FUNCTION NAME:
%   bpass
%
% DESCRIPTION:
%   Implements a real-space bandpass filter which surpressed pixel noise 
%   and long-wavelength image variations while retaining information about 
%   features of a characteristic size.
%
% INPUT (REQUIRED)
%           image: (int8) image to be filtered 
%          lnoise: High spatial frequnecy cutoff. Often 1 pixel for
%                  pixel-to-pixel noise
%         lobject: Low spatial frequency cutoff for slow varations in the
%                  image. Often a little larger than feature size. 
%
% INPUT (OPTIONAL)
%          noclip: Keep all values in output image, even negative values.  
%
% OUTPUT:
%             res: (single) filtered image
%
% CALLING SEQUENCE:
%   b = bpass(img,1,11)
%
% NOTES :
%   IDL Version - This code was translated from bpass.pro code provided on
%           Eric Weeks' website:
%           https://physics.emory.edu/faculty/weeks/idl/yykit/msd.pro
%
% REVISION HISTORY:
%   02/??/1993 - David G. Grier
%       * Wrote orginal version at UChicago
%   05/??/1995 - David G. Grier 
%       * Greatly revised version 
%   05/??/1997 - John C. Crocker 
%		* Added /field keyword 
%   10/23/2022 - K Aptowicz (WCU)
%       * Translated to MATLAB
%

%% Reading and setting parameters
% Set default values for optional parameters
default_noclip = [];

% Create fields for all optionals inputs
p = inputParser;
% Keywords
addOptional(p,'noclip', default_noclip)

% populate optional parameters from inputs
parse(p,varargin{:})
noclip = p.Results.noclip;

%% Setup kernels for convolution
nf = size(image,3);
b = single(lnoise);
w = round(max(lobject,(2.*b)))';
N = 2*w + 1;

r = (single([0:N-1])-w)/(2.*b);
xpt = exp(-r.^2);
xpt = xpt/sum(xpt);
factor = (sum(xpt.^2) - 1/N);

% Kernel to remove high spatial frequency
gx = xpt;
gy = gx';

% Kernel to remove low spatial frequency
bx = single(zeros(N,1)) - 1./N;
by = bx';


%% Calculate convolution on each image in stack
res = single(image);
for i = 1:nf 
    % Edges are padded with zeros for convn. To minimize artifacts, set 
    % mean pixel value of image to zero
    img=single(image(:,:,i));
    img=img-mean(img(:));

    % Remove high spatial frequency noise
    g = convn(img,gx,'same');
	g = convn(g,gy,'same');
    
    % Calculate low spatial frequency background 
	b = convn(img,bx,'same');
	b = convn(b,by,'same');
	
    % Remove low spatial frequency background
	res(:,:,i) = g-b;
end

res=255*res/max(res(:));  % Rescale max to 8-bit 

% By default, negative values are set to zero, but this can be skipped
% using the noclip parameter (noclip = 'y')
if isempty(noclip)
	res(res<0) = 0;
end














