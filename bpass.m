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
%   03/12/2024 - K Aptowicz (WCU)
%       * Tweaked to match IDL version. Padded zeros for result of 
%         convolution. Fixed an issue with factor. Removed subtracting the 
%         mean of the image (not sure why that was there). Now the result 
%         almost perfectly matches the Eric Weeks IDL bpass result.  There 
%         appears to be an issue with 'factor' being slightly different. I 
%         think it has do do with the single precision format.  Good
%         enough.
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
factor = (sum(xpt.^2) - round(1/N));
% KBA: Added round(1/N) to match IDL code.  Not sure why this is needed. 

% Kernel to remove high spatial frequency
gx = xpt;
gy = gx';

% Kernel to remove low spatial frequency
bx = single(zeros(N,1)) - 1./N;
bx = bx'; % orient it correctly
by = bx';


%% Calculate convolution on each image in stack
res = single(image);
for i = 1:nf 
    img=single(image(:,:,i));

    % Remove high spatial frequency noise
    g = convn(img,gx,'same');
    g(:,1:w)=0; g(:,end-w+1:end)=0; %KBA - Added to mirror IDL code

	g = convn(g,gy,'same');
    g(1:w,:)=0; g(end-w+1:end,:)=0; %KBA - Added to mirror IDL code

    % Calculate low spatial frequency background 
	b = convn(img,bx,'same');
    b(:,1:w)=0; b(:,end-w+1:end)=0; %KBA - Added to mirror IDL code

	b = convn(b,by,'same');
    b(1:w,:)=0; b(end-w+1:end,:)=0; %KBA - Added to mirror IDL code
	
    % Remove low spatial frequency background
	res(:,:,i) = g-b;
end

res=res/factor;

% By default, negative values are set to zero, but this can be skipped
% using the noclip parameter (noclip = 'y')
if isempty(noclip)
	res(res<0) = 0;
end















