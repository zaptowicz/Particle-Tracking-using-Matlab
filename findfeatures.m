function arr = findfeatures(image, extent, varargin)
% FUNCTION NAME:
%   findfeatures
%
% DESCRIPTION:
%   Finds and measures roughly circular 'features' within
%	an image.
%
% INPUT (REQUIRED)
%           image: (int8) image with 'blurbs' to be analyzed
%          extent: a parameter which should be a little greater than
%			       the diameter of the largest features in the image.
%                  Extent MUST BE ODD valued.
%
% INPUT (OPTIONAL)
%           noclip: ['y'] Keep all values in output image, even negative values.
%		separation: An optional parameter which specifies the
%			        minimum allowable separation between feature
%			        centers. The default value is diameter+1.
%		   masscut: Setting this parameter saves runtime by reducing the
%			        runtime wasted on low mass 'noise' features.
%		    minpix: Set this optional parameter to the minimum allowed
%			        value for the peak brightness of a feature. Useful
%			        for limiting the number of spurious features in
%			        noisy images.
%		     quiet:	['y'] Supress printing of informational messages.
%		   iterate: ['y'] If the refined centroid position is too far from
%			        the initial estimate, iteratively recalc. the centroid
%			        using the last cetroid to position the mask.  This
%			        can be useful for really noisy data, or data with
%			        flat (e.g. saturated) peaks.  Use with caution- it
%			        may 'climb' hills and give you multiple hits.
%
% OUTPUT:
%		    f(:,1):	this contains the x centroid positions, in pixels.
%		    f(:,2): this contains the y centroid positions, in pixels.
%		    f(:,3): this contains the integrated brightness of the
%			        features.
%		    f(:,4): this contains the square of the radius of gyration
%			        of the features.
%		    f(:,5): this contains the eccentricity, which should be
%			        zero for circularly symmetric features and order
%			        one for very elongated images.
%
% CALLING SEQUENCE:
%   f = feature(b,1,11)%
%   f = feature(b,1,11,minpix=150,quiet='y')%
%
% NOTES :
%   IDL VERSION
%           This code was translated from feature.pro code provided
%           on Eric Weeks' website:
%           https://physics.emory.edu/faculty/weeks/idl/kit/feature.pro
%   RESTRICTIONS
%		To work properly, the image must consist of bright,
%		circularly symmetric regions on a roughly zero-valued
%		background. To find dark features, the image should be
%		inverted and the background subtracted. If the image
%		contains a large amount of high spatial frequency noise,
%		performance will be improved by first filtering the image.
%		BPASS will remove high spatial frequency noise, and
%		subtract the image background and thus provides a good
%		complement to using this program. Individual features
%		should NOT overlap.
% PROCEDURE:
%		First, identify the positions of all the local maxima in
%		the image (defined in a circular neighborhood with diameter
%		equal to 'diameter'). Around each of these maxima, place a
%		circular mask, of diameter 'diameter', and calculate the x & y
%		centroids, the total of all the pixel values, and the radius
%		of gyration and the 'eccentricity' of the pixel values within
%		that mask. If the initial local maximum is found to be more
%		than 0.5 pixels from the centroid and iterate is set, the mask
%		is moved and the data are re-calculated. This is useful for
%		noisy data. If the restrictions above are adhered to, and the
%		features are more than about 5 pixels across, the resulting x
%		and y values will have errors of order 0.1 pixels for
%		reasonably noise free images.
%
% *********	       READ THE FOLLOWING IMPORTANT CAVEAT!        **********
%		'feature' is capable of finding image features with sub-pixel
%		accuracy, but only if used correctly- that is, if the
%		background is subtracted off properly and the centroid mask
%		is larger than the feature, so that clipping does not occur.
%		It is an EXCELLENT idea when working with new data to plot
%		a histogram of the x-positions mod 1, that is, of the
%		fractional part of x in pixels.  If the resulting histogram
%		is flat, then you're ok, if its strongly peaked, then you're
%		doing something wrong- but probably still getting 'nearest
%		pixel' accuracy.
%
%		For a more quantitative treatment of sub-pixel position
%		resolution see:
%		J.C. Crocker and D.G. Grier, J. Colloid Interface Sci.
%		*179*, 298 (1996).
%
% REVISION HISTORY:
%   ??/??/1992 - David G. Grier
%       * Wrote stats2 version ... a predecessor.
%   10/??/1993 - John C. Crocker
%       * Greatly revised version
%   05/??/1996 - John C. Crocker
%		* Added many keywords
%   06/08/2023 - K Aptowicz (WCU)
%       * Translated to MATLAB 
%       * Tranlation partly based on Blair and Dufresne code
%   08/14/2023 - K Aptowicz (WCU)
%       * Changed how minpix is handled. Now scaled to new value when image
%       is converted pixels values between 0-255 image

%% Reading and setting parameters
% Set default values for optional parameters
default_sep = extent+1;
default_masscut = [];
default_minpix = [];
default_quiet = [];
default_iterate = [];

% Create fields for all optionals inputs
p = inputParser;
% Variables
addParameter(p,'separation',default_sep,@isnumeric)
addParameter(p,'masscut',default_masscut,@isnumeric)
addParameter(p,'min',default_minpix,@isnumeric)

% Keywords
addOptional(p,'quiet', default_quiet)
addOptional(p,'iterate', default_iterate)

% populate optional parameters from inputs
parse(p,varargin{:});
sep = p.Results.separation;
masscut = p.Results.masscut;
minpix = p.Results.min;
quiet = p.Results.quiet;
iterate = p.Results.iterate;
%% *****************************

%% beginning of IDL feature function
extent = floor(extent);
if mod(extent,2) == 0
    disp('Requires an odd extent.  Adding 1...')
    extent = extent + 1;
end

sz = size(image);
ny = sz(1);
nx = sz(2);

%	Put a border around the image to prevent mask out-of-bounds
a = zeros(ny+extent+1,nx+extent+1);
a((extent+1)/2+1:((extent+1)/2)+ny,(extent+1)/2+1:((extent+1)/2)+nx) = image;

%	Finding local maxima
loc = lmx(a,sep,minpix);

if loc(1) == -1
    arr=-1;
    disp('FINDFEATURES: No features found.')
    return
end

[y,x] = ind2sub(size(a),loc);

nmax=numel(loc);
xl = x - fix(extent/2);
xh = xl + extent -1;
m  = zeros(nmax,1);

%	Set up some masks
rsq = rsqd(extent,extent);
t = thetarr( extent );
mask = (rsq <= (single(extent)/2).^2);
%mask2 = make_array( extent , extent , /single, /index ) mod (extent ) + 1.
mask2=(1:1:extent)'*ones(1,extent);
mask2 = (mask2.*mask)';
mask3= (rsq.*mask) + (1./6.);
cen = single(extent+1)/2.;
cmask = cos(2*t).*mask;
smask = sin(2*t).*mask;
cmask(cen,cen) = 0.;
smask(cen,cen) = 0.;

suba = zeros(extent, extent, nmax);
xmask = mask2;
ymask = transpose( mask2 );

yl = y - fix(extent/2);
yh = yl + extent -1;
yscale = 1;
ycen = cen;

%	Estimate the mass
for i=1:nmax
    m(i) = sum(a(yl(i):yh(i),xl(i):xh(i)).*mask,"all");
end

if ~isempty(masscut)
    [w] = find(m > masscut); nmax = numel(w);
    if nmax == 0
        arr=-1;
        disp('FINDFEATURES: No features found!');
        return
    end
    xl = xl(w);
    xh = xh(w);
    yl = yl(w);
    yh = yh(w);
    x = x(w);
    y = y(w);
    m = m(w);
end

if isempty(quiet)
    disp(['FINDFEATURES: ',int2str(nmax) ' features found.'])
end
%	Setup some result arrays
xc = zeros(nmax,1);
yc = zeros(nmax,1);
rg = zeros(nmax,1);
e  = zeros(nmax,1);

%	Calculate feature centers
for i=1:nmax
    xc(i) = sum(a(yl(i):yh(i),xl(i):xh(i)).*xmask,"all");
    yc(i) = sum(a(yl(i):yh(i),xl(i):xh(i)).*ymask,"all");
end

%	Correct for the 'offset' of the centroid masks
xc = xc./m - ((single(extent)+1.)/2.);
yc = (yc./m - ((single(extent)+1.)/2.))/yscale;

%	Iterate any bad initial estimate.
if ~isempty(iterate)
    counter = 0; nbad=-1;
    while (nbad ~= 0) || (counter ~= 10)
        counter = counter + 1;
        w = find(abs(xc) > 0.6); nbadx=numel(w);
        if nbadx > 0
            dx = round( xc(w) );
            xl(w) = xl(w) + dx;
            xh(w) = xh(w) + dx;
            x(w) = x(w) + dx;
        end
        w = find(abs(yc) > 0.6); nbady=numel(w);
        if nbady > 0
            dy = round(yc(w));
            yl(w) = yl(w) + dy;
            yh(w) = yh(w) + dy;
            y(w) = y(w) + dy;
        end

        w = find((abs(xc) > 0.6) | (abs(yc) > 0.6)); nbad = numel(w);
        if nbad > 0 	 % recalculate the centroids for the guys we're iterating
            for i=1:nbad
                m(w(i))  = sum(a(yl(w(i)):yh(w(i)),xl(w(i)):xh(w(i))).*mask,"all");
            end
            for i=1:nbad
                xc(w(i)) = sum(a(yl(w(i)):yh(w(i)),xl(w(i)):xh(w(i))).*xmask,"all");
                yc(w(i)) = sum(a(yl(w(i)):yh(w(i)),xl(w(i)):xh(w(i))).*ymask,"all");
            end
            xc(w) = xc(w)./m(w) - ((single(extent)+1.)/2.);
            yc(w) = (yc(w)./m(w) - ((single(extent)+1.)/2.))/yscale;
        end
    end
end

%	Update the positions and correct for the width of the 'border'
x = x + xc - (extent+1)/2;
y = ( y + yc - (extent+1)/2)*yscale;

%	Construct the subarray
for i=1:nmax
    suba(:,:,i) = fracshift(a(yl(i):yh(i),xl(i):xh(i)),-yc(i),-xc(i));
end

%	Calculate the 'mass'
for i=1:nmax
    m(i) = sum(suba(:,:,i).*mask,"all");
end

%	Calculate radii of gyration squared
for i=1:nmax
    rg(i) = sum(suba(:,:,i).*mask3,"all")/m(i);
end

%	Calculate the 'eccentricity'
for i=1:nmax
    e(i) = sqrt((sum(suba(:,:,i).*cmask,"all").^2) + ...
        (sum(suba(:,:,i).*smask,"all").^2))/(m(i)-suba(ycen,cen,i)+1e-6);
end

arr = [x,y,m,rg,e];
end

%% rsqd function
function r2 = rsqd(w,h)
r2 = zeros(w,h);
xc = single(w-1)/2;
yc = single(h-1)/2;
x = (0:(w-1)) - xc;
x = x.^2;
y = (0:(h-1)) - yc;
y = y.^2;

for j = 1:h
    r2(:,j) = x + y(j);
end
end

%
%	produce a 'theta' mask
%
function theta=thetarr(w)

theta = zeros(w);
xc = single(w-1)/2;
yc = single(w-1)/2;

x = (0:1:w-1) - xc;
x(xc+1) = 1e-5;
y = (0:1:w-1) - yc;

for j = 1:w
    theta(j,:) = atan2(y(j),x);
end
end

%
%	This routine returns the even or odd field of an image
%
% IDL fieldof function deleted ... not using field option
%
%	barrel "shifts" a floating point arr by a fractional pixel amount,
%		by using a 'lego' interpolation technique.
%
function res=fracshift(im,shifty,shiftx)

ipx = fix(shiftx);
ipy = fix(shifty);
fpx = shiftx - ipx;
fpy = shifty - ipy;
if fpx < 0
    fpx=fpx+1; ipx=ipx-1;
end
if fpy < 0
    fpy=fpy+1; ipy=ipy-1;
end

image = im;

imagex  = circshift(image,[ipy,ipx+1]);
imagey  = circshift(image,[ipy+1 ,ipx]);
imagexy = circshift(image,[ipy+1,ipx+1]);
image   = circshift(image,[ipy  ,ipx]);

res   = ( (1. - fpx) * (1. - fpy) * image   ) + ...
    ( (     fpx) * (1. - fpy) * imagex  ) + ...
    ( (1. - fpx) * (     fpy) * imagey  ) + ...
    ( (     fpx) * (     fpy) * imagexy );
end

%
%	John's version of local_max2, which supports the field keyword
%	and is otherwise identical.
%
function r=lmx(image, sep, minpix)

range = fix(sep/2);
a = uint8(rescale(image, 0, 255));
w = round(2*range+1);	% width of sample region
s = rsqd(w,w);			% sample region is circular
good = s <= range.^2;
mask = zeros(w,w,'uint8');
mask(good) = uint8(1);
yrange = range;
b = imdilate(a, mask);	% find local maxima in given range
% but don't include pixels from the
% background which will be too dim

% If not set, determine minpix
if isempty(minpix)
    h = histcounts(a,0:1:256);

    for i = 2:numel(h)
        h(i) = h(i) + h(i-1);
    end
    h = single(h)/max(h);
    minpix = 1;
    while h(minpix) < 0.64
        minpix = minpix + 1;
    end
else
minpix=minpix*255/(max(image(:))-min(image(:)));
end


r = find(a == b & a >= minpix);

% Discard maxima within range of the edge
sz = size(a);
nx = sz(2); ny = sz(1);
[y,x] = ind2sub(sz,r);
x0 = x - range; x1 = x + range;
y0 = y - yrange; y1 = y + yrange;
good = find(x0 >= 1 & x1 < nx+1 & y0 >= 1 & y1 < ny+1);
ngood = numel(good);
if ngood < 1
    r=-1;
else
    r = r(good);
    x = x(good); y = y(good);
    x0 = x0(good); x1 = x1(good);
    y0 = y0(good); y1 = y1(good);
    % There may be some features which get
    % found twice or which have flat peaks
    % and thus produce multiple hits.  Find
    % and clear such spurious points.
    c = 0*a;
    c(r) = a(r);
    center = w * yrange + range+1;

    for i = 1:numel(r)
        b = c(y0(i):y1(i),x0(i):x1(i));
        b =  b.*mask;		% look only in circular region
        [~,location] = max(b,[],"all");
        if location ~= center
            c(y(i),x(i)) = 0;
        end
    end
    % Ideally, the above routine would shrink
    % clusters of points down to their center.
    % As written, this will leave the lower
    % right (?) pixel of a cluster.

    r = find(c ~= 0);		% What's left are valid maxima.
end
end