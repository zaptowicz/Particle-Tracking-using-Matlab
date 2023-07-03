function output=fover2d(image,points,varargin)

% FUNCTION NAME:
%   fover2d
%
% DESCRIPTION:
%   Overlay points onto a 2d image.
%
% INPUT (REQUIRED)
%           image:	2D data set onto which the points should be overlaid
%		   points: (2,npoints) array of overlay points ... ignores columns
%		           above 2.
%
% INPUT (OPTIONAL)
%		    radius: sets size of circles
%              big: ['y'] doubles the picture in size in each direction
%            nodot: ['y'] turns off the black dot at the center of each circle
%           circle: ['y'] draws a circle around each point, rather than a disk
%            track: ['y'] re-orders a track array so it can be used as a pretrack
%                   array
%
% OUTPUT:
%		    output: returns an image ready for 'imagesc'.
%
% CALLING SEQUENCE:
%   fo=fover2d(a0,f,big='y')%
%
% NOTES :
%   IDL VERSION
%           This code was translated from feature.pro code provided
%           on Eric Weeks' website:
%           https://physics.emory.edu/faculty/weeks/idl/kit/fover2d.pro
% PROCEDURE:
%       Rescale the image to leave some color table indices free.
%		Make the rest of the color table into a grey ramp and
%		turn the 3 highest indices into 3 new colors.  Plot a
%       disk at each particle position.
% MODIFICATION HISTORY:
%       Lookup table code taken from David G. Grier's f_overlay
%       routine.  Mostly written by Eric Weeks in summer '98.
%       Ability to handle movies added 2-20-99.
%       Speed up 4-99
%       15 May 2006 Gianguido Cianci: added /track so it can use
%       positions from a track array.
%   06/08/2023 - K Aptowicz (WCU)
%       * Translated to MATLAB
%

%% Reading and setting parameters
% Set default values for optional parameters
default_radius = 5;
default_big = [];
default_nodot = [];
default_circle = [];
default_track = [];

% Create fields for all optionals inputs
p = inputParser;
% Variables
addParameter(p,'radius',default_radius,@isnumeric)

% Keywords
addOptional(p,'big', default_big)
addOptional(p,'nodot', default_nodot)
addOptional(p,'circle', default_circle)
addOptional(p,'track', default_track)

% populate optional parameters from inputs
parse(p,varargin{:});
radius = p.Results.radius;
big = p.Results.big;
nodot = p.Results.nodot;
circle = p.Results.circle;
track = p.Results.track;
%% *****************************

if ~isempty(track)
    ncols = numel(points(1,:));
    [~,s] = sort(points(:,ncols-1));         %sort times
    pnts = points(1:ncols-1,s);          %do not use id numbers
else
    pnts = points;
end

% KBA - ignored the color table stuff in the IDL code.
nc=255;

output = rescale(image,0,255);
s = size(output);
if (numel(s) == 3)
    % 3-D array
    nel=numel(pnts(1,:));
    t = pnts(:,nel);
    nzz=numel(output(1,1,:));
end
if isempty(big)
    x=round(pnts(:,2));
    y=round(pnts(:,1));
else
    if (numel(s) == 2)
        output=imresize(output,2);         %        2-D array
    else
        output = imresize3(output,[s(1)*2,s(2)*2,s(3)]); %   3-D array
    end
    x = round(pnts(:,2))*2; 
    y = round(pnts(:,1))*2;
    radius=2*radius;
end

x2=numel(output(:,1,1));    %		3rd array index just for safety
y2=numel(output(1,:,1));
xmax=max(x); xmin=min(x);
ymax=max(y); ymin=min(y);
if (xmin < 1 || ymin < 1 || xmax > x2 || ymax > y2)
    disp('FOVER2D: points outside of picture ... not plotted.')
    w=find((x-radius >= 1) & (x+radius <= x2))
    x=x(w); y=y(w);
    w=find((y-radius >= 1) & (y+radius <= y2))
    x=x(w); y=y(w);
end

extent=radius*2+1;
blob=fo_circ(zeros(extent,extent),circle);
blob = (blob > 0).*(nc);

for i = 1:numel(x)
    minx = x(i)-radius; bminx=1;
    miny = y(i)-radius; bminy=1;
    maxx = x(i)+radius; bmaxx=extent;
    maxy = y(i)+radius; bmaxy=extent;

    % fix if disc/circle goes outside of image
    if minx < 1; bminx=2-minx; minx=1; end
    if miny < 1; bminy=2-miny; miny=1; end
    if maxx > x2; bmaxx=extent-(maxx-x2); maxx=x2; end
    if maxy > y2; bmaxy=extent-(maxy-y2); maxy=y2; end

    if (numel(s) == 2)        % 2-D image
        output(minx:maxx,miny:maxy) = ...
            max(blob(bminx:bmaxx,bminy:bmaxy),output(minx:maxx,miny:maxy));
        if isempty(nodot)
            if (x(i) > 1) && (x(i) < (x2-1))
                output(x(i)-1:x(i)+1,y(i))=0;
            end
            if (y(i) > 1) && (y(i) < (y2-1))
                output(x(i),y(i)-1:y(i)+1)=0;
            end
        end
    else
        % 3-D image
        tnow=t(i);
        if (tnow < nzz)
            output(minx:maxx,miny:maxy,tnow) = ...
                max(blob(bminx:bmaxx,bminy:bmaxy),output(minx:maxx,miny:maxy,tnow));
            if isempty(nodot)
                if (x(i) > 1) & (x(i) < (x2-1))
                    output(x(i)-1:x(i)+1,y(i),tnow)=0;
                end
                if (y(i) > 1) & (y(i) < (y2-1))
                    output(x(i),y(i)-1:y(i)+1,tnow)=0;
                end
            end
        end
    end
end
if (numel(s) == 2)        % 2-D image
    colormap("gray")
    imagesc(output); axis equal tight
else
    colormap("gray")
    implay(output/255,10)
end
end
%First a utility function....

% circarray.pro,  started 6-22-98 by ERW
%   shortened into fo_circ 2-20-99 by ERW
%
function result=fo_circ(array,circle)
% returns an array, size equal to "array" variable, with value 1
% everywhere within a circle of diameter of the array size.  Circle
% is at center of array.
%
% 'radius' sets a radius different from the default radius (half the array size)

s=size(array);
result=array*0;

sx=s(1); sy=s(2);
minsize = sx;
cx=(sx+1)*0.5; cy=(sy+1)*0.5;
irad=(minsize-1)/2;

jrad=(irad-1.2)*(irad-1.2);
irad = irad*irad;
if (~isempty(circle))
    for j=1:sy
        rad2 = (cy-j)*(cy-j);
        for k=1:sx
            rad3 = rad2 + (cx-k)*(cx-k);
            result(k,j) = ((rad3 < irad) & (rad3 > jrad));
        end
    end
else
    for j=1:sy
        rad2 = (cy-j)*(cy-j);
        for k=1:sx
            rad3 = rad2 + (cx-k)*(cx-k);
            result(k,j) = (rad3 < irad);
        end
    end
end
end