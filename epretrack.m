function res = epretrack(stk, varargin)

% FUNCTION NAME:
%   epretrack
%
% DESCRIPTION:
%   Routine to analyze time-series stacks or single images to
%   determine where particles are located
%
% INPUT (REQUIRED)
%              stk: image or image stack with particles to track
%                   can also be avi file which is analyzed frame by frame
%
% INPUT (OPTIONAL)
%             bplo: ['y'] Keep all values in output image, even negative values.
%		      bphi: An optional parameter which specifies the
%			        minimum allowable separation between feature
%			        centers. The default value is diameter+1.
%		       dia: Setting this parameter saves runtime by reducing the
%			        runtime wasted on low mass 'noise' features.
%	           sep:
%             mass:
%              min: Set this optional parameter to the minimum allowed
%			        value for the peak brightness of a feature. Useful
%			        for limiting the number of spurious features in
%			        noisy images.
%		     quiet:	['y'] Supress printing of informational messages.
%           prefix:
%            first:
%            fskip:
%
%
% OUTPUT:
%		    pt(:,1): this contains the x centroid positions, in pixels.
%		    pt(:,2): this contains the y centroid positions, in pixels.
%		    pt(:,3): this contains the integrated brightness of the
%			        features.
%		    pt(:,4): this contains the square of the radius of gyration
%			        of the features.
%		    pt(:,5): this contains the eccentricity, which should be
%			        zero for circularly symmetric features and order
%			        one for very elongated images.
%		    pt(:,6): frame number
%
% CALLING SEQUENCE:
%   res = epretrack(a,bplo=1,bphi=11,dia=11,mass=10000)
%
% NOTES :
%   IDL VERSION
%           This code was translated from feature.pro code provided
%           on Eric Weeks' website:
%           https://physics.emory.edu/faculty/weeks/idl/kit/epretrack.pro
%   RESTRICTIONS
%       Read header of 'findfeatures.m'.
%   PROCEDURE
%		Runs bpass and featureon all frames and concatenates the results.
%
% REVISION HISTORY:
%  ppretrack -- Peter's version (begun 7/8/97)
%  jpretrack -- John's version (begun 7/8/98)
%  epretrack -- Eric's version (begun 2/15/05)
%        some minor errors fixed July 2017% improved "single" keyword
%  06/12/2023 - K Aptowicz (WCU)
%       * Translated to MATLAB
%  10/27/2023 - K Aptowicz (WCU)
%       * Added ability to read in .avi files (toolbox?)

%% Reading and setting parameters
% Set default values for optional parameters
default_bplo = [];
default_bphi = [];
default_dia = [];
default_sep = [];
default_mass = [];
default_min = [];
default_first = [];
default_quiet = [];

% Create fields for all optionals inputs
p = inputParser;
% Variables
addParameter(p,'bplo',default_bplo,@isnumeric)
addParameter(p,'bphi',default_bphi,@isnumeric)
addParameter(p,'dia',default_dia,@isnumeric)
addParameter(p,'sep',default_sep,@isnumeric)
addParameter(p,'mass',default_mass,@isnumeric)
addParameter(p,'min',default_min,@isnumeric)

% Keywords
addOptional(p,'first', default_first)
addOptional(p,'quiet', default_quiet)

% populate optional parameters from inputs
parse(p,varargin{:});
bplo = p.Results.bplo;
bphi = p.Results.bphi;
dia = p.Results.dia;
sep = p.Results.sep;
mass = p.Results.mass;
min = p.Results.min;
quiet = p.Results.quiet;
first = p.Results.first;
%% *****************************


%% Populate parameters not defined by user

msg='Defaults:';
if (isempty(bplo))
    bplo = 1; msg=msg+" bplo=1";
end
if (isempty(bphi))
    bphi = 5; msg=msg+" bphi=5";
end
if (isempty(dia))
    dia = 9; msg=msg+" dia=9";
end
if (isempty(sep))
    sep = dia+1;  %  this is what feature uses as a default
    msg=msg+" sep-unset";
end
if (isempty(min))
    min = 0;
    msg=msg+" min-unset";
end
if (isempty(mass))
    mass = 0;
    msg=msg+" mass-unset";
end

if (isempty(quiet))
    if strlength(msg) > 9
        disp(msg)
        disp("starting epretrack...")
    end
end

rep = 1;

if ischar(stk)
    stk = convertCharsToStrings(stk); % Convert the string if character array
end

if isstring(stk)
    disp("analyzing AVI video file frame by frame.")
    v = VideoReader(stk);
    ns = v.NumberOfFrames; % number of frames
    if ns >= 200, rep = 50; end
    if ~isempty(first), ns = 1; end  %handy for a quick looksee...
    res = ones(1,6)*(-1);
    for i = 1:ns
        if ((mod((i),rep) == 0) && isempty(quiet))
            disp(['processing frame ', int2str(i),' out of ',int2str(ns),'....'])
        end
        im = bpass(rgb2gray(read(v,i)),bplo,bphi);
        massTemp=mass;minTemp=min;quietTemp=quiet;
        f = findfeatures(im,dia,sep,mass=massTemp,min=minTemp,quiet=quietTemp,quiet='y');
        nf = numel(f(:,1));
        if (f(1) ~= -1)
            res=[[res];[f,ones(nf,1)*[i]]];
        end
    end

else
    ss=size(stk);
    ns = ss(3);
    if ns >= 200, rep = 50; end
    if ~isempty(first), ns = 1; end  %handy for a quick looksee...
    res = ones(1,6)*(-1);
    for i = 1:ns
        if ((mod((i),rep) == 0) && isempty(quiet))
            disp(['processing frame ', int2str(i),' out of ',int2str(ns),'....'])
        end
        im = bpass(stk(:,:,i),bplo,bphi);
        massTemp=mass;minTemp=min;quietTemp=quiet;
        f = findfeatures(im,dia,sep,mass=massTemp,min=minTemp,quiet=quietTemp,quiet='y');
        nf = numel(f(:,1));
        if (f(1) ~= -1)
            res=[[res];[f,ones(nf,1)*[i]]];
        end
    end
end
% If particles were found, remove blank row
if numel(res(:,1)) > 1
    res(1,:)=[];
end
end

