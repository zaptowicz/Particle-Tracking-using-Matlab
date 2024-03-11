function res = epretrack(stk, varargin)

% FUNCTION NAME:
%   epretrack
%
% DESCRIPTION:
%   Routine to analyze time-series stacks or single images to
%   determine where particles are located
%
% INPUT (REQUIRED)
%              stk: Three options are currently coded: 
%                   (1) stk is an array (2D if single image or 3D if image stack)
%                   (2) stk is the filename of an avi file to read in 
%                   (3) stk is a VIF file (EPIX).  If this is used, also
%                   need to input parameter into the variable VIF.Info.
%                   Discussed below. 
%
% INPUT (OPTIONAL)
%             ***** b = bpass(image, bplo, bphi) *****
%             bplo: High spatial frequnecy cutoff for bpass. Often 1 pixel for
%                   pixel-to-pixel noise. Sticking with the naming in IDL
%                   epretrack calling this 'lo'.
%		      bphi: Low spatial frequency cutoff for bpass. Eliminates slow 
%                   varations in the image. Often a little larger than 
%                   feature size.Sticking with the naming in IDL
%                   epretrack calling this 'hi'.
%
%             ****  findfeatures(image, extent, varargin) *****
%		    extent: a parameter which should be a little greater than
%			        the diameter of the largest features in the image.
%                   Extent MUST BE ODD valued.
%		separation: An optional parameter which specifies the
%			        minimum allowable separation between feature
%             mass: Setting this parameter saves runtime by reducing the
%			        runtime wasted on low mass 'noise' features.
%              min: Set this optional parameter to the minimum allowed
%			        value for the peak brightness of a feature. Useful
%			        for limiting the number of spurious features in
%			        noisy images.
%		     quiet:	['y'] Supress printing of informational messages.
%           prefix:
%            first:
%            fskip:
%         VIF_Info: structure containing information needed to read in VIF files.
%                   Example: 
%                       VIF_input.pix_w = 658;
%                       VIF_input.pix_h = 494;
%                       VIF_input.frame_N = 10;
%                       VIF_input.byte_offset=72;
%                       VIF_input.byte_spacing = 136;
%                       VIF_input.bit_10_unpacked='y';
%                    >> pt_all = epretrack('videofile.vif', VIF_Info=VIF_input);
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
%   res = epretrack(a,bplo=1,bphi=11,extent=11,mass=10000)
%   pt_all = epretrack('videofile.vif',bplo=1,bphi=9,extent=11,min=50,VIF_Info=VIF_input);
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
%  03/11/2024 - K Aptowicz (WCU)
%       * Made some edits to variables names and cleaned up header. 

%% Reading and setting parameters
% Set default values for optional parameters
default_bplo= [];
default_bphi = [];
default_extent = [];
default_separation = [];
default_mass = [];
default_min = [];
default_first = [];
default_quiet = [];
default_VIF_Info = [];

% Create fields for all optionals inputs
p = inputParser;
% Variables
addParameter(p,'bplo',default_bplo,@isnumeric)
addParameter(p,'bphi',default_bphi,@isnumeric)
addParameter(p,'extent',default_extent,@isnumeric)
addParameter(p,'separation',default_separation,@isnumeric)
addParameter(p,'mass',default_mass,@isnumeric)
addParameter(p,'min',default_min,@isnumeric)
addParameter(p,'VIF_Info',default_VIF_Info)

% Keywords
addOptional(p,'first', default_first)
addOptional(p,'quiet', default_quiet)

% populate optional parameters from inputs
parse(p,varargin{:});
bplo= p.Results.bplo;
bphi = p.Results.bphi;
extent = p.Results.extent;
separation = p.Results.separation;
mass = p.Results.mass;
min = p.Results.min;
quiet = p.Results.quiet;
first = p.Results.first;
VIF_Info = p.Results.VIF_Info;

%% *****************************

%% Populate parameters not defined by user
% populare VIF_Info if needed
if ~isempty(VIF_Info)
    if ~(isfield(VIF_Info,'pix_w') && isfield(VIF_Info,'pix_h'))
        disp('Need size of VIF video frame (pix_w and pix_h) in structure VIF_Info')
    end
    VIF_input.frame_N = 10;
    if ~(isfield(VIF_Info,'frame_N'))
        disp('VIF_Info.frame_N not set, so setting to 10')
        VIF_Info.frame_N=10; 
    end
    if ~(isfield(VIF_Info,'byte_offset')), VIF_Info.byte_offset=[]; end
    if ~(isfield(VIF_Info,'byte_spacing')), VIF_Info.byte_spacing=[]; end
    if ~(isfield(VIF_Info,'frame_skip')), VIF_Info.frame_skip=0; end
    if ~(isfield(VIF_Info,'bit_10_packed')), VIF_Info.bit_10_packed=[]; end
    if ~(isfield(VIF_Info,'bit_10_unpacked')), VIF_Info.bit_10_unpacked=[]; end
    if ~(isfield(VIF_Info,'bit_8')), VIF_Info.bit_8=[]; end
    if ~(isfield(VIF_Info,'byte_offset')), VIF_Info.byte_offset=[]; end
end

msg='Defaults:';
if (isempty(bplo))
    bplo= 1; msg=msg+" bplo=1";
end
if (isempty(bphi))
    bphi = 5; msg=msg+" bphi=5";
end
if (isempty(extent))
    extent = 9; msg=msg+" extent=9";
end
if (isempty(separation))
    separation = extent+1;  %  this is what feature uses as a default
    msg=msg+" separation-unset";
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
    stk = string(stk); % Convert the string if character array
end

% Analysis of AVI videos
if isstring(stk)
    if endsWith(stk,'avi','IgnoreCase',true)
        disp("Analyzing AVI video file frame by frame. Converting to grayscale if RGB")
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
            f = findfeatures(im,extent,separation,mass=massTemp,min=minTemp,quiet=quietTemp,quiet='y');
            nf = numel(f(:,1));
            if (f(1) ~= -1)
                res=[[res];[f,ones(nf,1)*[i]]];
            end
        end
    elseif endsWith(stk,'vif','IgnoreCase',true)
        disp("Analyzing VIF video file frame by frame.")
        ns = VIF_Info.frame_N; % number of frames
        if ns >= 200, rep = 50; end
        if ~isempty(first), ns = 1; end  %handy for a quick looksee...
        res = ones(1,6)*(-1);
        for i = 1:ns
            if ((mod((i),rep) == 0) && isempty(quiet))
                disp(['processing frame ', int2str(i),' out of ',int2str(ns),'....'])
            end
            im = read_vif(stk,VIF_Info.pix_w,VIF_Info.pix_h, ...
                frame_start=i,frame_N=1, ...
                byte_offset=VIF_Info.byte_offset, ...
                byte_spacing = VIF_Info.byte_spacing, ...
                frame_skip = VIF_Info.frame_skip, ...
                bit_10_packed = VIF_Info.bit_10_packed, ...
                bit_10_unpacked = VIF_Info.bit_10_unpacked, ...
                bit_10_unpacked = VIF_Info.bit_8);
            im = bpass(im,bplo,bphi);
            massTemp=mass;minTemp=min;quietTemp=quiet;
            f = findfeatures(im,extent,separation,mass=massTemp,min=minTemp,quiet=quietTemp,quiet='y');
            nf = numel(f(:,1));
            if (f(1) ~= -1)
                res=[[res];[f,ones(nf,1)*[i]]];
            end
        end
    else
        disp('Cannot process: Expecting .vif or .avi file formats')
    end
else

    % Analysis of arrays
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
        f = findfeatures(im,extent,separation,mass=massTemp,min=minTemp,quiet=quietTemp,quiet='y');
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

