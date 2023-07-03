function plottr(tarray,varargin)

% goodenough=goodenough, noconnect=noconnect

% IDL keywords not coded: ,tv=tv,, x=x,y=y,dot=dot,timecut=timecut,_extra=eee
%
% FUNCTION NAME:
%   plottr
%
% DESCRIPTION:
%   Plots found tracks as different colors
%
% INPUT (REQUIRED):
%        tarray: (double) 
%
% INPUT (OPTIONAL):
%         maxtime: (double) is the maximum 'dt' to calculate, in seconds.
%                   If this is not set, MSD computes out in dt until fewer than
%                   1000 independent measurements are being averaged.  Setting
%                   'maxtime' or 'mydts' overrides this limitation, of course.
%         outfile: (string) will save the result to a gdf file named
%                  'outfile', which is useful for batch file operation.
%       micperpix: (double) set to the pixel size to convert the input
%                  data file to microns if that has not been done already.
%                  If it is single number, the magnification is assumed to
%                  be isotropic, otherwise, provide a dim-dimensional vector.
%        timestep: (double) the number of seconds between frames (e.g.
%                  1/60.)
%             dim: (double) the spatial dimensionality of the data (default 2).
%           mydts: (array) allows the user to supply a vector of (frame
%                  number) dt's if he/she doesn't want my log-spaced points
%                  in time.
%           erode: (double) drops all data within 'erode' time steps of an
%                  end or gap in a trajectory, to eliminate mistaken track
%                  contributions. Very useful-- usually erode=2 or 3 works
%                  pretty well.
%            minN: (double) sets minimum number of independent measurements;
%                  only useful maxtime is not used
%           quiet: set to 'y' to NOT display messages
%       keepdrift: set to 'y' to NOT remove <x>^2
%          noplot: set to 'y' to NOT plot anything at the end
%
% OUTPUT:
%             data: (double)tracked data array where
%                   if dim =2 then
%                   [time,<x>,<y>,<x^2>,<y^2>,<x^2 + y^2>,<N>],
%                   or if dim=3 then
%                   [time,<x>,<y>,<z>,<x^2>,<y^2>,<z^2>,<x^2 + y^2 + z^2>,<N>],
%                   and N is the approx. independent number of data points
%                   in the average.
%	                NB: for all variances, <x^2> has had <x>^2 subtracted off.
%
% CALLING SEQUENCE:
%   m = msd(t)
%   m = msd(t,micperpix=0.137, dim=2, timestep=1/60,maxtime=1,minN=100,quiet='y',noplot='y');
%
% NOTES :
%   IDL Version - This code was translated from msd.pro code provided on
%           Eric Weeks' website:
%           https://physics.emory.edu/faculty/weeks/idl/yykit/msd.pro
%   Variances - The variances are presumably what you are interested in. The
%           means are there merely to give you an idea of the average 'drift'
%           in the data (and because <x>^2 has to be subtracted from the
%           variances).  'N' can be used for error estimation: the error
%           in the variance < x(dt)^2 > ~ < x(dt)^2 >/sqrt( N(dt) ).
%
% REVISION HISTORY:
%   06/??/1999 - John Crocker
%       * Wrote orginal version at U. Penn
%   05/??/2002 - Victor Breedveld
%       * Included minN keyword
%   12/??/2003 - Victor Breedveld
%       * Fixed sorting bug in LTRINTERP.PRO under DOS (Windows)
%       * interpolation routine now named LTRINTERP_VB
%   06/??/2004 - Victor Breedveld
%       * Adapted to work for single particle tracks!!!
%   07/??/2004 - Victor Breedveld
%       * Adapted to work with results of poor tracking!!!
%   10/23/2022 - K Aptowicz
%       * Translated to MATLAB
%
%% Reading and setting parameters
% Set default values for optional parameters


%% Reading and setting parameters
% Set default values for optional parameters
default_ID = [];
default_goodenough = [];

% Create fields for all optionals inputs
p = inputParser;
% Variables
addParameter(p,'ID',default_ID,@isnumeric)
addParameter(p,'goodenough',default_goodenough,@isnumeric)
% Keywords
%addOptional(p,'quiet', default_quiet)


% populate optional parameters from inputs
parse(p,varargin{:})
ID = p.Results.ID;
goodenough = p.Results.goodenough;

if ~isempty(ID)
    w=tarray(:,end)==ID;
    tarray=tarray(w,:);
end

% Save only tracks of length equal to or greater than good enough
if ~isempty(goodenough)
    ndat=numel(tarray(1,:));
    [~,u,~] = unique(tarray(:,ndat));
    u=[u;length(tarray(:,ndat))+1];
    track_length = u(2:end)-u(1:end-1);
    w=find(track_length >= goodenough);
    temp=zeros(1,ndat);
    for i=1:numel(w)
        temp=vertcat(temp,tarray(u(w(i)):u(w(i)+1)-1,:));
    end
    tarray=temp(2:end,:);
end

% Plot each track
ndat=numel(tarray(1,:));
[~,u,~] = unique(tarray(:,ndat));
ntracks=numel(u);
u=[u;length(tarray(:,ndat))+1];
for i=1:ntracks
    x=tarray(u(i):u(i+1)-1,1);
    y=tarray(u(i):u(i+1)-1,2);
    plot(x,y,'-')
    hold on
end
hold off
axis equal tight
end