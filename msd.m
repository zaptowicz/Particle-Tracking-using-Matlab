function data = msd(tracks, varargin)
% FUNCTION NAME:
%   msd
%
% DESCRIPTION:
%   Calculates the mean-squared displacements (MSD) of particle tracks
%
% INPUT (REQUIRED):
%        tracks: (double) is EITHER a filename/wildcard that matches gdf
%                   files from 'track.pro' OR the track output array itself. The
%                   average is calculated for all valid pairs of positions--i.e. it is
%                   overcounted.
%                   NB: if 'tracks' is a wildcard matching several files, the
%	                routine will combine the results into one MSD.
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
%   06/30/2023 - K Aptowicz
%       * Fixed many programming issues to get minN, maxtime, etc. to work
%       correctly. 
%
%% Reading and setting parameters
% Set default values for optional parameters
default_maxtime = [];
default_outfile = [];
default_micperpix = 1;
default_timestep = 1;
default_dim = 2;
default_mydts = [];
default_erode = [];
default_minN = 1000;
default_quiet = [];
default_keepdrift = [];
default_noplot = [];

% Create fields for all optionals inputs
p = inputParser;
% Variables
addParameter(p,'maxtime',default_maxtime,@isnumeric)
addParameter(p,'outfile',default_outfile,@isstring)
addParameter(p,'micperpix',default_micperpix,@isnumeric)
addParameter(p,'timestep',default_timestep,@isnumeric)
addParameter(p,'dim',default_dim,@isnumeric)
addParameter(p,'mydts',default_mydts,@isnumeric)
addParameter(p,'erode',default_erode,@isnumeric)
addParameter(p,'minN',default_minN,@isnumeric)
% Keywords
addOptional(p,'quiet', default_quiet)
addOptional(p,'keepdrift', default_keepdrift)
addOptional(p,'noplot', default_noplot)

% populate optional parameters from inputs
parse(p,varargin{:})
maxtime = p.Results.maxtime;
outfile = p.Results.outfile;
timestep = p.Results.timestep;
dim = p.Results.dim;
mydts = p.Results.mydts;
erode = p.Results.erode;
minN = p.Results.minN;
quiet = p.Results.quiet;
keepdrift = p.Results.keepdrift;
noplot = p.Results.noplot;
micperpix = p.Results.micperpix.*ones(1,dim);  % x and y

if isempty(maxtime)
    maxdt=1e6;
else
    maxdt=maxtime;
end

% if track is file name
if isstring(tracks)
    f = dir(tracks);
    if length(a) == 0
        disp(['No Match: ',tracks])
    end
    nf = length(a)
    filebased = logical(1)
else
    nf=1;  % Number of 'files'
    filebased = logical(0);
end

% Setup time steps to calculate MSD
if isempty(mydts)
    % generate the time partition-- about 10 points per decade
    dt = round(1.15.^[0:1:99]');
    dt = unique(dt);
    w = find(dt <= round(maxdt/timestep));
    if isempty(w)
        disp('Invalid maximum dt!');
    else
        dt = dt(w);
    end
    ndt = length(dt);
else
    if size(mydts,1) == 1
        mydts=mydts';
    end
    dt = mydts;
    ndt = length(dt);
end

% Set up some arrays -- the running sums of the data
data = zeros(ndt,(2*dim)+3);    % x,y,x^2,y^2,r^2,N(r,t) and t

%% MAIN code for building the output file
% do the big nested loop
for i=1:nf
    if filebased
        disp(['Reading file: ', f(i).name])
        trj = read_gdf(f(i)).name
    else
        trj = tracks;
    end
    % convert to um.
    for j=1:dim
        trj(:,j) = trj(:,j)*micperpix(j);
        trj=ltrinterp_vb(trj);				%; fill in gaps in track data
    end
    j = 1; Ncounts=inf;
    while ~((j == ndt+1) | ((isempty(maxtime) & isempty(mydts) & (Ncounts < minN))))
        if filebased && isempty(quiet)
            disp(['Processing file: ',f(i).name,' -- dt equals: ', num2str(dt(j))]);
        elseif isempty(quiet)
            disp(['dt equals: ', num2str(dt(j))])
        end
        data(j,1:(2*dim)+2) = data(j,1:(2*dim)+2) + get_drs(trj,dt(j),erode,dim);
        Ncounts = 2*data(j,2*dim+2)/dt(j);      % KBA - divide by dt for independent measurements
        j=j+1;
    end
end

% truncate the data as needed
if j ~= ndt+1
    w = find((2.*data(:,2*dim+2)./dt >= minN));
    nw=length(w);
    if nw == 0
        disp(['Fewer than ',int2str(minN),' points for all dt, set MAXTIME manually!'])
    end
    data = data(w,:);
    dt = dt(w);
end

% make the running sums into means and variances.
for i = 1:2*dim
    data(:,i) = data(:,i)./data(:,2*dim+2);
end
if ~isempty(keepdrift)
    data(:,dim+1:(2*dim)) = data(:,dim+1:(2*dim)) - (data(:,1:dim)^2);
end
data(:,2*dim+1) = sum(data(:,dim+1:(2*dim)),2);

% put on dt vectors
data(:,2*dim+3) = dt*timestep;

% make up the 'effective N's, expect a factor of two improvement due to overcounting.
data(:,2*dim+2) = 2*data(:,2*dim+2)./dt;
% shift the time from last row to first row to make Eric happy.
data = circshift(data,1,2);
% optionally write out a file
if ~isempty(outfile)
    disp(['Writing file: ',outfile])
    write_gdf(data,outfile);
end

% added by ERW on 8-15-04 just 'cuz I like to plot stuff:
if isempty(noplot)
    plot(data(:,1), data(:,dim+2),'o')
    xlabel('\tau')
    ylabel('msd')
end

end

%% Code for calculating each data point using get_drs
function res = get_drs(t,dt,erode,dim)

% declare some arrays
ncol = length(t(1,:));
res = zeros(1,(2*dim)+2);

% get the rows, we do all for dt le 3, and 'triple count' for bigger dts.
st = circshift(t,-dt,1);
if ~isempty(erode)
    mask = zeros((2*int8(erode))+1,1)+1;
    bad = imdilate((t(:,1) == 0),mask);
    bad = bad | circshift(bad,-dt);		% oink!
    % next few lines rewritten 8-21-00 ERW to be less memory intensive
    w1=find((st(:,ncol-1)-t(:,ncol-1) == dt) & (st(:,ncol)-t(:,ncol)) == 0 & (~bad));
    nw1=length(w1);
    %we require that the point not be in a dilated 'bad'
    %mask (for gaps) nor within erode of an end:
    if (nw1 > 0)
        st2 = circshift(t,erode,1);
        w2=find(st2(w1,ncol) - t(w1,ncol) == 0);
        nw2=length(w2);
        st2 = 0;
        if (nw2 > 0)
            st3 = circshift(t,-dt-erode,1);
            w3 = find((st3(w1(w2),ncol)-t(w1(w2),ncol)) == 0);
            nw3 = length(w3);
            if (nw3 > 0)
                w = w1(w2(w3));
                nw = nw3;
                w1=0; w2=0; w3=0;
            else
            end
        else
            w=w2; w2=0; w1=0;
            nw=nw2;
        end
    else
        w=w1; w1=0;
        nw=nw1;
    end
else
    w = find(((st(:,ncol-1)-t(:,ncol-1)) == dt) & (t(:,1) ~= 0) & (st(:,1) ~= 0) & (st(:,ncol)-t(:,ncol)) == 0);
    nw = length(w);
end

% get the dx's
if (nw > 0)
    dxyz = t(w,1:dim) - st(w,1:dim);
    dxyzsq = dxyz.^2;
    if nw == 1 % DEAL WITH CASE OF VERY SPARSE TRACKS
        disp('WARNING: tracks with very sparse data --> Check ERODE and GOODENOUGH')
        %dxyz = reform(dxyz,dim,1)  %needed?
        %dxyzsq = reform(dxyzsq,dim,1) %needed?
    end
    % do the running totals
    res(1:dim) = sum(dxyz);
    res(dim+1:2*dim) = sum(dxyzsq);
    res(2*dim+2) = nw;
end
end

%% Fill in missing rows of track file; sets x & y values to zero
function tarray = ltrinterp_vb(tarray)
% ;	interpolates gaps in tracked data (time,id only), zeros elsewhere
% ; 	tarray is an array of tracked data
ndat=length(tarray(1,:));
dp=tarray-circshift(tarray,1,1);
w=find((dp(:,ndat) == 0) & (dp(:,ndat-1) ~= 1));
ngood=length(w);
if dp(1,ndat) == 0 %in case tarray is single-particle track
    if ngood > 1
        w = w(2:end);
    end
    ngood = ngood-1;
end
if (ngood >= 1)
    totdt=sum(dp(w,ndat-1))-ngood;
    if totdt > 0  % a subtle case, but it happens!
        dp=0;
        storeres=zeros(totdt,ndat);
        count = 1;
        for i=1:ngood
            dt = tarray(w(i),ndat-1) - tarray(w(i)-1,ndat-1);
            timer = tarray(w(i)-1,ndat-1) + [0:1:dt-2] + 1;
            storeres(count:count+dt-2,ndat) = tarray(w(i),ndat);
            storeres(count:count+dt-2,ndat-1) = timer;
            count = count + dt - 1;
        end
        tarray = [[tarray];[storeres]];
        tarray = sortrows(tarray,[ndat,ndat-1]);
    end
end
end