function res = track(xyzs, maxdisp, varargin)

% FUNCTION NAME:
%   track
%
% DESCRIPTION:
%	Constructs n-dimensional trajectories from a scrambled list of
%	particle coordinates determined at discrete times (e.g. in
%	consecutive video frames).
%
% INPUT (REQUIRED)
%	  positionlist: an array listing the scrambled coordinates and data
%		            of the different particles at different times, such that:
%	                positionlist[:,1:d] contains the d coordinates and
%		            data for all the particles, at the different times.
%	                positionlist[:,d+1] contains the time t that the position
%		            was determined (e.g. frame number). These values must
%                   be monotonically increasing and uniformly gridded in time.
%	       maxdisp: an estimate of the maximum distance that a particle
%		            would move in a single time interval.(see Restrictions)
%
% INPUT (OPTIONAL)
%	        inipos: if the user wants to track only a subset of the
%		            particles, this argument is set to an array (d,n)
%		            which contains the d dimensional initial positions of
%		            the n particles to be tracked. Other 'new' particles
%		            will then NOT be added.
%	        memory: this is the number of time steps that a particle can be
%		            'lost' and then recovered again.  If the particle reappears
%		            after this number of frames has elapsed, it will be
%		            tracked as a new particle. The default setting is zero.
%		            this is useful if particles occasionally 'drop out' of
%		            the data.
%	           dim: if the user would like to unscramble non-coordinate data
%		            for the particles (e.g. apparent radius of gyration for
%		            the particle images), then positionlist should
%		            contain the position data in positionlist[:,1:dim]
%		            and the extra data in positionlist[:,dim+1:d]. It is then
%		            necessary to set dim equal to the dimensionality of the
%		            coordinate data to so that the track knows to ignore the
%		            non-coordinate data in the construction of the
%		            trajectories. The default value is two.
%	       verbose: [y] set this keyword for more informational messages.
%	    goodenough: set this keyword to eliminate all trajectories with
%		            fewer than goodenough valid positions.  This is useful
%		            for eliminating very short, mostly 'lost' trajectories
%		            due to blinking 'noise' particles in the data stream.
%            quiet: [y] use to not print any messages
%
% OUTPUT:
%	        result: a list containing the original data rows sorted
%		            into a series of trajectories.  To the original input
%		            data structure there is appended an additional column
%		            containing a unique 'id number' for each identified
%		            particle trajectory.  The result array is sorted so
%		            rows with corresponding id numbers are in contiguous
%		            blocks, with the time variable a monotonically
%		            increasing function inside each block.  For example:
%
%	 	            For the input data structure (positionlist):
%           			    (x)	         (y)	      (t)
%       	     	pos = 3.60000      5.00000      0.00000
%       		          15.1000      22.6000      0.00000
%       		          4.10000      5.50000      1.00000
%       		          15.9000      20.7000      2.00000
%       		          6.20000      4.30000      2.00000
%
% MATLAB>> res = track(pos,5,mem=2)
%
%		            track will return the result 'res'
%	                		(x)	         (y)	      (t) 	       (id)
%           		res = 3.60000      5.00000      0.00000      0.00000
%		                  4.10000      5.50000      1.00000      0.00000
%		                  6.20000      4.30000      2.00000      0.00000
%		                  15.1000      22.6000      0.00000      1.00000
%		                  15.9000      20.7000      2.00000      1.00000
%
%		NB: for t=1 in the example above, one particle temporarily
%		vanished.  As a result, the trajectory id=1 has one time
%		missing, i.e. particle loss can cause time gaps to occur
%		in the corresponding trajectory list. In contrast:
%
% MATLAB>> res = track(pos,5)
%
%		            track will return the result 'res'
%		                 	(x)	          (y)	       (t) 	       (id)
%	             	res = 15.1000      22.6000      0.00000      0.00000
%      		              3.60000      5.00000      0.00000      1.00000
%  		                  4.10000      5.50000      1.00000      1.00000
% 		                  6.20000      4.30000      2.00000      1.00000
% 		                  15.9000      20.7000      2.00000      2.00000
%
%		where the reappeared 'particle' will be labelled as new
%		rather than as a continuation of an old particle since
%		mem=0.  It is up to the user to decide what setting of
%		'mem' will yield the highest fidelity tracking.
%
% CALLING SEQUENCE:
%   res = track(pos,5,mem=2)
%
% NOTES:
%   IDL VERSION
%          *This code was translated from feature.pro code provided
%           on Eric Weeks' website:
%           https://physics.emory.edu/faculty/weeks/idl/kit/track.pro
%          *Produces informational messages.  Can be memory intensive for
% 	        extremely large data sets.
%   RESTRICTIONS
%	        maxdisp should be set to a value somewhat less than the mean
%	        spacing between the particles. As maxdisp approaches the mean
%	        spacing the runtime will increase significantly. The function
%	        will produce an error message: "Excessive Combinatorics!" if
%	        the run time would be too long, and the user should respond
%	        by re-executing the function with a smaller value of maxdisp.
%	        Obviously, if the particles being tracked are frequently moving
%	        as much as their mean separation in a single time step, this
%	        function will not return acceptable trajectories.
%   PROCEDURE
%	     Given the positions for n particles at time t(i), and m possible
%	     new positions at time t(i+1), this function considers all possible
%	     identifications of the n old positions with the m new positions,
%	     and chooses that identification which results in the minimal total
%	     squared displacement. Those identifications which don't associate
%	     a new position within maxdisp of an old position ( particle loss )
%	     penalize the total squared displacement by maxdisp^2. For non-
%	     interacting Brownian particles with the same diffusivity, this
%	     algorithm will produce the most probable set of identifications
%	     (provided maxdisp >> RMS displacement between frames ).
%	     In practice it works reasonably well for systems with oscillatory,
%	     ballistic, correlated and random hopping motion, so long as single
%	     time step displacements are reasonably small.  NB: multidimensional
%	     functionality is intended to facilitate tracking when additional
%	     information regarding target identity is available (e.g. size or
%	     color).  At present, this information should be rescaled by the
%	     user to have a comparable or smaller (measurement) variance than
%	     the spatial displacements.
%
% REVISION HISTORY:
%	 2/93 Written by John C. Crocker, University of Chicago (JFI).
%	 7/93 JCC fixed bug causing particle loss and improved performance
%		for large numbers of (>100) particles.
%	11/93 JCC improved speed and memory performance for large
%		numbers of (>1000) particles (added subnetwork code).
%	 3/94 JCC optimized run time for trivial bonds and d<7. (Added
%		d-dimensional raster metric code.)
%	 8/94 JCC added functionality to unscramble non-position data
%		along with position data.
%	 9/94 JCC rewrote subnetwork code and wrote new, more efficient
%		permutation code.
%	 5/95 JCC debugged subnetwork and excessive combinatorics code.
%	12/95 JCC added memory keyword, and enabled the tracking of
%		newly appeared particles.
%	 3/96 JCC made inipos a keyword, and disabled the adding of 'new'
%		particles when inipos was set.
%	 3/97 JCC added 'add' keyword, since Chicago users didn't like
%		having particle addition be the default.
%	 9/97 JCC added 'goodenough' keyword to improve memory efficiency
%		when using the 'add' keyword and to filter out bad tracks.
%       10/97 JCC streamlined data structure to speed runtime for >200
%               timesteps.  Changed 'quiet' keyword to 'verbose'. Made
%               time labelling more flexible (uniform and sorted is ok).
%	 9/98 JCC switched trajectory data structure to a 'list' form,
%		resolving memory issue for large, noisy datasets.
%  09/17/1998 Eric Weeks, Emory University, luberize code.
%	 2/99 JCC added Eric Weeks's 'uberize' code to post-facto
%		rationalize the particle id numbers, removed 'add' keyword.
%  03/30/2010 David G. Grier, New York University: Modernized array
%    notation.  Small code modernizations.  Formatting.
%    Moved luberize code into main procedure.  Replaced UNQ with
%    IDL system routine,
%  06/12/2023 - K Aptowicz (WCU)
%       * Translated to MATLAB
%  07/21/2023 - K Aptowicz (WCU)
%       * Fixed issues caused by SUBREF.M when an array is indexed with an
%       array. Comments in code (search SUBREF)
%  05/06/2024 - K Aptowicz (WCU)
%       * Fixed multiple bugs that arose with dense packings. Output now
%       appears to match perfectly with IDL version.
%
%	This code 'track.pro' is copyright 1999, by John C. Crocker.
%	It should be considered 'freeware'- and may be distributed freely
%	(outside of the military-industrial complex) in its original form
%	when properly attributed.
%
% LICENSE:
%    This program is free software% you can redistribute it and/or
%    modify it under the terms of the GNU General Public License as
%    published by the Free Software Foundation% either version 2 of the
%    License, or (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY% without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program% if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
%    02111-1307 USA
%
%    If the Internet and WWW are still functional when you are using
%    this, you should be able to access the GPL here:
%    http://www.gnu.org/copyleft/gpl.html

%% Reading and setting parameters
% Set default values for optional parameters
default_inipos = [];
default_memory = [];
default_dim = [];
default_goodenough = [];
default_verbose = [];
default_quiet = [];

% Create fields for all optionals inputs
p = inputParser;
% Variables
addParameter(p,'inipos',default_inipos,@isnumeric)
addParameter(p,'memory',default_memory,@isnumeric)
addParameter(p,'dim',default_dim,@isnumeric)
addParameter(p,'goodenough',default_goodenough,@isnumeric)

% Keywords
addOptional(p,'verbose', default_verbose)
addOptional(p,'quiet', default_quiet)

% populate optional parameters from inputs
parse(p,varargin{:})
inipos = p.Results.inipos;
memory_b = p.Results.memory;
dim = p.Results.dim;
goodenough = p.Results.goodenough;
verbose = p.Results.verbose;
quiet = p.Results.quiet;
%% *****************************

% calculate parameters about time
dd =numel(xyzs(1,:));               % # of column for time
t = xyzs(:,dd);                     % time column
st = circshift(t,1);        
st = t(2:end)-st(2:end);            % change in frame number in adjacent rows
w = find(st>0); z = length(w);      % z is the number of nonzero time jumps 
z = z + 1;                          % total number of frames (# of frame jumps + 1)

% check the input time vector is ok, i.e. sorted and uniform
if  sum(st(:)<0)
    disp('****************************************')
    disp('ERROR: Some time vectors are decreasing!')
    disp('****************************************')
    return
end

% Fix for missing frames (KBA)
% Original track program wasn't written to handle missing frames. An easy
% fix is to populate the pt file with dummy particles for the missing frames
% for dummy particle with positions x = y = 0, then to delete these particles 
% in the tracked file before outputting. 

if (sum(st(w) ~= st(w(1))) ~= 0)
    disp(' *** Warning - frames are missing. ***')
    disp(' Adding dummy particles in missing frames: ')
    disp(['   Assuming frames are an arithmetic sequential with a '])
    disp(['   common difference of one  (8, 9, 10, 11, ...'])
    disp(' *************************************')

    % KBA: Fix for videos with missing frames
    tmiss = setdiff([min(t):1:max(t)], t);  % times for missing frames. 
    pt_pad = zeros(length(tmiss),dd);     % make dummy pt file
    pt_pad(:,dd) = tmiss;                 % add time stamps of missing frames
    pt_pad(:,1:2) = 0;                      % set x and y positions to zero
    xyzs=sortrows([xyzs;pt_pad],dd);      % add in original pt file and sort by time

    % Update time variables
    t = xyzs(:,dd);
    st = circshift(t,1);
    st = t(2:end)-st(2:end);
    w = find(st>0); z = length(w);
    z = z + 1;                      % number of time steps
    dummy=1;                        % track if dummy particles were inserted
end

if isempty(dim)
    dim = min(2,dd-1);
    if isempty(quiet)
        disp(['Setting dim = ',int2str(dim)+' by default']);
    end
end

if isempty(memory_b); memory_b = 0; end

if z == 0
    disp('****************************************')
    disp('ERROR: All data have the same time stamp!')
    disp('****************************************')
    return
end

% partition the input data by unique times
% res is the index of the beginning of unique time in the track
[~,res,~] = unique(t);
res = [res', length(t)+1]'; % KBA: IDL has the last index twice, not sure why


% get the initial positions
ngood = res(2) - res(1);  % Number of particles in first frame.
eyes = 1:ngood;

if  ~isempty(inipos)
    pos = inipos(:,1:dim);
    istart = 1;
    n = numel(pos(:,0));
else
    pos = xyzs(eyes,1:dim);
    istart = 2;                  %we don't need to track t=1.
    n = ngood;
end

% How long are the 'working' copies of the data?
% More particles in a first frame, shorter 'working' copies in time
zspan = 50;
if n > 200, zspan = 20; end
if n > 500, zspan = 10; end

resx = zeros(zspan,n) - 1;  % Size = [Time of "working" copy, N particle in first frame];
bigresx = zeros(z,n) - 1;   % Size = [Total number of frames, N particle in first frame];
mem = zeros(n,1);           % Size = [N particle in first frame,1];
uniqid = 1:n;
maxid = n;
olist = [0.,0.];

if  ~isempty(goodenough)
    dumphash = zeros(n,1);
    nvalid = zeros(n,1);    % KBA changed from ones to zeros
end

% we may not need to track the first time step!
if isempty(inipos)
    resx(1,:) = eyes;
    if  ~isempty(goodenough), nvalid=nvalid+1; end
end

% set up some nice constants
maxdisq = maxdisp^2;
verbose =  ~isempty(verbose);

%Use fancy code for large n, small d
notnsqrd = (sqrt(n*ngood) >= 200) && (dim < 7);

if notnsqrd
    %   construct the vertices of a 3x3x3... d-dimensional hypercube
    % for dim=2 ... 3x3 grid; for dim = 3 ... 3x3x3 grid; etc
    cube = zeros(3^dim,dim);
    for d=0:dim-1
        numb = 0;
        for j=0:(3^d):(3^dim)-1
            cube(j+1:j+(3^(d)),d+1) = numb;
            numb = mod(numb+1,3);
        end
    end

    %   calculate a blocksize which may be greater than maxdisp, but which
    %   keeps nblocks reasonably small.
    %   KBA: Not sure what this does.

    volume = 1;
    for d = 0:dim-1
        minn = min(xyzs(w,d+1));
        maxx = max(xyzs(w,d+1));
        volume = volume * (maxx-minn);
    end
    blocksize = max( [maxdisp,((volume)/(20*ngood))^(1.0/dim)] );
end

% Setup some variables for displaying updates
rep = 1;
remove = 0;

%   Start the main loop over the frames.
for i = istart:z
    ispan = mod(i-1,zspan)+1;

    %   Get the new particle positions.
    m = res(i+1) - res(i);      % Number of particles in new frame
    eyes = (1:m) + res(i)-1;    % New possible unique particle IDs

    if m  > 0
        xyi = xyzs(eyes,1:dim);
        found = zeros(m,1);
        %%   THE TRIVIAL BOND CODE BEGINS
        %   Trivial bond code has two approaches:
        %       notnsqrd == 1
        %       Looks at each frame (i and i-1) and determines if there are
        %       any particles that are only 'bonded' to each other. These
        %       are trivial cases and identified.
        %
        %       notnsqrd == 0
        %       Not sure how this one works, although it is called 'fancy'

        if notnsqrd
            %   Use the raster metric code to do trivial bonds

            %   construct "s", a one dimensional parameterization of the space
            %   ( which consists of the d-dimensional raster scan of the volume.)
            abi = floor(xyi./blocksize);      % Frame i positions (round towards zero)
            abpos = floor(pos./blocksize);    % Frame i-1 positions
            si = zeros(m,1);
            spos = zeros(n,1);
            dimm = zeros(dim,1);
            coff = 1.;

            for j=1:dim
                minn = min([abi(:,j);abpos(:,j)]);  % minimum reduced value of dimension dim of frames i and i-1.
                maxx = max([abi(:,j);abpos(:,j)]);  % maximum reduced value of dimension dim of frames i and i-1.
                abi(:,j) = abi(:,j) - minn;         % subtract minimum value from frame i
                abpos(:,j) = abpos(:,j) - minn;     % subtract minimum value from frame i-1
                dimm(j) = maxx-minn + 1;            % Find max size of dimension dim.
                si = si + abi(:,j).*coff;           % For each particle, sum reduced coordinates for frame i
                spos = spos + abpos(:,j).*coff;     % For each particle, sum reduced coordinates for frame i-1
                coff = dimm(j).*coff;               % Product of reduced dimensions ( max volume)
            end

            si = si+1;      % Start indexing at 1 for Matlab
            spos = spos+1;  % Start indexing at 1 for Matlab

            nblocks = coff;
            % trim down (intersect) the hypercube if its too big to fit in the
            % particle volume. (i.e. if dimm(j) lt 3)

            cub = cube;
            deg = find( dimm < 3);
            if ~isempty(deg)
                for j = 0:length(deg)-1
                    cub = cub((cub(:,deg(j+1)) < dimm(deg(j+1))),:);
                end
            end

            % calculate the "s" coordinates of hypercube (with a corner @ the origin)
            scube = zeros(length(cub(:,1)),1);
            coff = 1;
            for j=1:dim
                scube = scube + cub(:,j).*coff;
                coff = coff*dimm(j);
            end

            % shift the hypercube "s" coordinates to be centered around the origin
            coff = 1;
            for j=1:dim
                if dimm(j) > 3
                    scube = scube - coff;
                end
                coff = dimm(j).* coff;
            end
            scube = mod((scube + nblocks),nblocks);

            % si:
            %   - a vector equal to the number of the particles in frame i
            %   - each element of the vector is the location of that
            %     particle in AB space (a number between 1 to nblocks)
            % isort:
            %   - a vector equal to the number of the particles in frame i
            %   - each element convert the sorted si (i.e. strt and fnsh)
            %     back to the unsorted si
            % strt:
            %   - a vector equal to the size of AB space (nblocks)
            %   - index equal to -1 if no particle is in the superpixel
            %   - index equal # of the FIRST particle (sorted list) in that superpixel
            % fnsh:
            %   - a vector equal to the size of AB space (nblocks)
            %   - index equal to 0 if no particle from frame i is in the superpixel
            %   - index equal # of the LAST particle (sorted list) from frame i in that superpixel

            % table
            % si(isort) makes sure they are in increasing order
            [~,isort] = sort(si);
            strt = zeros(nblocks,1) -1;
            fnsh = zeros(nblocks,1);
            for j=1:m
                if strt(si(isort(j))) == -1
                    strt(si(isort(j))) = j;
                    fnsh(si(isort(j))) = j;
                else
                    fnsh(si(isort(j))) = j;
                end
            end

            % loop over the old particles, and find those new particles in the 'cube'.
            %
            % coltot:
            %    - a vector equal to the number of particles in frame i
            %    - each element is the number of particles in frame i-1 that
            %    overlap with particles in frame i
            % rowtot:
            %    - a vector equal to the number of particles in frame i-1
            %    - each element is the number of particles in frame i that
            %    overlap with particles in frame i-1
            % which1:
            %    - a vector equal to the number of particles in frame i-1
            %    - each element is the number of the FIRST particles in frame i that
            %    overlap with particles in frame i-1

            coltot = zeros(m,1);
            rowtot = zeros(n,1);
            which1 = zeros(n,1);

            % n is the number of particle in previous frame (old particles)
            for j = 1:n
                map = fix(-1);
                scub_spos = scube + spos(j);

                % s is cells near or at particle j in frame i-1
                % find particle in frame i that is near particle j in frame
                % i-1 ... these are given by strt(s(w))

                % MATLAB: need to remove 0 indexing and replace with max
                % value
                s = mod(scub_spos,nblocks); s(s==0) = nblocks;

                w = find(strt(s) ~= -1); ngood=length(w);
                if (ngood ~= 0)
                    s = s(w);

                    % map is row number of new particles (frame i) near
                    % particle j (frame i-1)

                    % s is the location in the ab that has a particle near
                    % particle j

                    for k=1:ngood
                        map = [map;isort(strt(s(k)):fnsh(s(k)))];
                    end
                    map = map(2:end);

                    % find those trivial bonds
                    distq = zeros(length(map),1);
                    for d=1:dim
                        distq = distq + (xyi(map,d) - pos(j,d)).^2;
                    end
                    ltmax = distq < maxdisq;
                    rowtot(j) = sum(ltmax);

                    if rowtot(j) >= 1
                        w = find(ltmax == 1);
                        coltot( map(w) ) = coltot( map(w)) +1;
                        which1(j) = map( w(1) );
                    end
                end
            end
            % Record the particles that are a one-to-one match.

            % resx:
            %   - Size = [Time of "working" copy (i.e. zspan), N particle in frame i-1]
            %   - Identifies the particle number that matchs frame to
            %   frame.

            ntrk = fix(n - sum(rowtot == 0));       % Particles in frame i-1 that overlaps with a particle in frame i
            w = find(rowtot == 1); ngood = length(w);
            if ngood ~= 0
                ww = find(coltot( which1(w) ) == 1); ngood = length(ww);
                if ngood ~= 0
                    ndx = w(ww);
                    resx(ispan,ndx) = eyes(which1(ndx));
                    found(which1(ndx)) = 1;
                    rowtot(ndx) = 0;
                    coltot(which1(ndx)) = 0;
                end
            end
            labely = find( rowtot > 0); ngood = length(labely);
            if ngood ~= 0
                labelx = find( coltot > 0);
                nontrivial = 1;
            else
                nontrivial = 0;
            end
            if verbose
                disp(['Frame ' int2str(i),': ',int2str(sum(found)),...
                    ' trival bonds out of ', int2str(m), ' particles.']);
            end
        else
            if verbose
                disp('Use simple N^2 time routine to calculate trivial bonds')
            end
            %   or: Use simple N^2 time routine to calculate trivial bonds

            % let's try a nice, loopless way!
            % don't bother tracking permanently lost guys.
            wh = find( pos(:,1) >= 0); ntrack = length(wh);
            if ntrack == 0
                disp('Warning: No valid particles to track!')
                break
            end

            xmat=mod(reshape(0:1:m*ntrack-1,m,ntrack),m)'+1;
            ymat=mod(reshape(0:1:m*ntrack-1,ntrack,m),ntrack)+1;

            for d = 1:dim
                x2 = xyi(:,d);  % Postion of particles in new frame
                x1 = pos(wh,d); % Postion of particles in previous frame

                % Matlab has an issue with indexing matrices with matrices.
                % This is from subref.m:
                % SUBSREF Subscripted reference.
                %   A(I) is an array formed from the elements of A specified by the
                %   subscript vector I.  The resulting array is the same size as I except
                %   for the special case where A and I are both vectors.  In this case,
                %   A(I) has the same number of elements as I but has the orientation of A.
                % FIX: Adding an if statement for the case where I (xmat or ymat) is a vector

                if min(size(xmat)) == 1 % Special case of xmat is a vector
                    if d == 1
                        dq = (reshape(x2(xmat),size(xmat)) - reshape(x1(ymat),size(ymat))).^2;
                    else
                        dq = dq + (reshape(x2(xmat),size(xmat)) - reshape(x1(ymat),size(ymat))).^2;
                    end
                else % normal case when it is not a vector
                    if d == 1
                        dq = (x2(xmat) - x1(ymat)).^2;
                    else
                        dq = dq + (x2(xmat) - x1(ymat)).^2;
                    end
                end
            end
            % dq is a displacement matrix between all particles in one frame and the
            % next
            ltmax = dq < maxdisq;

            % figure out which trivial bonds go with which
            rowtot = zeros(n,1);
            rowtot(wh) = sum(ltmax,2);   % Number of matching particles in from i+1 for each particle in frame i
            if ntrack > 1
                coltot = sum(ltmax,1);   % Number of matching particles in from i for each particle in frame i+1
            else
                coltot = ltmax;
            end
            which1 = zeros(n,1);
            for j=1:ntrack
                [mx, w] = max(ltmax(j,:));
                which1(wh(j)) = w;      % Array of best guess matches
            end

            ntrk = fix( n - sum(rowtot == 0));
            w= find( rowtot == 1) ;
            ngood = length(w);
            if ngood ~= 0
                ww = find(coltot(which1(w)) == 1);
                ngood = length(ww);
                if ngood ~= 0
                    resx( ispan, w(ww) ) = eyes( which1( w(ww)));
                    found(which1( w(ww))) = 1;
                    rowtot(w(ww)) = 0;
                    coltot(which1(w(ww))) = 0;
                end
            end
            labely = find( rowtot > 0);
            ngood = length(labely);

            if ngood ~= 0
                labelx = find( coltot > 0);
                nontrivial = 1;
            else
                nontrivial = 0;
            end
        end

        %%
        %   THE TRIVIAL BOND CODE ENDS

        % labelx:
        %   - Particles in frame i that still need to be matched
        % labely:
        %   - Particles in frame i-1 that still need to be matched

        if nontrivial
            if verbose
                disp('Performing nontrivial bond calculation')
                disp(['MATCHING: ',int2str(length(labelx)),' particles (frame ',int2str(i), ...
                    ') TO ', ...
                    int2str(length(labely)),' particles (frame ',int2str(i-1),')']);
            end
            xdim = length(labelx);
            ydim = length(labely);
            %  make a list of the non-trivial bonds
            bonds = zeros(1,2);
            bondlen = 0;
            for j=1:ydim
                distq = zeros(xdim,1);
                for d=1:dim
                    distq = distq + (xyi(labelx,d) - pos(labely(j),d)).^2;
                end
                w= find(distq <  maxdisq)'; % KBA removed a -1, not sure why it was there
                ngood = length(w);
                newb = [w;(zeros(1,ngood)+j)];
                bonds = [bonds;newb'];
                bondlen = [ bondlen;distq(w) ];
            end
            bonds = bonds(2:end,:);
            bondlen = bondlen(2:end);
            numbonds = length(bonds(:,1));
            mbonds = bonds;

            % bonds:
            %   - Particles in frames i and i-1 that have a bond (are less than max disp)
            %   - Column 1 is particles in frame i (reduced numbering)
            %   - Column 2 is particles in frame i-1 (reduced numbering)
            % bondlen
            %   - the length of each bond
            % xdim & ydim
            %   - number of unmatched particles in frames i and i-1

            if max([xdim,ydim]) < 4
                nclust = 1;
                maxsz = 0;
                mxsz = xdim;
                mysz = ydim;
                bmap = zeros(length(bonds(:,1)),1) - 1;
            else

                %   THE SUBNETWORK CODE BEGINS
                lista = zeros(numbonds,1);
                listb = zeros(numbonds,1);
                nclust = 0;
                maxsz = 0;
                thru = xdim;
                while thru ~= 0
                    % thru:
                    %  - number of bonds still to be accounted for in frame i
                    %
                    %   the following code extracts connected sub-networks of the non-trivial
                    %   bonds.  NB: lista/b can have redundant entries due to
                    %   multiple-connected subnetworks.

                    % Identify first bond pair (using bonds) not matched.
                    % Record the particle number
                    % lista:
                    %   - particle in frame i part of the cluster
                    % listb:
                    %   - particle in frame i-1 part of the cluster

                    w = find(bonds(:,2) >= 1);
                    lista(1) = bonds(w(1),2);
                    listb(1) = bonds(w(1),1);
                    bonds(w(1),:) = -(nclust+1);    % Identify bond as being part of a cluster
                    adda = 1; addb = 1;
                    donea = 1; doneb = 1;
                    repeat1=true;

                    % Identify clusters
                    % Identify each bonds (particle pair) in a cluster and
                    % label them with the cluster numbers.
                    % When a new particle is identified as being in the
                    % cluster, search for particles it is bonded to and add
                    % them to the cluster.

                    while repeat1
                        if (donea ~= adda+1)
                            w = find(bonds(:,2) == lista(donea));
                            ngood = length(w);
                            if ngood ~= 0
                                listb(addb+1:addb+ngood,1) = bonds(w,1);
                                bonds(w,:) = -(nclust+1);
                                addb = addb+ngood;
                            end
                            donea = donea+1;
                        end

                        if (doneb ~= addb+1)
                            w = find(bonds(:,1) == listb(doneb));
                            ngood = length(w);
                            if ngood ~= 0
                                lista(adda+1:adda+ngood,1) = bonds(w,2);
                                bonds(w,:) = -(nclust+1);
                                adda = adda+ngood;
                            end
                            doneb = doneb+1;
                        end
                        repeat1 = ~((donea == adda+1) && (doneb == addb+1));
                    end

                    % Determine size of cluster
                    xsz = numel(unique(listb(1:doneb-1),'sorted'));
                    ysz = numel(unique(lista(1:donea-1),'sorted'));

                    % Record largest size
                    if xsz*ysz > maxsz
                        maxsz = xsz*ysz;
                        mxsz = xsz;
                        mysz = ysz;
                    end

                    thru = thru -xsz;   % Keep going until all particles in frame i are part of a cluster.
                    nclust = nclust + 1;
                end
                bmap = bonds(:,1);
            end
            % THE SUBNETWORK CODE ENDS
            %
            % mbonds
            %   - Old version of bonds listing the particles in each bond
            %   - Particles in frames i and i-1 that have a bond (are less than max disp)
            %   - Column 1 is particles in frame i (reduced numbering)
            %   - Column 2 is particles in frame i-1 (reduced numbering)
            % bonds
            %   - NOW, lists the cluster each bond is associated with (as a negative number).
            % bmap
            %   - List of cluster number for each bond

            if verbose
                disp(['Frame ',int2str(i),': ','Permuting ', int2str(nclust),' networks']);
                disp(['  Max. network: ',int2str(mxsz),' x ',int2str(mysz)]);
            end

            %   THE PERMUTATION CODE BEGINS

            for nc =1:nclust
                w = find( bmap == -1*(nc)); nbonds = length(w);
                bonds = mbonds(w,:);
                lensq = bondlen(w);
                uold = unique(bonds(:,1),'sorted');
                nold = numel(uold);
                unew = unique(bonds(:,2),'sorted'); % KBA - IDL didn't sort?
                nnew = numel(unew);

                % check that runtime is not excessive
                if nnew  > 5
                    rnsteps = 1;
                    for ii = 1:nnew
                        rnsteps = rnsteps*...
                            numel(find(bonds(:,2) == unew(ii)));
                        if rnsteps  > 5.e4
                            disp(...
                                ' Warning: difficult combinatorics encountered.')
                        end
                        if rnsteps  > 2.e5
                            disp(...
                                ' Excessive Combinatorics! Try reducing maxdisp.')
                        end
                    end
                end

                st = zeros(nnew,1);
                fi = zeros(nnew,1);
                h = zeros(nbonds,1);
                ok = ones(nold,1);
                nlost = (nnew - nold) > 0;

                for ii=1:nold
                    h((bonds(:,1) == uold(ii))) = ii;
                end
                st(1) = 1 ;
                fi(nnew) = nbonds; % check this later

                if nnew > 1
                    sb = bonds(:,2);
                    sbr = circshift(sb,1);
                    sbl = circshift(sb,-1);
                    st(2:end) = find( sb(2:end) ~= sbr(2:end)) + 1;
                    fi(1:nnew-1) = find( sb(1:nbonds-1) ~= sbl(1:nbonds-1));
                end
                checkflag = 0;
                repeat2=true;
                while repeat2

                    pt = st -1;
                    lost = zeros(nnew,1);
                    who = 0;
                    losttot = 0;
                    mndisq = nnew*maxdisq;
                    repeat3=true;
                    while repeat3
                        if pt(who+1) ~= fi(who+1)
                            w = find( ok( h( pt( who+1 )+1:fi( who+1 ) ) ) ); % check this -1
                            ngood = length(w);
                            if ngood > 0
                                if pt(who+1) ~= (st(who+1)-1)
                                    ok(h(pt(who+1))) = 1;
                                end
                                pt(who+1) = pt(who+1) + w(1);
                                ok(h(pt(who+1))) = 0;
                                if who == nnew -1
                                    ww = find( lost == 0);
                                    dsq = sum(lensq(pt(ww))) + losttot*maxdisq;

                                    if dsq < mndisq
                                        minbonds = pt(ww);
                                        mndisq = dsq;
                                    end
                                else
                                    who = who+1;
                                end
                            else
                                if ~lost(who+1) && (losttot ~= nlost)
                                    lost(who+1) = 1;
                                    losttot = losttot + 1;
                                    if pt(who+1) ~= st(who+1) -1
                                        ok(h(pt(who+1))) = 1;
                                    end
                                    if who == nnew-1
                                        ww = find( lost == 0);
                                        dsq = sum(lensq(pt(ww))) + losttot*maxdisq;
                                        if dsq < mndisq
                                            minbonds = pt(ww);
                                            mndisq = dsq;
                                        end
                                    else
                                        who = who + 1;
                                    end

                                else
                                    if pt(who+1) ~= (st(who+1) -1)
                                        ok(h(pt(who+1))) = 1;
                                    end
                                    pt(who+1) = st(who+1) -1;
                                    if lost(who+1)
                                        lost(who+1) = 0;
                                        losttot = losttot -1;
                                    end
                                    who = who -1;
                                end
                            end
                        else
                            if ~lost(who+1) && (losttot ~= nlost)
                                lost(who+1) = 1;
                                losttot = losttot + 1;
                                if pt(who+1) ~= st(who+1)-1
                                    ok(h(pt(who+1))) = 1;
                                end
                                if who == nnew -1
                                    ww = find( lost == 0);
                                    dsq = sum(lensq(pt(ww))) + losttot*maxdisq;

                                    if dsq < mndisq
                                        minbonds = pt(ww);
                                        mndisq = dsq;
                                    end
                                else
                                    who = who + 1;
                                end
                            else
                                if pt(who+1) ~= st(who+1) -1
                                    ok(h(pt(who+1))) = 1;
                                end
                                pt(who+1) = st(who+1) -1;
                                if lost(who+1)
                                    lost(who+1) = 0;
                                    losttot = losttot -1;
                                end
                                who = who -1;
                            end
                        end
                        repeat3=~(who == -1);
                    end

                    checkflag = checkflag + 1;
                    if checkflag == 1
                        %   we need to check that our constraint on nlost is not forcing us away from
                        %   the minimum id's
                        plost = min([fix(mndisq/maxdisq) , (nnew -1)]);
                        if plost > nlost
                            nlost = plost;
                        else
                            checkflag = 2;
                        end
                    end
                    repeat2 = ~(checkflag == 2);
                end

                %   update resx using the minimum bond configuration
                resx(ispan,labely(bonds(minbonds,2))) = eyes(labelx(bonds(minbonds,1)));
                found(labelx(bonds(minbonds,1))) = 1;
            end
            %   THE PERMUTATION CODE ENDS
        else
            if verbose
                disp([int2str(i),': Only trivial networks'])
            end
        end
        %     here we want to update our initial position estimates
        w = find(resx(ispan,:) >= 0);
        nww = length(w);

        if nww > 0
            pos(w,:) = xyzs( resx(ispan,w) , 1:dim);
            if  ~isempty(goodenough)
                nvalid(w) = nvalid(w) + 1;
            end
        else
            disp(' Warning, tracking zero particles!')
        end

        %     we need to add new guys, as appropriate.
        newguys = find(found == 0); nnew = length(newguys);
        if (nnew  > 0) && isempty(inipos)
            newarr = zeros(zspan,nnew) -1;
            resx = [resx,newarr];
            resx(ispan,n+1:end) = eyes(newguys);
            pos = [[pos];[xyzs(eyes(newguys),1:dim)]];
            nmem = zeros(nnew,1);
            mem = [mem;nmem];
            nun = 1:nnew;
            uniqid = [uniqid,((nun) + maxid)];
            maxid = maxid + nnew;
            if  ~isempty(goodenough)
                dumphash = [dumphash;zeros(1,nnew)'];
                nvalid = [nvalid;zeros(1,nnew)'+1];
            end
            n = n + nnew;
        end
    else
        disp([' Warning- No positions found for t=',int2str(i),"!"])
    end

    %   update the 'memory' array
    w = find( resx(ispan,:) ~= -1); nok = length(w);
    if nok ~= 0, mem(w) =0; end
    mem = mem + (resx(ispan,:)' == -1);

    %  if a guy has been lost for more than memory times, mark him as permanently
    %  lost.  For now, set these guys to pos = ( -maxdisp, -maxdisp, ... ),
    %  so we can never track them again. It would be better to make a smaller
    %  pos, but then we'd have to change 'n', which would be gnarly.
    wlost = find(mem == memory_b+1); nlost =length(wlost);
    if nlost > 0
        pos(wlost,:) = -maxdisp;
        % check to see if we should 'dump' newly lost guys
        if  ~isempty(goodenough)
            wdump = find(nvalid(wlost) < goodenough); ndump = length(wdump);
            if ndump > 0
                dumphash(wlost(wdump)) = 1;
            end
        end
    end

    %  we need to insert the working copy of resx into the big copy bigresx
    %  do our house keeping every zspan time steps (dumping bad lost guys)

    if (ispan == zspan) || (i == z)

        %  if a permanently lost guy has fewer than goodenough valid positions
        %  then we 'dump' him out of the data structure- this largely alleviates
        %  memory problems associated with the 'add' keyword and 'noise' particles
        %  To improve speed- do it infrequently.
        % in case we've added some we need to pad out bigresx too
        nold = length(bigresx(1,:));
        nnew = n-nold;
        if nnew > 0
            newarr = zeros(z,nnew) -1;
            bigresx = [bigresx,newarr];
        end

        if  ~isempty(goodenough)
            if (sum(dumphash(:)) > 0)
                if verbose, disp('Dumping bad trajectories...'); end
                wkeep = find(dumphash == 0);
                nkeep = length(wkeep);
                resx = resx(:,wkeep);
                bigresx = bigresx(:,wkeep); % this really hurts runtime
                pos = pos(wkeep,:);
                mem = mem(wkeep);
                uniqid = uniqid(wkeep);
                nvalid = nvalid(wkeep);
                n = nkeep;
                dumphash = zeros(nkeep,1);
            end
        end

        if ~verbose && isempty(quiet)
            if remove == 1
                reverseStr = repmat(sprintf('\b'), 1, N_rm);
            else
                reverseStr = [];
                remove = 1;
            end
            disp_string = ['Frame ', int2str(i),' of ', int2str(z), ' done. ', ...
                'Current frame has ',int2str(ntrk),' particles. '];
            disp([reverseStr,disp_string]);
            N_rm = length(disp_string)+1;
        end
        bigresx(i-(ispan)+1:i,:) = resx(1:ispan,:);
        resx = zeros(zspan,n) - 1;

        %  We should pull permanently lost guys, parse them and concat them
        %  onto the 'output list', along with their 'unique id' number to
        %  make scanning the data files a little easier.  Do infrequently.
        wpull = find(pos(:,1) == -maxdisp); npull = length(wpull);

        if npull > 0
            lillist = zeros(1,2);
            for ipull=1:npull
                wpull2 = find(bigresx(:,wpull(ipull)) ~= -1);
                npull2 = length(wpull2);
                thing = [bigresx(wpull2,wpull(ipull)),zeros(npull2,1)+uniqid(wpull(ipull))];
                lillist = [lillist;thing];

            end
            olist = [[olist];[lillist(2:end,:)]];
        end

        %     now get rid of the guys we don't need any more....
        %     but watch out for when we have no valid particles to track!
        wkeep = find(pos(:,1) >= 0); nkeep = length(wkeep);
        if nkeep == 0
            disp('Were going to crash now, no particles....')
        end

        resx = resx(:,wkeep);
        bigresx = bigresx(:,wkeep);

        pos = pos(wkeep,:);
        mem = mem(wkeep);
        uniqid = uniqid(wkeep);
        n = nkeep;
        dumphash = zeros(nkeep,1);
        if  ~isempty(goodenough)
            nvalid = nvalid(wkeep);
        end

    end  % the big loop over z time steps....
end

% Main output from "big loop over z time steps"
%
% olist
%   - stands for output list
%   - list of all particles that should be part of the output result
%   - second column is the particle ID; is using goodenough, IDs will jump
%   around since some tracks are removed.
%   - first column is the is the row number in xyzz (i.e. pt) for the particle.
%
% For final sweep, also need ...
% bigresx
%   - Each row is a frame in the video
%   - Each column is a particle that has been linked over multiple frames.
%   The number is the row number of that particle in xyzs.
%   - '-1' means is wasn't found in that frame.
%
% uniqid
%   - Unique ID numbers of trajectories in bigresx
%   - Arbitrarily set
%
% pos
%   - xy positions of the final frame analyzed
%   - used for truncated bigresx
%
% n
%   - the number of particles being tracked in the last frame

%  make final scan of bigresx for trajectories that are long enough
%  only save those trajectories

if  ~isempty(goodenough)
    nvalid = sum(bigresx >= 0 ,1);
    wkeep = find(nvalid >= goodenough); nkeep = length(wkeep);

    if nkeep < n    % KBA: Not sure why this is needed.
        bigresx = bigresx(:,wkeep);
        n = nkeep;
        uniqid = uniqid(wkeep);
        pos = pos(wkeep,:);
    end
end

%  make the final scan to 'pull' everybody else into the olist.
wpull = find( pos(:,1) ~= -2*maxdisp); npull = length(wpull);

if npull > 0
    lillist = zeros(1,2);
    for ipull=1:npull
        wpull2 = find(bigresx(:,wpull(ipull)) ~= -1);   % Find elements in column that aren't -1
        npull2 = length(wpull2);
        thing = [bigresx(wpull2,wpull(ipull)),zeros(npull2,1)+uniqid(wpull(ipull))];
        lillist = [lillist;thing];
    end
    olist = [olist;lillist(2:end,:)];
end
olist = olist(2:end,:);

%  free up a little memory for the final step!
bigresx = [];
resx = [];

% need to make up a result array!
if verbose
    disp('Preparing result array...')
end
nolist = length(olist(:,1));
res = zeros(nolist,dd+1);
for j=1:dd
    res(:,j) = xyzs(olist(:,1),j);
end
res(:,dd+1) = olist(:,2);

% Remove dummy particles if inserted for missing frames 
if ~isempty(dummy)
    w=find(res(:,1)==0);
    res(w,:) = [];
end

% Renumber IDs incase there are any gaps
ndat=numel(res(1,:));
[~,u,~] = unique(res(:,ndat),'stable');
ntracks=numel(u);
u=[u;length(res(:,ndat))+1];
for i=1:ntracks
    res(u(i):u(i+1)-1,ndat) = i;
end

if isempty(quiet)
    disp(['****** TRACKING STATS ******']);
    disp(['Total number of tracks found: ', int2str(ntracks)]);
    trkLength=u(2:end)-u(1:end-1);
    [maxLength,ID]=max(trkLength);
    disp(['Maximum track length: ', int2str(maxLength), ' (ID: ',int2str(ID),')']);
    disp(['Mean track length: ', int2str(mean(trkLength))]);
    disp(['Median track length: ', int2str(median(trkLength))]);
end

%end
% ***************** end of track.pro