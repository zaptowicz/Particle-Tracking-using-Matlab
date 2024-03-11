function rho = artifact(img, bphi, extent, varargin)
% FUNCTION NAME:
%   artifact
%
% DESCRIPTION:
%   Determines optical artifact when particles approach each other. Assumed
%   images are particles approaching each other can be superimposed.
%
% INPUT (REQUIRED)
%             img: Image of dilute particles to be analyzed
%            bphi: Approximate feature size for bpass
%          extent: Size of centroid in findfeatures
%
% INPUT (OPTIONAL)
%         masscut: Determines minimum brightness a blob should be to be considered
%                  a particle.
%      separation: Minimum separation of a particle from other particles to
%                  be analyzed.
%    plot_results: ['y'] Show plot of results. 
%
% OUTPUT:
%         rho(:1): Actual separation of particle images
%         rho(:2): Apparent separation of particle images
%
% CALLING SEQUENCE:
%   rho = artifact(img, 11, 13, masscut = 3000, separation = 30, plot_results='y');
%
% NOTES:
%   IDL Version - A version of this code was originally written by Yilong
%   Han (Hong Kong University of Science and Technology).
%
% REVISION HISTORY:
%   03/11/2024 - K Aptowicz (WCU)
%       * Translated to MATLAB

%% Reading and setting parameters

% Set default values for optional parameters
default_masscut = 3000; % Default value
default_separation = []; % Default value
default_plot_results = []; % Default value

% Create fields for all optionals inputs
p = inputParser;

% Keywords
addParameter(p,'masscut',default_masscut,@isnumeric)
addParameter(p,'separation',default_separation,@isnumeric);
addOptional(p,'plot_results', default_plot_results)

% populate optional parameters from inputs
parse(p,varargin{:});
masscut = p.Results.masscut;
separation = p.Results.separation;
plot_results = p.Results.plot_results;

%%

% If 'separation' not set, make it 5 times bphi
if isempty(separation)
    minsep = 5*bphi;
else
    minsep=separation;
end

a=img;
height = size(a,1);
width = size(a,2);

% Find candidate particles
b = bpass(a,1,bphi);
f = findfeatures(b,extent,masscut=masscut, iterate = 'y',quiet = 'y');
np = size(f,1);

% Find the shortest nearest-neighbor distance for each particle
X = f(:,1)*ones(1,np);
XX = (X-X').^2;
Y = f(:,2)*ones(1,np);
YY = (Y-Y').^2;
dis = sqrt(XX+YY);
dis(dis==0) = inf;
dis=min(dis);

% Select particles that are distant from other particles
w=find(dis>minsep);
if isempty(w)
    disp("No isolated particles!")
else
    f=f(w,:);
end

% Set minimum distance from edge of image. minsep is too conversative. Make
% slightly smaller. 
minsep = minsep/sqrt(2);
minsep = round(minsep);

% Region of interest around particle
xmin = floor(f(:,1) - minsep);
xmax = floor(xmin + 2.* minsep);
ymin = floor(f(:,2) - minsep);
ymax = floor(ymin + 2.* minsep);

% Pick particles far from the edges
w =  find(xmin > 1 & xmax < width & ymin > 1 & ymax < height);
if isempty(w)
    disp("No isolated particles far from edges!")
else
    f=f(w,:);
    np = length(w);
    disp(['Analyzing images of ', num2str(np), ' particles'])
end

if ~isempty(plot_results)
    figure
    fo=fover2d(a,f,radius=bphi,circle='y');
    axis off
    title('Isolated Particles')
end

xmin = xmin(w);
xmax = xmax(w);
ymin = ymin(w);
ymax = ymax(w);

% calculate apparent separations
aa = a(ymin(1):ymax(1),xmin(1):xmax(1));
rho = artifact_rho_avg(aa, bphi, masscut);
for i = 2:np
    aa = a(ymin(i):ymax(i),xmin(i):xmax(i));
    rho = rho + artifact_rho_avg(aa, bphi, masscut);
end
rho(:,2) = rho(:,2)./rho(:,1);
rho(:,1) = [1:1:size(rho,1)];

if ~isempty(plot_results)
    figure
    plot(rho(:,1), (rho(:,2)-rho(:,1)),'ko','MarkerFaceColor', 'k')
    xlabel('spacing of particles (pixels)') %\it to make x italics
    ylabel('measurement discrepancy (pixels)')
    title('Discrepancy in Measured Spacing of Particles')
    grid on
end

end

% ************************************************************************
% *** artifact_rho_avg *** calculate the average rho for an image of a particle
% and its tranpose.  Thus averages both x and y image translations.

function rho = artifact_rho_avg(a, bphi, masscut)
rho = artifact_rho(a, bphi, masscut);
rho = rho + artifact_rho(a', bphi, masscut);
end

% ************************************************************************
% *** artifact_rho *** calculates the distance between two images of a
% particle using bpass and findfeatures as the images are brough closer and
% closer together.
% rho(:,1) = 1 means particles were identified at both images, otherwise 0
% rho(:,2) = absolute distance between found particles
function rho = artifact_rho(a, bphi, masscut)
a=double(a);
w = size(a,1);
h = size(a,2);
a = a - mean(a(:));

rho = zeros(w,2);

for r = 1:w
    e = zeros(h,2.*w);
    e(:,1:w) = a;
    e(:,r+1:r+w) = e(:,r+1:r+w) + a;

    eb = bpass(e, 1, bphi);
    ep = findfeatures(eb, bphi+2, masscut = masscut, iterate = 'y', quiet = 'y');
    if size(ep,1) == 2
        rho(r, 1) = 1;
        rho(r, 2) = abs(ep(1,1) - ep(2,1));
    end
end
end