function plottr(tarray,varargin)
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
%              ID: (interger) set equal to the particle ID number to plot.
%                  Plots just this particle track
%      goodenough: (interger) Set the minimum length of a track to plot.
%                   Shorter tracks are ignored. 
%
% OUTPUT: NONE
%
% CALLING SEQUENCE:
%   plottr(t)
%   plottr(t, ID=4)
%   plottr(t, goodenough=10)
%
% NOTES :
%   IDL Version - This code was translated from plottr.pro code provided on
%           Eric Weeks' website:
%           https://physics.emory.edu/faculty/weeks/idl/kit/plottr.pro
%   Matlab Version - Added the ability to pick ID, but didn't inlcude a
%           bunch of other stuff from the IDL version
%
% REVISION HISTORY:
%   05/15/1998 - Eric Weeks
%       * Started the code
%   07/05/2023 - K Aptowicz
%       * Translated to MATLAB
%

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