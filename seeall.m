function seeall(tarray,varargin)
% FUNCTION NAME:
%   seeall
%
% DESCRIPTION:
%   Makes a multi-panel figure exploring the properties of each track. 
%
% INPUT (REQUIRED)
%          tarray: Track array.
%
% INPUT (OPTIONAL)
%              ID: (interger) set equal to the particle ID number to plot.
%                  Plots just this particle track
%      goodenough: (interger) Set the minimum length of a track to plot.
%                   Shorter tracks are ignored. 
%
% OUTPUT: None
%
% CALLING SEQUENCE:
%   >> seeall(t)
%   >> seeall(t,ID=10) 
%
% NOTES:
%   IDL VERSION
%           This code was translated from getdx.pro code provided
%           on Eric Weeks' website:
%           https://physics.emory.edu/faculty/weeks/idl/kit/seeall.pro
% 
%
% MODIFICATION HISTORY:
%   10/29/1998 - Eric Weeks
%       * Wrote the first version
%   07/05/2023 - Kevin Aptowicz (WCU)
%       * Translated latest version to MATLAB

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
    tarray=tarray(tarray(:,end)==ID,:)
end

% Save only tracks of length equal to or greater than goodenough
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

% Renumber IDs to be sequential ...
ndat=numel(tarray(1,:));
[~,u,~] = unique(tarray(:,ndat),'stable');
ntracks=numel(u);
u=[u;length(tarray(:,ndat))+1];
l=u(2:end)-u(1:end-1);
[~,ls] = sort(l,'descend');
IDs=tarray(u(ls),ndat);

% Plot each track
ndat=numel(tarray(1,:));
[~,u,~] = unique(tarray(:,ndat));
ntracks=numel(u);
figure
set(gcf,'position',[10,10,1200,600])
for i=1:ntracks
    t1=tarray((tarray(:,7)==IDs(i)),:);
    subplot(2,3,1)
    plot(t1(:,1),t1(:,2),'-o')
    axis equal
    title(['ID: ',int2str(IDs(i))])
    xlabel('x')
    ylabel('y')
    subplot(2,3,2)
    plot(t1(:,1))
    xlabel('time')
    title('x position')
    subplot(2,3,3)
    plot(t1(:,2))
    xlabel('time')
    title('y position')
    subplot(2,3,4)
    plot(t1(:,3))
    xlabel('time')
    title('brightness')
    subplot(2,3,5)
    plot(t1(:,4))
    xlabel('time')
    title('radius')
    subplot(2,3,6)
    plot(t1(:,5))
    xlabel('time')
    title('eccentricity')
    drawnow
    prompt = "What is the original value? ";
    disp("Next track? Enter to continue. Q to quit.");
    x = input(">> ","s");
    if x == 'Q' | x == 'q'
        return
    end
end
figure;
close(gcf)