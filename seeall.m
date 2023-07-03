function seeall(tarray,varargin)

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
for i=1:ntracks
    tarray(u(i):u(i+1)-1,ndat) = i;
end

% Plot each track
ndat=numel(tarray(1,:));
[~,u,~] = unique(tarray(:,ndat));
ntracks=numel(u);
u=[u;length(tarray(:,ndat))+1];
set(gcf,'position',[10,10,1200,600])
for i=1:ntracks
    t1=tarray((tarray(:,7)==i),:);
    subplot(2,3,1)
    plot(t1(:,1),t1(:,2),'-o')
    axis equal
    title('track')
    xlabel('x')
    ylabel('y')
    subplot(2,3,2)
    plot(t1(:,1))
    xlabel('time')
    ylabel('x')
    subplot(2,3,3)
    plot(t1(:,2))
    xlabel('time')
    ylabel('y')
    subplot(2,3,4)
    plot(t1(:,3))
    xlabel('time')
    ylabel('brightness')
    subplot(2,3,5)
    plot(t1(:,4))
    xlabel('time')
    ylabel('radius')
    subplot(2,3,6)
    plot(t1(:,5))
    xlabel('time')
    ylabel('eccentricity')
    pause
end
close(gcf) 
end