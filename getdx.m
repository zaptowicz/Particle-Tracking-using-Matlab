function result=getdx(tr,thedt,varargin)


% getdx                   5-22-00   Eric R. Weeks
% patched 8-14-01 to return -1 for nonvalid results
% major bug fixed 6-3-05 for samples with unevenly gridded time stamps
% gdxtrinterp added internally 6-13-05
% Gianguido Cianci: 24 Jan 07. Now rounds timestamp to remedy old
% interpolation bug.

% function getdx,tr,dt,dim=dim
%
% see http://www.physics.emory.edu/~weeks/idl/getdx.html
% for more details.
% =================================================================

% gdxtrinterp.pro -- from "trinterp.pro"
%
% started 6-17-98 by ERW
%
% interpolates gaps in tracked data

%% Reading and setting parameters
% Set default values for optional parameters
default_dim = 2;
default_flag = [];

% Create fields for all optionals inputs
p = inputParser;

% Variables
addParameter(p,'dim',default_dim,@isnumeric)

% Keywords
addOptional(p,'flag', default_flag)

% populate optional parameters from inputs
parse(p,varargin{:})
dim = p.Results.dim;
flag = p.Results.flag;

trin=gdxtrinterp(tr,flag);

nel=numel(trin(1,:));
trins=circshift(trin,-round(thedt(1),1));
fl1 = trin(:,nel-2);
fl2 = trins(:,nel-2);
result=trins-trin;
% ww=where((result(nel-1,*) eq 0) and (fl1 lt 0.5) and (fl2 lt 0.5))
w=find((result(:,nel) ~= 0) | (round(result(:,nel-1)) ~= round(thedt(1))) | ...
       (fl1 > 0.5) | (fl2 > 0.5));
nw = numel(w);
result=trins(:,1:dim+1)-trin(:,1:dim+1);
trins=0; % free up memory
result(:,dim+1)=sum(result(:,1:dim).^2,2);
result(:,dim+1)=sqrt(result(:,dim+1));
if (nw > 0) 
    result(w,:)=-1;
end

w2=find(trin(:,nel-2) < 0.5);
result=result(w2,:);
w3=find(result(:,dim+1) < 0); nw3=numel(w3);

if (nw3 > 0) 
    result(w3,1:dim) = 0;
end
end
%%
function result = gdxtrinterp(tarray,flag)
% tarray is an array of tracked data
% returns a new tracked data array
%
% flag=1 to append a column of 1's for interpolated values, 0's for
%       original values

s=size(tarray);
result=tarray;
if ~isempty(flag) 
    result=[result,zeros(s(1),1)];
end

ndat=length(tarray(1,:)); 

%old interpolation algorithm could cause non integer value in timestamp
%column. This should fix it. ERW's suggestion
tarray(:,ndat-1) = round(tarray(:,ndat-1));

dp=tarray-circshift(tarray,1,1);
% changed at REC's suggestion by ERW: 1-10-03
% w=where((dp(ndat-1,*) eq 0) and (dp(ndat-2,*) ne 1),ngood)
w=find((dp(:,ndat) == 0) & (dp(:,ndat-1) > 1));
ngood = numel(w);
count = 1;
storeres=zeros(1,ndat);
if (ngood >= 1) 
    totdt=sum(dp(w,ndat-1))-ngood;
    clear dp;% free up memory
    storeres=zeros(totdt,ndat);
    for i=1:ngood 
       x0=tarray(w(i)-1,:);
       x1=tarray(w(i),:);
       dt=x1(ndat-1)-x0(ndat-1);
       t=1.0-(1:1:dt-1)/dt;
       xx = x0'*t + x1'*(1.0-t);
       storeres(count:count+dt-2,:)=xx';
       count = count + dt - 1;
    end
    if ~isempty(flag) % append a column of 1's for interpolated values
       storeres(:,ndat+1)=1;
    end
    result=vertcat(result,storeres);
end

% put in Victor Breedveld's patch here, to make sort work in DOS
[~,ind]=sort(result(:,ndat));
result=result(ind,:);
[~,ind]=unique(result(:,ndat));
u=vertcat(ind,numel(result(:,1))+1);
nu=numel(u);
for i=1:nu-1 
    tpart=result(u(i):u(i+1)-1,:);
    [~,ind]=sort(tpart(:,ndat-1));
    result(u(i):u(i+1)-1,:) = tpart(ind,:);
end
% end of Victor's patch

result(:,ndat-1:ndat)=round(result(:,ndat-1:ndat));

if ~isempty(flag) 
    result(:,ndat-1:end)=result(:,[ndat+1,ndat-1,ndat]);
end
end
% function gdxtrinterp ends