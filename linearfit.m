function res=linearfit(x,y,varargin)
% FUNCTION NAME:
%   linearfit
%
% DESCRIPTION:
%   Fits (x,y) data to a line using least-squares fitting.  See Chapter 8 
%   of the John Taylor book on Error Analysis.   
%
% INPUT (REQUIRED)
%          x: x-values array
%          y: y values array
%
% INPUT (OPTIONAL)
%      quiet: ['y'] Nothing is displayed.  
%
% OUTPUT: 
%        res: row 1 --> [Slope, Uncertainty]
%             row 2 --> [y-intercept, Uncertainty]
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
%   PROCEDURE
%           Expression copied from the text Introduction to Error Analysis
%           by John R. Taylor.
%
% MODIFICATION HISTORY:
%   07/05/2023 - Kevin Aptowicz (WCU)
%       * Wrote program using other matlab scripts. 

%% Reading and setting parameters
% Set default values for optional parameters
default_quiet = [];

% Create fields for all optionals inputs
p = inputParser;

% Variables
% addParameter(p,'dim',default_dim,@isnumeric)

% Keywords
addParameter(p,'quiet',default_quiet)

% populate optional parameters from inputs
parse(p,varargin{:})
quiet = p.Results.quiet;

if length(x) ~= length(y)
    disp('Abort: Vectors not the same length')
    return
end

N=length(x); % Number of datapoints
Delta = N*(sum(x.^2))-(sum(x).^2);  % Determine Delta
A = ((sum(x.^2))*(sum(y))-(sum(x))*sum(x.*y))/Delta; % Calculate Constant A
B = (N*sum(x.*y)-(sum(x))*(sum(y)))/Delta; % Calculate Constant B

sigma_y = sqrt((1/(N-2))*sum((y-A-B*x).^2)); % Calculate sigma y
sigma_A = sigma_y*sqrt(sum(x.^2)/Delta); % Calculate sigma A
sigma_B = sigma_y*sqrt(N/Delta); % Calculate sigma B

res = [[B,sigma_B];[A,sigma_A]];
if isempty(quiet)
    plot(x,y,'ko','MarkerFaceColor', 'k')
    drawnow
    hold on
    plot(x,A+B*x,'-k')
    hold off
    drawnow
    xlabel('{\it x} value') %\it to make x italics
    ylabel('{\it y} value')
    title('Linear fit')
    disp('*** Linear fit parameters ***')
    fprintf('\n')
    disp('Slope | uncertainty') 
    disp([B,sigma_B])
    disp('Intercept | uncertainty') 
    disp([A,sigma_A])
end
end