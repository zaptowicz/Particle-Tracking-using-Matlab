function write_gdf(data, file, varargin)
%
% COMMENT: This code was translated from write_gdf.pro code provided on
%   David Grier' github page:
%   https://github.com/davidgrier/utility/blob/master/write_gdf.pro
%   It has NOT been thoroughly tested. As bugs are dealt with, updates will be
%   posted.
%
% IDL FORMAT: 
%   Header consists of the following lines: 
%       long MAGIC number (082991) (ise to check system format)
%       long array of IDL function size of data (number of dimensions, length of each dimension, number type, number of elements)
%       data using data type of data (eg. byte, long, float) 
% 		
%       If the file is ASCII, the header tag is an identifying
% 		statement, followed by all of the same information in
% 	    ASCII format.
%
% MODIFICATION HISTORY:
%   Translated by Kevin Aptowicz, WCU 10/'22
%       EDITS:
%       03/17/2024 - KBA (WCU)
%           - Edits made so that format more closely matches IDL
%           - Can write IDL formats: BYTE, INT, LONG, FLOAT, and DOUBLE
%
% ;+
% ; NAME:
% ;		write_gdf
% ; PURPOSE:
% ;		Writes IDL-style data to disk in a format which can be
% ;		read back in easily.
% ;
% ; CATEGORY:
% ;		General Purpose Utility
% ; CALLING SEQUENCE:
% ;		write_gdf,data,file
% ; INPUTS:
% ;		data:	Data structure to be written to disk.
% ;		file:	Complete pathname of the file to be written.
% ; KEYWORD:
% ;		ascii:	Produce an ASCII file rather than the default
% ;			binary file.
% ; SIDE EFFECTS:
% ;		Creates a file.
% ; RESTRICTIONS:
% ;		Current version does not support structures or arrays of
% ;		structures.
% ; PROCEDURE:
% ;		Writes a header consisting of a long MAGIC number followed
% ;		by the long array produced by SIZE(DATA), followed by the
% ;		data itself.
% ;		If the file is ASCII, the header tag is an identifying
% ;		statement, followed by all of the same information in
% ;		ASCII format.
% ; MODIFICATION HISTORY:
% ; Written by David G. Grier, AT&T Bell Laboratories, 09/01/1991
% ; 03/17/2010 DGG: Code and documentation cleanups.
% ;
% ; Copyright (c) 1991-2010 David G. Grier
% ;-

%% Reading and setting parameters
% Set default values for optional parameters
default_ascii = [];

% Create fields for all optionals variables
p = inputParser;
addOptional(p,'ascii', default_ascii);

% populate optional parameters from inputs
parse(p,varargin{:});
ascii = p.Results.ascii;

%% Prepare header
% Fix issues with Data
if numel(size(data)) == 2
    data=data'; % to match IDL preferred format of array(coln, row)
end

% Prepare header
MAGIC = int32(082991);
HEADER = "GDF v. 1.0";
fileID = fopen(file,'w');
sz = size(data);
type = class(data);

% Determine number if IDL data type
switch type
    case 'uint8'
        i=1; 
    case 'int16'
        i=2;
    case 'int32'
        i=3;
    case 'single'
        i=4;
    case 'double'
        i=5;
    otherwise
        disp('WRITE_GDF: Data type not programmed for write_gdf')
end

IDL_size = [length(sz), sz, i, numel(data)]; % Use to identify Matlab created array
IDL_size = int32(IDL_size); % Make sure it is in LONG format for IDL

%% Write to file
if ~isempty(ascii) % If ASCII
    fprintf(fileID, '%s\n',HEADER);
    fprintf(fileID,'%i\n',IDL_size(1));
    fprintf(fileID,'%i ',IDL_size(2:end));
    fprintf(fileID,'\n');
    fprintf(fileID, '%f\n', data);
else
    fwrite(fileID, MAGIC,'long');
    fwrite(fileID, IDL_size,'long');
    fwrite(fileID, data, type);
end
fclose(fileID);
end