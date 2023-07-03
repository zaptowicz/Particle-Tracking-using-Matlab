function write_gdf(data, file, varargin)
%
% COMMENT: This code was translated from write_gdf.pro code provided on
%   David Grier' github page:
%   https://github.com/davidgrier/utility/blob/master/write_gdf.pro
%   It has NOT been thoroughly tested. As bugs are dealt with, updates will be
%   posted.
%
% MODIFICATION HISTORY:
%   Translated by Kevin Aptowicz, WCU 10/'22
%   Adapted by XXXX: XX/'XX; XXXXXXX
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

MAGIC = int32(082991);
HEADER = "GDF v. 1.0";
fileID = fopen(file,'w');
sz = size(data);

% Tweak to match IDL format of array(coln,row) rather than array(row,coln)
if length(sz) == 2 % 2D array
    data=data'; % to match IDL preferred format of array(coln, row)
    sz2 = [fliplr(sz),[-1,0]]; % Use to identify Matlab created array
else
    sz2 = [sz,[-1,0]]; % Use to identify Matlab created array
end

if ~isempty(ascii)
    fprintf(fileID, '%s\n',HEADER);
    fprintf(fileID,'%i\n',length(sz));
    fprintf(fileID,'%i ',sz2);
    fprintf(fileID,'\n');
    fprintf(fileID, '%f\n', data);
else
    fwrite(fileID, MAGIC,'long');
    fwrite(fileID, length(sz),'long');
    fwrite(fileID, sz2,'long');
    fwrite(fileID, data,'single');
end
fclose(fileID);
end