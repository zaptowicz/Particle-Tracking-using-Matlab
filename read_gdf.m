function res = read_gdf(filespec)

% COMMENT: This code was translated from read_gdf.pro code provided on
%   David Grier' github page:
%   https://github.com/davidgrier/utility/blob/master/read_gdf.pro
%   It has NOT been thoroughly tested. And was not written to handle different OS writing the file.
%   As bugs are dealt with, updates will be posted.
%
% MODIFICATION HISTORY:
%   Translated by Kevin Aptowicz, WCU 10/'22
%   EDITS:
%       03/15/2024 - KBA (WCU)
%           - Fixed issue with reading different types of data formats. Now
%           can handle uint8, int16, int32, single, and double
%
% ;+
% ; NAME:
% ;		read_gdf
% ; PURPOSE:
% ;		Read in data files created by WRITE_GDF.
% ;
% ; CATEGORY:
% ;		General Purpose Utility
% ; CALLING SEQUENCE:
% ;		data = read_gdf(file)
% ; INPUTS:
% ;		file:	Complete pathname of the file to be read.
% ; OUTPUTS:
% ;		data:	Data structure.  For example, if the original
% ;			data was stored as an array of bytes, then
% ;			DATA will be returned as an array of bytes also.
% ; RESTRICTIONS:
% ;		Current implementation does not support structures or
% ;		arrays of structures.
% ; PROCEDURE:
% ;		Reasonably straightforward.
% ;		Determines if the file is ASCII or binary, reads the size
% ;		and dimension info from the header, and reads in the data
% ; MODIFICATION HISTORY:
% ; Written by David G. Grier, AT&T Bell Laboratories, 9/91
% ; 12/01/1995 DGG: Figures out how to deal with data from different-endian
% ;   machines.
% ; 03/17/2010 DGG: Use file_search to find files.  Code and documentation
% ;   cleanups.
% ;
% ; Copyright (c) 1991-2010 David G. Grier
% ;
% ; UPDATES:
% ;    The most recent version of this program may be obtained from
% ;    http://physics.nyu.edu/grierlab/software.html
% ;
% ; LICENSE:
% ;    This program is free software; you can redistribute it and/or
% ;    modify it under the terms of the GNU General Public License as
% ;    published by the Free Software Foundation; either version 2 of the
% ;    License, or (at your option) any later version.
% ;
% ;    This program is distributed in the hope that it will be useful,
% ;    but WITHOUT ANY WARRANTY; without even the implied warranty of
% ;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% ;    General Public License for more details.
% ;
% ;    You should have received a copy of the GNU General Public License
% ;    along with this program; if not, write to the Free Software
% ;    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
% ;    02111-1307 USA
% ;
% ;    If the Internet and WWW are still functional when you are using
% ;    this, you should be able to access the GPL here:
% ;    http://www.gnu.org/copyleft/gpl.html
% ;-


file = dir(filespec);

if size(file,1) == 0
    disp(['No files matched specification ',filespec])
else
    fid = fopen(file(1).name);

    % Check if ASCII
    gdf=fscanf(fid,'%s',1);
    ascii = (gdf=="GDF");

    if ascii
        fscanf(fid,'%s',2); % Not needed
        sz=fscanf(fid,'%i',1); % number of dimensions
        dim = fscanf(fid,'%i',sz); % array of dimension
        check = fscanf(fid,'%i',2); % Check if Matlab array
        res = fscanf(fid,'%f',Inf); % single vector of data
        if sz > 1
            res = reshape(res,dim');
        end
    else
        frewind(fid)
        mgc = fread(fid,1,'long'); % Grier MAGIC
        sz = fread(fid,1,'long'); % number of dimensions
        dim = fread(fid,sz,'long'); % array of dimension
        type = fread(fid,1,'long'); % format
        N = fread(fid,1,'long'); % number of elements
        if type == 1 % IDL 'byte' -  8-bit unsigned interger
            res = fread(fid,N,'uint8'); % single vector of data
        elseif type == 2 % IDL 'int'  -  16-bit signed interger
            res = fread(fid,N,'int16'); % single vector of data
        elseif type == 3 % IDL 'long'  -  32-bit signed interger
            res = fread(fid,N,'int32'); % single vector of data
        elseif type == 4 % IDL 'float'  -  32-bit single-precision
            res = fread(fid,N,'single'); % single vector of data
        elseif type == 5 % IDL 'double' -  64-bit single-precision
            res = fread(fid,N,'double'); % single vector of data
        else
            disp('READ_GDF: Data type from IDL not programmed')
        end
        if sz > 1
            res = reshape(res,dim');
        end
    end
    if sz == 2 % 2D array
        res=res'; % flip from idl array(col,row) to matlab array(row,col)
    end
    fclose(fid);
end