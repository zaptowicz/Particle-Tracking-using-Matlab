function vid = read_vif(filename, pix_w, pix_h, varargin)
% FUNCTION NAME:
%   read_vif
%
% DESCRIPTION:
%   Reads video file into an array
%trash
% INPUT (REQUIRED):
%        filename: (string) filename of video to be imported
%           pix_w: (double) pixel width of video frames
%           pix_h: (double) pixel height of video frames
%
% INPUT (OPTIONAL):
%     byte_offset: (double) Offset to first image (bytes to skip)
%    byte_spacing: (double) Gap between images (bytes to skip)
%     frame_start: (double) Starting frame to import (default is 1)
%         frame_N: (double) Number of frames to import (detault is all)
%      frame_skip: (double) Frames to skip (1 will skip every other frame)
%   bit_10_packed: set to 'y' if data is 10-bit packed.
% bit_10_unpacked: set to 'y' if data is 10-bit unpacked.
%           bit_8: set to 'y' to output video to 8-bit
%
% OUTPUT:
%             vid: Output video (pix_h, pix_w, frame_N)
%
% CALLING SEQUENCE:
%   a = read_vif(filename, pix_w, pix_h)
%   a = read_vif('test.vif',656,491,byte_offset=72,byte_spacing=136,bit_10_unpacked='y')
%   a = read_vif('test.vif',700,700,byte_offset=8,byte_spacing=496, frame_start=100, frame_N=200)
%
% NOTES :
%   ImageJ - ImageJ (free software) is useful for determining the
%           parameters needed to read in a raw video format. With ImageJ
%           open, use File->Import->Raw.  Playaround with parameters as
%           needed.
%   VIF DataType - XCAP software saves video in a VIF format using three
%           possible data types: 8-bit unsigned, packed (10-bit unsigned
%           packed into bytes, where 5 bytes [40 bits] is 4 data points [40
%           bits]), and unpacked (10-bit unsigned interger saved in 16-bit
%           unsigned interger, so 6 bits are unused).
%   VIF pix_w,pix_h - In the .fmt file, saved with the VIF file,
%           is information about these values (xviddim and yviddim).
%   VIF byte_spacing - To determine byte_spacing, use the the frameIOSize
%           in the .ini file. To determine byte_spacing, subtract number of
%           bytes for each image from this vale. Example: 700x700 8-bit (1
%           byte) image with a frameIOsize of 490496 will result in
%           byte_spacing = 490496-700*700*1 byte or 496 bytes.
%
% REVISION HISTORY:
%   05/30/2013 - K Aptowicz
%       * Wrote original version in IDL
%   10/23/2022 - K Aptowicz
%       * Translated to MATLAB
%   06/08/2023 - K Aptowicz
%       * Cleaning up code and renaming variables
%
%% Reading and setting parameters
% Set default values for optional parameters
default_byte_offset = 72;  
default_byte_spacing = []; % If not set, calculated below. 
default_frame_start = 1;
default_frame_skip = 0;
default_frame_N = [];
default_packed = [];
default_unpacked = [];
default_bit_8 = [];

% Create fields for all optionals inputs
%Variables
p = inputParser;
addParameter(p,'byte_offset',default_byte_offset,@isnumeric)
addParameter(p,'byte_spacing',default_byte_spacing,@isnumeric)
addParameter(p,'frame_start',default_frame_start,@isnumeric)
addParameter(p,'frame_skip',default_frame_skip,@isnumeric)
addParameter(p,'frame_N',default_frame_N,@isnumeric)
% Keywords
addOptional(p,'bit_10_packed',default_packed)
addOptional(p,'bit_10_unpacked',default_unpacked)
addOptional(p,'bit_8',default_bit_8)

% populate optional parameters from inputs
parse(p,varargin{:})
byte_offset = p.Results.byte_offset;
byte_spacing = p.Results.byte_spacing;
frame_start = p.Results.frame_start;
frame_skip = p.Results.frame_skip;
frame_N = p.Results.frame_N;
bit_10_packed = p.Results.bit_10_packed;
bit_10_unpacked = p.Results.bit_10_unpacked;
bit_8 = p.Results.bit_8;

% Image size
m = pix_w;
n = pix_h;

% Calculating byte_spacing (padding)
% Frames are saved to disk sectors of 512 bytes
if isempty(byte_spacing)
    disp('Setting byte_spacing to pad frames to interger*512 bytes.')
    if ~isempty(bit_10_packed)
        byte_spacing = 512*ceil((m*n*5/4+72)/512)-m*n*5/4;
    elseif ~isempty(bit_10_unpacked)
        byte_spacing = 512*ceil((m*n*2+72)/512)-m*n*2;
    else
        byte_spacing = 512*ceil((m*n+72)/512)-m*n;
    end
    disp(['...... byte_spacing = ',int2str(byte_spacing),' bytes'])
end

% Read in file
finfo=dir(filename);
if size(finfo,1) == 0
    disp(['No files matched specification ',filespec])
else
    fid = fopen(finfo(1).name);
end

% Determine maximum number of frames
Nbytes=finfo.bytes;
if ~isempty(bit_10_packed)
    disp('Reading a packed video file was never tested!')
    if mod(m*n*5/4,1) ~= 0
        disp('Issue with bitpacking. pix_w*pix_h*10bit/8bit result in non-interger number of bytes.')
    end
    fmax=floor((Nbytes-byte_offset-m*n*5/4)./(m*n*5/4+byte_spacing))+1;
    fsize=m*n*5/4+byte_spacing;
    data_type = 'uint16';
elseif ~isempty(bit_10_unpacked)
    fmax=floor((Nbytes-byte_offset-m*n*2)./(m*n*2+byte_spacing))+1;
    data_type = 'uint16';
    fsize=m*n*2+byte_spacing;
else
    fmax=floor((Nbytes-byte_offset-m*n)./(m*n+byte_spacing))+1;
    data_type = 'uint8';
    fsize=m*n+byte_spacing;
end

% Determine last frame number (fmax) and number of frames
if ~isempty(frame_N)
    frame_end = frame_start + frame_N;
    frame_N = floor((frame_N-1)/(frame_skip+1))+1;
else
    frame_N = fmax-frame_start+1;
    frame_end=fmax;
    frame_N = floor((frame_N-1)/(frame_skip+1))+1;
end

% Goto first frame to be read in
fseek(fid, byte_offset+int32(frame_start-1)*fsize, 'bof');

%% Read in frames
vid=zeros(n,m,frame_N,data_type); % pre-allocation for video
for count=1:frame_N
    % Skip byte spacing between frames
    if ~isempty(byte_spacing) && count > 1
        fseek(fid,byte_spacing,'cof'); % Remove extra bytes between frames
    end
    % Skip frames if requested
    if (frame_skip~=0) && count > 1
        fseek(fid,fsize*frame_skip,'cof'); % Remove extra bytes between frames
    end
    if ~isempty(bit_10_packed)
        rawimg = fread(fid,m*n,'ubit10=>uint16');
        rawimg = reshape(rawimg,[m,n]);
    else
        rawimg = fread(fid,[m,n],data_type);
    end
    vid(:,:,count)=rawimg';
    count=count+1;
end

if ~isempty(bit_8) && (data_type ~= "uint8")
    vid=uint8(vid/4);
end

fclose(fid);

end
