# Particle-Tracking-using-Matlab
MATLAB code for analyzing colloidal particles captured with video microscopy. Based on IDL code written by David Grier (NYU), John Crocker (UPenn), and Eric Weeks (Emory). 

## Part 1: Reading in video file - reaf_vif.m
In our research lab, we use XCAP software by EPIX to grab video files (VIF). These can be converted into images using XCAP, but the Matlab function 'read_vif.m' can also read these VIF files, although it often needs to be tweaked. Read the more header of read.vif for more information.  

Examples of use **read_vif.m**: 
```
a = read_vif(filename, pix_w, pix_h)
a = read_vif('test.vif',656,491,byte_offset=72,byte_spacing=136,unpacked='y')
a = read_vif('test.vif',700,700,byte_offset=8,byte_spacing=496, frame_start=100, frame_N=200)
```
 INPUT (REQUIRED):
 ```
        filename: (string) filename of video to be imported
           pix_w: (double) pixel width of video frames
           pix_h: (double) pixel height of video frames
```
INPUT (OPTIONAL):
```
     byte_offset: (double) Offset to first image (bytes to skip
    byte_spacing: (double) Gap between images (bytes to skip)
     frame_start: (double) Starting frame to import (default is 1)
         frame_N: (double) Number of frames to import (detault is all)
      frame_skip: (double) Frames to skip (1 will skip every other frame)
   bit_10_packed: set to 'y' if data is 10-bit packed.
 bit_10_unpacked: set to 'y' if data is 10-bit unpacked.
           bit_8: set to 'y' to output video to 8-bit
``` 
Viewing image in Matlab: 
```
imagesc(a(:,:,*))
```
