IMCORR Release 1.1

------------------------------------------------------------------------
2020-12-21

This is unsupported software.  Information below, particularly contact
information, may be out-of-date.
------------------------------------------------------------------------

(Note:  This file, the makefiles, and a couple of other details about the
distribution have been slightly edited and put under revision control by Bruce
Raup.  The software is functionally unchanged (except for one minor
compile-time option), but the release number is no longer 1.0.  If you have
comments or suggestions for changes or additions, please send them to B.
Raup and T.  Scambos.)

Notes on using IMCORR image-to-image greylevel correlation software

T. Scambos       (1)
M. Fahnestock    (2)
R. Bindschadler  (3)
B. Raup          (4)

(1) Cooperative Institute for Research in Environmental Sciences
University of Colorado Boulder   CO      80309-0449      USA
tascambos@Colorado.EDU

(2) Joint Center for Earth System Science, Department of Meterology,
University of Maryland, College Park, MD  USA
voice: 301/405-5384 net:mark@firn.gsfc.nasa.gov

(3)Code 971 NASA/Goddard Space Flight Center  Greenbelt MD 20771   USA
Voice:  301/286-4717 or 301/286-2142     Fax:  301/286-2717
Internet:   iceflow@laural.gsfc.nasa.gov (IP:  128.183.47.79)

(4) National Snow and Ice Data Center, CIRES; Campus Box 449; 1540 30th St.
University of Colorado Boulder   CO      80309-0449      USA
303-492-8814
braup@nsidc.org


INTRODUCTION

The following documentation is intended to help interested users acquire
and begin using IMCORR, a simple image cross-correlation program which can
be used to find offsets of features in a pair of images.  This software
uses subroutines found within LAS 5.0, an extensive image processing
package written in FORTRAN at NASA/Goddard.  IMCORR consists of these
subroutines and a simple C routine that acts as a front end.  All the
subroutines required for displacement analysis are within the directory
associated with this document file.

IMCORR

IMCORR takes two images and a series of input parameters and attempts to
match small subscenes (called 'chips') from the two images.  The program
uses a fast fourier transform - based version of a normalized
cross-covariance method (see Berenstein, 1983; also Scambos et al., 1992).
The most common use of this type of algorithm in image processing is to
accurately locate tie-point pairs in two images to coregister them.
However, if the images are already coregistered by other means, the
algorithm may be used to find the displacements of moving features,
provided that the features show little change in their appearance, and that
the motion is strictly translational.  IMCORR takes as input the image
names and sizes, parameters determining search chip size, reference chip
size, grid spacing, and output filename.  Further, preset offsets of search
chip centers may be specified, and subareas of the full image files may be
used to restrict the area over which IMCORR attempts to find displacements.
At each of the gridpoints IMCORR calculates a correlation index for every
location at which the reference chip will entirely fit within the search
chip.  IMCORR takes the correlation values in the vicinity of the best
integer-pixel match and interpolates a peak correlation location to
sub-pixel precision.  The program returns a file containing the locations
of the grid centers for the reference chips, the displacements required to
best match the chip pairs (or indicates that none could be found), and
several quality control parameters that may be used to evaluate the
validity of the match.  We use this program to measure glaicer velocities;
however, the same program may be useful for other applications.

It is important for the reader to understand that all of the correlation,
peak finding, and error estimation is done by the FORTRAN subroutines from
LAS, which was written by others at NASA long before the authors of this
document were involved in image analysis.  We have simply added the C code
which makes the use of these routines more straightforward and automated
for velocity-mapping applications.

A paper describing the details of this technique and the preprocessing of
images required to optimize the displacement measurements is published in
Remote Sensing of Environment (Scambos et al., 1992; see bibliography).  We
will gladly provide a preprint upon request.

GETTING THE SOFTWARE

The latest distribution for imcorr is at
http://www-nsidc.colorado.edu/NSIDC/ANTARCT_VELOC/imcorr.html.

LOADING THE SOFTWARE

Transfer the entire contents of the 'distribution' directory to your
workstation (the tutorial discussion below will assume you have transferred
the four test images we have included in the directory).  This directory is
roughly 1.25 Mbytes, and contains the following files:

conv.y87         test image file
conv.y89         test image file
correlate.h      Header file for imcorr.c
cross.f          FORTRAN source
esterr.f         FORTRAN source
eval.f           FORTRAN source
fang.y87         test image file
fang.y89         test image file
fft2d.f          FORTRAN source
fitreg.f         FORTRAN source
gcorr.f          FORTRAN source
gnorm.f          FORTRAN source
imcorr.c         C source
kvert.f          FORTRAN source
makefile.Linux   make file for Linux machines
makefile.SGI     make file for Silicon Graphics machines
makefile.SUN     make file for Sun Microsystems machines
README.txt       this file
sums.f           FORTRAN source
test1.out        Example output file
test2.out        Example output file
test3.out        Example output file
TO_DO.txt        list of to-do's

COMPILE-TIME OPTION

Before compiling, you can uncomment (delete the '#' character) the line in
the makefile that says

#DEFS = -DCOMMANDLINE

This will cause imcorr to print the command-line that was issued at the
shell prompt as the first line in the output file, preceded by a '#'
character.

HOW TO COMPILE

The IMCORR source code is roughly 65 kilobytes.  It requires both C and FORTRAN
compilers.  It contains no graphics commands in either language, and should be
relatively easy to port to a variety of systems.  We use Silicon Graphics
workstations; the software has also been compiled on SUN Sparcstations and
Linux machines.  We have included makefile files for SGI, Linux, and SUN
workstations.  To generate an executable file from the source code, type:

prompt>make -f  makefile.[SUN or SGI or Linux] imcorr

On any system, to remove the object files, type:
prompt>make -f makefile.[SUN or SGI or Linux] clean

If you want or need to make any changes to the makefile, we recommend:

prompt> cp makefile.[SUN or SGI or Linux] Makefile
prompt> vi Makefile (make changes with whatever editor you like)
prompt> make

By default, make uses "Makefile" as input, if it exists.

All of these makefiles produce an executable file named imcorr.

INPUT PARAMETERS

This IMCORR version takes either 5, 10, or 14 input parameters, which are
typed on the same line as the 'imcorr' command separated by spaces.  Simply
typing 'imcorr' and return gives a list of the required input parameters
for each of these three formats.

	List of Input Parameters:
reference image filename (generally the older image)
search image filename (generally the newer image)
number of pixels per line in images (no default)
number of lines in images (no default)
output filename (does not need to exist prior to a run)
search chip size in pixels
      (default : 64; must be divisible by 16, less than or equal to 256;
      entering '32' means a search chip of 32 x 32 pixels)
reference chip size in pixels
      (default : 32; must be divisible by 16, less than or equal to 128;
      entering '32' means a ref chip of 32 x 32 pixels)
grid spacing in pixels
      (default: 25; no limits except those imposed by the images;
      grids are always square; entering '25' means that the search and
      reference image chips will be centered on points 25 pixels apart
      in both line and sample directions)
x offset of reference chip relative to search chip center (default: 0)
y offset of reference chip relative to search chip center (default: 0)
x offset of subimage to be analysed (default: 0)
y offset of subimage to be analysed (default: 0)
x extent of subimage (default: full x extent)
y extent of subimage (default: full y extent)

The reference image will have the smaller, 'reference chip', subimages derived
from it, at regular grid spacings, and these chips will be compared to larger
'search chip' subimages derived from the search image.  The  images must be the
same size, and a single (sample, line) coordinate system will be used for both
of them.  In the following discussion, it is assumed that the images are
coregistered, i.e., that a feature that does not move would be found at the same
(sample, line) coordinates in both images.  Regarding maximum image size,
machine memory is dynamically allocated during an IMCORR run.  Therefore, images
may be very large; however, larger images will result in slower processing due
to memory swapping.

The x and y offsets of the subimage answer the question, "where does the
reference chip come from?".  For example, if the ice has moved in the positive
x direction (to the right) by 20 pixels, then the x-offset should be -20.

The coordinate system has x positive to the right, and y positive downward:

    --------------> x
    |
    |
    |
    |
    |
    |
    V
    y

Large search and reference chip sizes slow the runtime considerably, and a
large reference chip may cause problems if any distortion of the features
occurs between the two images.  Note that if the search and reference chips are
derived from the same (sample, line) grid centers in their respective images
(the default case), the maximum displacement of a feature that can be measured
is: (search chip size/2 - reference chip size/2) * C2, and this maximum applies
only to exactly diagonal motion.

  Grid spacing parameter determines the density of attempted matches, and
therefore also has a strong effect on runtime.  Note that reducing the grid
spacing parameter by a factor of 2 increases the number of points attempted by
4.  The grid of attempted matches will be more or less centered within the
reference image area; since search chips and reference chips are derived from
the same (sample, line) centers in the default case, reference chips are derived
from centers no closer than (search chip size/2) from the reference image edge.
In the search image, matches may be found as close as (reference chip size/2)
from the search image edge.
	For a variety of reasons, it is often desirable to give a preset offset
of the search chip center relative to the reference chip center.  For example,
in cases where the displacement is approximately known beforehand, one may reduce
the size of the search chip if a preset offset of the search chip area in the
approximate direction of the displacement is entered.  This is the purpose of
the next two parameters.  The offsets are in image coordinates, relative to the
SEARCH chip; therefore, entering a positive value for the x, or sample, offset
coordinate will place the center of the reference chip to the right of the
search chip (by moving the search chip center leftward).  A positive value for
the y, or line, offset will place the center of the reference chip down, or in
the positive line direction.
	A subarea of the full images may be specified, limiting the search for
displacement matches to a specific region within the scenes.  The run proceeds
as if the subareas were the full image size.  However, the output locations of
any matches are given relative to full scene coordinates.

THE OUTPUT PARAMETERS

The output file is in ASCII format, and consists of a series of lines with nine
entries each.  A line is generated for every gridpoint, whether a successful
match was found or not, and no matter what the quality of the match was.

List of Output Parameters:

x coordinate of center of reference chip	(integer)
y coordinate of center of reference chip	(integer)
total displacement in pixels	(real)
strength of correlation (arbitrary units)	(real)
result flag (no units)	(integer)
x displacement to best match within search chip in pixels	(real)
y displacement to best match within search chip in pixels	(real)
x estimated error in pixels	(real)
y estimated error in pixels	(real)

	The first two parameters give the location of the center of the reference
chip.  The total displacement is exactly the distance formula applied to the
x and y displacements listed later; it is zero if no successful match was found.

	The correlation strength parameter is a function of the character of
a correlation index 'map' that the greycorr subroutines create and evaluate -
it is a combination of peak height of the correlation map, height of peak to
second-highest peak, and height of peak to background value of correlation
index.  You want the formula? OK, here:
correl. strength = (peak correl. value - mean background value)/
		std.dev. of background values  +
(peak correl. value - highest value more than 3 pixels away from peak)/
		std.dev. of background values  +
0.2 * (number of "large" values more than 3 pixels away from peak - 1.0)

Its magnitude changes somewhat with the input parameters for any given run -
larger search and reference chips tend to produce larger strength of correlation
values for a given region of the images.

	The result flag parameter is an integer returned by imcorr indicating
if a good match was found, or if not, what went wrong.  Flag values mean:
 1: successfully correlated
 2: failed to correlate; Correlation peak too near edge of search chip
 3: failed to correlate; Subsidiary peak comparable in strength to
	main peak
 4: failed to correlate; Strength of peak below minimum specified by
	user (we have set minimum to zero, so this one will not	appear)
 5: failed to correlate;  Diagonal displacement from nominal location
	exceeds maximum specified by the user
Case 2 is indicated if the match was found within 2 pixels of the limit of where
the reference chip can fit within the search chip - such proximity to the edge
does not allow for good statistical determination of whether the match is valid.
If case 2, 3, or 4 occurs, zero values are output for peak strength,
displacements, and errors.
	The next two output parameters are the offsets required to best match the
reference chip pixels with the search chip pixels, again in image coordinates
(positive y is downward).  Preset offsets in the input parameters are included
in the reported displacement measurement.   The x and y error estimates are
derived from a peak-height-to-peak-width comparison.  The values are sensitive
to the size of the reference chip; good matches with larger reference chips
yeild smaller x and y errors, in general.

RUNNING A TEST ANALYSIS

	Within the same directory as the IMCORR software there are four images
that have been preprocessed for IMCORR analysis, and a series of results files
labeled test1.out, test2.out, and test3.out.  (The first lines of the
example output files contain the command line.)  The images, called
conv.y87, conv.y89, fang.y87, and fang.y89, are Landsat TM subscenes of a
portion of Ice Stream D in West Antarctica (roughly 81! 20' S, 135! W).
They are byte images 512 x 512 pixels in size.  The image pairs are already
coregistered relative to basement features and are high-pass filtered to
enhance small, sharp ice surface features that move with the ice's surface
velocity.  See the preprocessing notes below for more information on
preparing images for IMCORR.  In the sample runs suggested below, we will
use these images as input for IMCORR.

The following is an example of an IMCORR run over an area where a high
percentage of the matches are 'correct'.  It uses all the default parameters.

prompt>imcorr conv.y87 conv.y89 512 512 test1.out

you will see the following as the run proceeds:

opening conv.y87 as ref image and conv.y89 as search image, size 512 by 512
writing to test1.out
xoff=0 yoff=0 xext=512 yext=512
opening conv.y87 size 512 by 512
reading conv.y87
done reading image
read 262144 pixels
opening conv.y89 size 512 by 512
reading conv.y89
done reading image
read 262144 pixels
l 64 r 64 t 64 b 64
64	64	8.582	3.793	1	0.753	8.549	0.093	0.065
64	89	0.000	0.000	3	0.000	0.000	0.000	0.000
64	114	9.233	6.515	1	1.842	9.048	0.589	0.045
64	139	9.785	5.911	1	0.713	9.759	2.193	0.128
.
.
.
439 	439	29.983	11.734	1	-13.120	21.260	0.157	0.197

The output beginning with the first 64, 64 line will go into a file labeled
test1.out.  This data may then be edited in a number of ways, using the strength
parameter, result flag, x and y errors, etc.,  to remove bad or poor quality
points.

Another example of a IMCORR run, using the fang image pair :

prompt>imcorr fang.y87 fang.y89 512 512 test2.out 128 32 20 0 0

opening fang.y87 as ref image and fang.y89 as search image, size 512 by 512
writing to test2.out
xoff=0 yoff=0 xext=512 yext=512
opening fang.y87 size 512 by 512
reading fang.y87
done reading image
read 262144 pixels
opening fang.y89 size 512 by 512
reading fang.y89
done reading image
read 262144 pixels
l 64 r 64 t 64 b 64
64	64	0.000	0.000	3	0.000	0.000	0.000	0.000
64	84	23.411	7.392	1	-11.579	20.347	4.308	4.425
64	104	0.000	0.000	3	0.000	0.000	0.000	0.000
64	124	23.283	10.332	1	-10.797	20.628	0.388	0.219
.
.
.
444	444	5.912	116.630	1	-3.613	4.680	8.011	4.871

This area is more difficult for a number of reasons. You will note that a
diagonal region (corresponding with a chaotically fractured zone in the image)
has no successful matches.  This is due to high shear in this area.  Further,
there are some areas where relatively stationary shadowing has caused matches
which are displaced relative to the probable direction of flow.  This area would
require careful editing, using the quality parameters written in the output
file, as well as hand editing by inspection.

An example of a run which attempts to get data in the shear area is


prompt>imcorr fang.y87 fang.y89 512 512 test3.out 32 16 8 5 -11 50 0 200 512

This run has a reduced reference chip size, thus encompassing less velocity
gradient, and uses a reduced search chip area and offsets to select a smaller,
more likely target area.  A denser grid spacing gives more chances that usuable
matches will be found.  It also uses the subimaging input parameters to limit
the region over which this parameter set will be applied (over most of the
image, this parameter set will yield poorer results than the first - only in the
high-strain region will it be better).

The results on the screen from the second run should be:

opening fang.y87 as ref image and fang.y89 as search image, size 512 by 512
writing to test3.out
xoff=50 yoff=0 xext=200 yext=512
opening fang.y87 size 512 by 512
reading fang.y87
done reading image
read 102400 pixels
opening fang.y89size 512 by 512
reading fang.y89
done reading image
read 102400 pixels
l 24 r 8 t 8 b 34
71	8	0.000	0.000	2	0.000	0.000	0.000	0.000
71	16	0.000	0.000	3	0.000	0.000	0.000	0.000
71	24	0.000	0.000	3	0.000	0.000	0.000	0.000
71	32	0.000	0.000	3	0.000	0.000	0.000	0.000
.
.
.
239	480	0.000	0.000	2	0.000	0.000	0.000	0.000

(other lines had successful matches)

The above 'style' of running IMCORR should be used with caution, as it presumes
that a specific, narrow, range of displacements is occurring in the analyzed
area.

PREPROCESSING OF IMAGERY

Since most users of this software will be concerned with the measurement of
ice velocity using digital satellite imagery, or, perhaps, digitized aerial
or satellite photographs, this section will address preprocessing from that
standpoint.

The objective of the preprocessing is to generate images of the surface features
of  the ice, with as little noise, sensor effects, and solar illumination
effects as possible.  To have the surface features appear as similar as possible
in the sequential imagery, the imagery should be taken at the same time of day,
and if possible, at same time of year.  With multiband data, such as Landsat, a
good first step is to generate a first principal component image of the visible
and near-infrared bands.  This will greatly reduce the noise and has the effect
of giving greater brightness resolution to the image (see Orheim and Lucchitta,
1987).  Scan-line striping and swathing, if present, should be removed, as its
presence will tend to generate matches based on the stripes and not the
features.  One technique which should be applicable to a number of striping
problems is discussed in Crippen, 1989.  Sun angle variations across the image
may be addressed by dividing the pixel brightness values by the cosine of the
sun elevation as it varies across the image area, or by simply high-pass
filtering.  High-pass filtering may also be used to remove bed-related
topographic features which remain stationary as the ice flows over them.  The
presence of such features will distort the displacements measured because the
IMCORR routine will attempt to shift matches based on the surficial features so
that the larger-scale pixel brighness variations associated with the bed
topography also tend to match.

BIBLIOGRAPHY

Bernstein, R. (1983), Image geometry and rectification, In Manual of Remote
Sensing (R. N. Colwell, ed.), American Society of Photogrammetry, Falls Church,
VA, pp.881-884.

Crippen, R. E. (1989), A simple filtering routine for the cosmetic removal of
scan-line noise from Landsat TM P-tape imagery, Photogrammetric Engineering and
Remote Sensing, 55, 327-331.

Orheim, O., and Lucchitta, B. K. (1987), Snow and ice studies by Thematic Mapper
and Multispectral Scanner Landsat images, Annals of Glaciology, 9, 109-118.

Scambos, T. A., M. J. Dutkiewicz, J. C. Wilson, and R. A. Bindschadler, 1992.
Application of image cross-correlation to the measurement of glacier velocity
using satellite image data.  Remote Sensing of Environment, vol. 42, 177 - 186.
