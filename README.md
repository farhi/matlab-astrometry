# matlab-astrometry
  astrometry: A Matlab class to annotate astrophotography images (identify objects/astrometry)
  
  ![Image of Astrometry](https://github.com/farhi/matlab-astrometry/blob/master/examples/M13-solved.jpg?raw=true)
 
Purpose
-------
 
This Matlab class allows to use the astrometry.net software, either installed 
locally, or through internet connection, in order to solve (annotate) 
astrophotography images.
 
Syntax/Usage
------------
 
First navigate to the matlab-astrometry directory or type:
 
   ```matlab
   addpath /path/to/matlab-astrometry
   ```
   
Then use:
 
   **as = astrometry;**
   
Create a solver, but does not solve. 
Use annotation(as, file) or web(as, file) afterwards.
 
   **as = astrometry(file, ...); image(as);**
   
Solve the given astrophotography image with local or web method. 
Then plot the result. Additional arguments may include name/value pairs
(see example below):
 
       - ra:      approximate RA coordinate  (e.g. deg or  'hh:mm:ss')
       - dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
       - radius:  approximate field size     (in deg)
       - scale-low:   lower estimate of the field coverage (in [deg], e.g. 0.1)
       - scale-high:  upper estimate of the field coverage (in [deg], e.g. 180)
 
These two syntaxes will try first any local astrometry.net installation, and
if failed, the http://nova.astrometry.net/ service.

   **image(as)**
   
Plot the astrometry solution. The identified objects are indicated in green for 
stars, cyan circle for extended deep sky objects, and cyan squares for localized
deep sky objects. Each indicated object has a contextual menu which gives more 
information (right click). The central coordinate of the field is shown in red.
 
Going further
-------------
 
   **as = astrometry.load(dir); image(as);**
   
Read an existing Astrometry.net set of files stored in a given directory.
The directory may contain WCS, CORR, RDLS, JSON, and image.
Then plot the result. This allows to get previous data files, or obtained
externally, and label them. The 'as' astrometry object must have been used
to solve or import astrometry data.
 
   **[x,y] = sky2xy(as, ra, dec)**
   
Convert a RA/DEC set of coordinates (in [deg] or 'hh:mm:ss'/'deg::mm:ss')
into pixel coordinates on the image. The 'as' astrometry object must have 
been used to solve or import astrometry data.
 
   **[ra, dec] = xy2sky(as, x,y)**
   
   **[ra, dec] = xy2sky(as, x,y, 'string')**
   
Convert pixel coordinates on the image into a RA/DEC set of coordinates 
(in [deg]). When given a 'string' argument, the result is given in 
'hh:mm:ss'/'deg:mm:ss'. The 'as' astrometry object must have been used
to solve or import astrometry data.
 
   **f=astrometry.findobj('object name')**
   
Return information about a named object (star, deep sky object) from the 
data base. Example: astrometry.findobj('M33')
 
   **as = astrometry.annotation(file, ...);**
   
Explicitly use the local 'solve-field' astrometry.net installation.
See above for the additional arguments.
 
   **as = astrometry.web(file, ...);**
   
Explicitly use the http://nova.astrometry.net/ web service.
See above for the additional arguments.

**as.result.RA** and **as.result.Dec** provide the center coordinates of the 
field (in [deg]), while **as.result.RA_hms** and **as.result.Dec_dms** provide 
the 'HH:MM:SS' and 'Deg:MM:SS' coordinates. The field rotation wrt sky is stored
in **as.result.rotation**. The pixel scale is given in [arcmin/pixel] 
as **as.result.pixel_scale**. The field extension is given with its bounds as
**as.result.RA_min**, **as.result.RA_max**, **as.result.Dec_min**, and 
**as.result.Dec_min**. The constallation name is stored in **as.result.Constellation**.
 
Improving the plate-solve efficiency
------------------------------------
 
To facilitate the plate-solve/annotation of images, you may:
 
* specify the field size with additional arguments such as: 
  astrometry(..., 'scale-low', 0.5, 'scale-high',2)

* provide an initial guess for the location, and its range, such as:
  astrometry('examples/M33-2018-05-19.jpg','ra','01:33:51','dec','30:39:35','radius', 2)

* add more star data bases (e.g. 2MASS over Tycho2).
 
Examples
--------
 
```matlab
   as=astrometry('examples/M33-2018-08-15.jpg','scale-low', 0.5, 'scale-high',2);
   image(as);
```

  You will then get, in about 30 sec, 

Methods
-------
   
    - as=astrometry(filename)
    - image(as)
    - load(astrometry, dir)
    - annotation(astrometry, filename, ...)
    - web(astrometry, filename, ...)
    - sky2sx(as, ra, dec)
    - xy2sky(as, x, y)
    - findobj(as, 'name')
 
Installation
------------
 
   **Local installation (recommended)**
 
On Linux systems, install the 'astrometry.net' package, as well as the 'tycho2' data base. On Debian-class systems, this is achieved with:
     
```bash
        sudo apt install astrometry.net astrometry-data-tycho2 sextractor
```

On other systems, you will most probably need to compile it.
See: http://astrometry.net/doc/build.html
RedHat/Arch/MacOSX have specific installation instructions.
 
If you have images spanning on very tiny areas (e.g. much smaller than a 
degree), you will most probably need to install the '2MASS' data base.
 
   **Using the web service**
 
 You will need Python to be installed, and to have a 'NOVA astrometry API' key.
 Enter the API_KEY when prompt, or set it with:
 
 ```matlab
     as = astrometry;
     as.api_key = 'blah-blah';
     as.web(file, ...)
 ```
 
   **Matlab files**
   
First navigate to the matlab-astrometry directory or type:
 
   ```matlab
   addpath /path/to/matlab-astrometry
   ```
 
  Credit
  ------
 
**sky2xy** and **xy2sky** from E. Ofek http://weizmann.ac.il/home/eofek/matlab/
 
  (c) E. Farhi, 2018. GPL2.
