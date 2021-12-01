# ttd
Interpolation tephra thickness distribution with detrending and kriging. 
The method presented here is used to interpolate the thickness distribution of tephra fall deposits from point measurements. It first reconstructs the general thickness decay pattern, or "trend," and then models the residuals that cannot be explained by the trend with ordinary kriging. The sum of these two is assumed to be the final thickness distribution of a tephra fall deposit. It is able to better characterize the influence of wind on tephra sedimentation with the help of a term, downwind distance, which is the distance to the source vent projected on to the wind dispersal axis. This method deals with log thickness data, and allows the plotting of the interpolated surface as an isopach map, and the resultant file is compatible with most GIS software.  See Yang and Bursik (2016) for a complete description:

Yang, Qingyuan, and Marcus Bursik. "A new interpolation method to model thickness, isopachs, extent, and volume of tephra fall deposits." Bulletin of Volcanology 78.10 (2016): 68.

The downloadable file is the R code that executes this method.  The way to infer the extent of tephra fall deposits is done by a "plus-one" transformation, which is introduced in the paper (Yang and Bursik, 2016), and therefore there is no specific function; a brief description is mentioned in the notes of the R script.

The files include:

  (1) The r scripts for this method that includes a short demo of how to implement it.

  (2) The thickness dataset of the North Mono Bed 1 tephra fall deposits. The dataset is digitized from the work of Sieh and Bursik (1986). "x" and "y" stand for the coordinates of sample sites, "rz" stands for the real thickness measured in millimeters, and "z" stands for the logarithm (base=10) of "rz".

  (3) The thickness dataset of the Fogo A deposits, originally and generously provided by Dr. Samantha Engwell. The dataset is digitized from the work of Walker and Croasdale (1971). It has four columns, "x" and "y" stand for the coordinates of the sample sites, "rz" stands for the measured thickness in millimeter. "notes" are for the sample sites with specific notes: LT and GT stand for "the real thickness should be less than or greater than the measurement". These notes are based on the original fieldwork (Walker and Croasdale, 1971). The user can delete or include these uncertain measurements for the isopach mapping, and then calculate the "z", i.e., the log-transformed thickness, to make the input thickness dataset compatible with the R script. The sample sites with no deposits found have to be deleted because of the log-transformation.

The source vent locations and wind directions that we used are (in utm):

  (322992, 4196119) for North Mono bed 1 deposit, with in put wind direction 20;

  (634109, 4180765) for Fogo A deposit, with input wind direction 160.

Users are welcomed to try different source vent locations and wind directions based on the knowledge on these deposits and personal preferences.  We hope the method will be helpful in others'  work. If you have any questions, comments or suggestions on this code and the method itself, please feel free to contact the authors.

-------------------------------------------------------

Frequently Asked Questions (FAQ):

Question One: How do I work with the method if (1) my measurements were located in different UTM zones or (2) the data are in geographic coordinates?

Answer:  The algorithm requires the coordinates of the input sample sites to be in the same projected coordinate system (such that the unit is in meters or kilometers). Therefore, users should transform the coordinate system of the sample sites before using the algorithm. For example, if you have measurements of a tephra fall deposit with sample sites located within different UTM zones in the U.S., you could transform them into the Albers Equal Area Conic coordinate system (http://desktop.arcgis.com/en/arcmap/latest/map/projections/albers-equal-area-conic.htm).  This could be done by any spatial data analysis program or package. Since Arcmap is the most commonly-used and presumably the most accessible spatial data analysis program for universities and research institutes, we suggest the following procedures to deal with this issue:

    For the first case, we recommend you to:

     (a) Import your data including UTM zone into Arcmap using Add XY coordinates in Arcmap;

            (http://desktop.arcgis.com/en/arcmap/10.3/map/working-with-layers/adding-x-y-coordinate-data-as-a-layer.htm)

     (b) Transform your data points with different UTM zones into a projected coordinate system that covers the whole sampled area;

            (http://desktop.arcgis.com/en/arcmap/10.3/map/working-with-arcmap/specifying-a-coordinate-system.htm#GUID-7D31FB5A-FFEF-4CCF-AB6C-A745DD138C1B)

     (c) Add the new coordinates to your data table/points/feature class (https://pro.arcgis.com/en/pro-app/tool-reference/data-management/add-xy-coordinates.htm).

     (d) Export your data table/points/feature class (which should include the old coordinates, your thickness or maximum clast size measurements, and the new transformed coordinates):

            d1. Right click the point shape file from the Table of Contents/layer tab, and click "open attribute table";

            d2. Click "table options" tab in upper right of Table;

            d3. Export.

For the second case, you could just transform the geographic coordinates to a projected, metric coordinate system that covers you whole sampled area following the same procedure. 

-------------------------------------------------------

    *** We are glad to receive our first question from the users. Any feedback, questions, comments, and suggestions are always welcome! 

References:

Sieh, Kerry, and Marcus Bursik. "Most recent eruption of the Mono Craters, eastern central California." Journal of Geophysical Research: Solid Earth91.B12 (1986): 12539-12571.

Walker, George Patrick Leonard, and Ronald Croasdale. "Two plinian-type eruptions in the Azores." Journal of the Geological Society127.1 (1971): 17-55.

Yang, Qingyuan, and Marcus Bursik. "A new interpolation method to model thickness, isopachs, extent, and volume of tephra fall deposits." Bulletin of Volcanology 78.10 (2016): 68.

