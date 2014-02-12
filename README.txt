Requirements:
============
*Software that needs to be installed (if any) with URL's to download and instructions to install them.* 
Packages are required for R: igraph, compiler, doSNOW, foreach, parallel

_Environment variable settings (if any) and OS it should/could run on._
=======================================================================
None.

_Instructions on how to run the program._
=========================================
We need to change the path where the graph to be investigated lies (setwd). Update line setwd("<Full_Path_Where_Graph_Is_Present>"), for example
setwd("C:\\Users\\Shayan\\Desktop\\VirusPropagation").
For sample code run 
```````````````````````EXAMPLE``````````````````````````
1. unzip the file Virus_propagation
2. cd Virus_propagation/
3. update virusprop.R with full path where the static graph lies.
4. execute
````````````````````````````````````````````````````````

_Instructions on how to interpret the results._
===============================================
Graphs are generated for displaying the results for each section.

Sample input and output files. 
==============================
Output Files : will be generated once the code is executed.
Sample input : can be run on static.network directly

Tested On 
=========
static.network