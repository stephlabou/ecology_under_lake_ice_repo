#Packages required

The following packages are needed to run the scripts in this repository:

####Data organization/formatting:
data.table <br>
reshape2 <br>
plyr <br>
dplyr <br>
tidyr <br>
devtools <br>

####Time series and LME:
zoo <br>
forecast <br>
AICcmodavg <br>
nlme <br>
MuMIn <br>
pls <br>

####Regression trees:
rpart <br>
rpart.plot <br>
mvpart <br>

####Other analysis:
exactRankTests <br>
robustbase <br>
vegan <br>

####Graphing:
colorspace <br>
ggplot2 <br>
RColorBrewer <br>
grid <br>
gridExtra <br>
scales <br>
plotrix <br>

####Maps:
maps <br>
mapdata <br>
mapproj <br>

#Package notes

The mvpart package can be downloaded from the CRAN archive: https://cran.r-project.org/src/contrib/Archive/mvpart/
Once downloaded, the package can be installed using: `install.packages(".../Downloads/mvpart_1.6-2.tar.gz", repos = NULL, type = "source")`. 
For help downloading/installing mvpart, see the stackoverflow page: http://stackoverflow.com/questions/29656320/r-mvpart-package-any-option-to-use-in-r-3-1-x  

Please use caution when loading packages, as some may conflict and/or mask calls from one another (i.e. in cases where the same call name is in multiple packages). 
We recommend running each script individually and restarting the R session before running the next script. This ensures that function calls within the script are 
called from the correct library and all analyses work as intended

#System requirements

Scripts included in this repository were developed on a variety of systems. All scripts have been tested and run as expected on the following system with the following package versions:

###Output of sessionInfo()
####version
R version 3.2.5 (2016-04-14) <br>
Platform: x86_64-w64-mingw32/x64 (64-bit) <br>
Running under: Windows >= 8 x64 (build 9200) <br>

####locale:
LC_COLLATE=English_United States.1252  
LC_CTYPE=English_United States.1252    
LC_MONETARY=English_United States.1252 
LC_NUMERIC=C                          
LC_TIME=English_United States.1252    

####attached base packages:
grid      
stats     
graphics  
grDevices utils     
datasets  
methods   
base     

####other attached packages:
mvpart_1.6-2          
rpart.plot_2.0.1      
rpart_4.1-10          
vegan_2.4-0           
lattice_0.20-33       
permute_0.9-0         
robustbase_0.92-6     
exactRankTests_0.8-28
pls_2.5-0             
MuMIn_1.15.6          
nlme_3.1-128          
AICcmodavg_2.0-4      
forecast_7.1          
timeDate_3012.100     
zoo_1.7-13            
devtools_1.12.0      
tidyr_0.5.1           
dplyr_0.5.0           
plyr_1.8.4            
reshape2_1.4.1        
data.table_1.9.6      
mapproj_1.2-4         
mapdata_2.2-6         
maps_3.1.0           
plotrix_3.6-3         
scales_0.4.0          
gridExtra_2.2.1       
RColorBrewer_1.1-2    
ggplot2_2.1.0         
colorspace_1.2-6     

####loaded via a namespace (and not attached):
splines_3.2.5 
stats4_3.2.5    
mgcv_1.8-13     
chron_2.3-47    
survival_2.39-5 
withr_1.0.2     
DBI_0.4-1       
sp_1.2-3        
stringr_1.0.0   
munsell_0.4.3   
gtable_0.2.0   
raster_2.5-8    
VGAM_1.0-2      
memoise_1.0.0   
tseries_0.10-35 
parallel_3.2.5  
DEoptimR_1.0-6  
Rcpp_0.12.6     
xtable_1.8-2    
unmarked_0.11-0 
fracdiff_1.4-2  
digest_0.6.9   
stringi_1.1.1   
quadprog_1.5-5  
tools_3.2.5     
magrittr_1.5    
tibble_1.1      
cluster_2.0.4   
MASS_7.3-45     
Matrix_1.2-6    
assertthat_0.1  
reshape_0.8.5   
R6_2.1.2       
nnet_7.3-12    
