#This module contains the classes and functions for manipulating MODIS land cover data.
#It is highly recommended that users leave the MODIS file naming convention in place.
#(I.e. don't rename your files.)
#For full documentation of the dataset see the MCD12Q1 and MCD12C1 data pages at LPDAAC.
#
#Written by Christopher Phillips
#University of Alabama in Huntsville
#Department of Atmospheric and Earth Sciences
#
#This code may be distributed freely provided this header remains intact during distribution.
#
#Requirements:
#Numpy       For array support
#PyHDF       For HDF reading
#PyProj      For the map projection
#
#Importing required modules
from matplotlib.colors import ListedColormap
from pyhdf import SD
from pyproj import Proj
import numpy

#This class is for accessing the MODIS MCD12Q1 product
#This class has one method: self.get()
#self.get(varname) is used to retreive numpy arrays from the HDF file
#Example -
#
#from pymodis.modislc import MCD12Q1
#datafile = MCD12Q1(path_to_modis_file)
#landcover = datafile.get("LC_Type1")
#
#/Example
#
#Class attributes are:
#nhtiles - number of MODIS tiles in the horizontal globally
#nvtiles - number of MODIS tiles in the vertical globally
#pixel_size - edge length of MODIS pixel in meters
#npixels - number of pixels along edge of the MODIS tile
#proj - dictionary containing MODIS map projection information
#tile - tuple containing MODIS tile number (Horz, Vert)
#hfile - HDF file object
#structure - HDF file structure string
#meta - HDF file metadata string
#lons - 2D array of longitudes corresponding to MODIS pixel locations [-180, 180]
#lats - 2D array of latitudes corresponding to MODIS pixel locations [-90, 90]
#lcN_legend - List of strings corresponding to land cover classes for classification scheme N
#lpN_legend - List of strings corresponding to land cover properties for classification scheme N
#qc_legend - List of strings corresponding to quality control flags for the dataset
class MCD12Q1():
    ### Initialization function that attaches HDF file object to class object
    #It also defines attributes with land cover classification legends
    def __init__(self, filepath):
        ### Creating attributes with metadata from the documentation
        self.nhtiles = 36
        self.nvtiles = 18
        self.pixel_size = 463.312716525
        self.npixels = 2400
        self.proj = {"proj":"sinusoidal", "a":6371007.181, "b":6371007.181, "units":"meters"}        
        
        ### Locating MODIS tile (H/V)
        #This depends on the file naming convention
        self.tile = (int(filepath[-27:-25]), int(filepath[-24:-22]))
                
        ### Create file reference and set as attribute
        self.hfile = SD.SD(filepath, SD.SDC.READ)
        self.structure = self.hfile.attributes()["StructMetadata.0"]
        self.meta = self.hfile.attributes()["CoreMetadata.0"]
        
        ###Creating lat/lon grid and storing as attributes
        #Creating projection object for converting the data to lat/lon
        modis_grid = Proj("+proj=sinu +R={} +nadgrids=@null +wktext".format(self.proj["a"]))
        
        #Locating upper left point in meters (on the projection)
        dummy_ind1 = self.structure.find("UpperLeftPointMtrs=(")
        dummy_ind2 = self.structure.find(")", dummy_ind1)
        [x_ul, y_ul] = numpy.array(self.structure[dummy_ind1+20:dummy_ind2].split(","), dtype="float")
        
        #Locating lower right point in meters (on the projection)
        dummy_ind1 = self.structure.find("LowerRightMtrs=(")
        dummy_ind2 = self.structure.find(")", dummy_ind1)
        [x_lr, y_lr] = numpy.array(self.structure[dummy_ind1+16:dummy_ind2].split(","), dtype="float")
        
        #Locating every other point
        x1 = numpy.linspace(x_ul, x_lr, self.npixels, endpoint=True)
        y1 = numpy.linspace(y_ul, y_lr, self.npixels, endpoint=True)
        x2,y2 = numpy.meshgrid(x1,y1)
        
        #Convertng x, y to lat/lon
        self.lons, self.lats = modis_grid(x2, y2, inverse=True)
                
        ### Creating attributes with data legends and matching colormaps
        ### Note that colormaps are my own suggestions and are not official in any capacity
        #IGBP legend (LC_Type1)
        self.lc1_legend = ["1 - Evergreen Needleleaf Forests", "2 - Evergreen Broadleaf Forests",
            "3 - Deciduous Needleleaf Forests", "4 - Deciduous Broadleaf Forests", 
            "5 - Mixed Forests", "6 - Closed Shrublands", "7 - Open Shrublands",
            "8 - Woody Savannas", "9 - Savannas", "10 - Grasslands",
            "11 - Permanent Wetlands", "12 - Croplands", "13 - Urban and Built-up",
            "14 - Cropland/Natural Vegetation Mosaics", "15 - Permanent Snow and Ice",
            "16 - Barren", "17 - Water Bodies", "255 - Unclassified"]
        
        #IGBP colormap (LC_Type1)        
        self.lc1_cmap = ListedColormap(["darkgreen", "forestgreen", "darkolivegreen",
            "olivedrab", "greenyellow", "olive", "darkkhaki", "darkorange", "orange",
            "yellow", "navy", "mediumorchid", "red", "mediumpurple", "lightcyan",
            "sienna", "blue", "grey"])
            
        #UMD legend (LC_Type2)
        self.lc2_legend = ["0 - Water bodies", "1 - Evergreen Needleleaf Forests",
            "2 - Evergreen Broadleaf Forests", "3 - Deciduous Needleleaf Forests", "4 - Deciduous Broadleaf Forests", 
            "5 - Mixed Forests", "6 - Closed Shrublands", "7 - Open Shrublands",
            "8 - Woody Savannas", "9 - Savannas", "10 - Grasslands",
            "11 - Permanent Wetlands", "12 - Croplands", "13 - Urban and Built-up",
            "14 - Cropland/Natural Vegetation Mosaics", "15 - Non-Vegetated Lands",
            "255 - Unclassified"]
            
        #LAI legend (LC_Type3)
        self.lc3_legend = ["0 - Water bodies", "1 - Grasslands", "2 - Shrublands",
            "3 - Broadleaf Croplands", "4 - Savannas", "5 - Evergreen Broadleaf Forests",
            "6 - Deciduous Broadleaf Forests", "7 - Evergreen Needleleaf Forests",
            "8 - Deciduous Needleleaf Forests", "9 - Non-vegetated Lands",
            "10 - Urban and Built-up Lands", "255 - Unclassified"]
            
        #BGC legend (LC_Type4)
        self.lc4_legend = ["0 - Water bodies", "1 - Evergreen Needleleaf Vegetation",
            "2 - Evergreen Broadleaf Vegetation", "3 - Deciduous Broadleaf Vegetation",
            "4 - Deciduous Needleleaf Vegetation", "5 - Annual Broadleaf Vegetation",
            "6 - Annual Grass Vegetation", "7 - Non-vegetated Lands",
            "8 - Urban and Built-up Lands", "255 - Unclassified"]
            
        #PFT legend (LC_Type5)
        self.lc5_legend = ["0 - Water bodies", "1 - Evergreen Needleleaf Trees",
            "2 - Evergreen Broadleaf Trees", "3 - Deciduous Broadleaf Trees",
            "4 - Deciduous Needleleaf Trees", "5 - Shrub",
            "6 - Grass", "7 - Cereal Croplands", "8 - Broadleaf Croplands"
            "9 - Urban and Built-up Lands", "10 - Permanent Snow and Ice",
            "11 - Barren", "255 - Unclassified"]
            
        #FAO-LCCS1 legend (LC_Prop1)
        self.lp1_legend = ["1 - Barren", "2 - Permanent Snow and Ice", "3 - Water bodies",
            "11 - Evergreen Needleleaf Forests", "12 - Evergreen Broadleaf Forests",
            "13 - Deciduous Needleleaf Forests", "14 - Deciduous Broadleaf Forests",
            "15 - Mixed Broadleaf/Needleleaf Forests",
            "16 - Mixed Evergreen Broadleaf/Needleleaf Forests", "21 - Open Forests",
            "22 - Sparse Forests", "31 - Dense Herbaceous", "32 - Sparse Herbaceous",
            "41 - Dense Shrublands", "42 - Shrubland/Grassland Mosaics",
            "43 - Sparse Shrublands", "255 - Unclassified"]
            
        #FAO-LCCS2 legend (LC_Prop2)
        self.lp2_legend = ["1 - Barren", "2 - Permanent Snow and Ice", "3 - Water Bodies",
            "9 - Urban and Built-up Lands", "10 - Dense Forests", "20 - Open Forests",
            "25 - Forest/Cropland Mosaics", "30 - Natural Herbaceous",
            "35 - natural Herbaceous/Croplands Mosaics", "36 - Herbaceous Croplands",
            "40 - Shrublands", "255 - Unclassified"]
            
        #FAO-LCCS3 legend (LC_Prop3)
        self.lp3_legend = ["1 - Barren", "2 - Permanent Snow and Ice", "3 - Water Bodies",
            "10 - Dense Forests", "20 - Open Forests", "27 - Woody Wetlands",
            "30 - Grasslands", "40 - Shrublands", "50 - Herbaceous Wetlands",
            "51 - Tundra", "255 - Unclassified"]
            
        #Quality Control legend (QC)
        self.qc_legend = ["0 - Classified Land", "1 - Unclassified Land",
            "2 - Classified Water", "3 - Unclassified Water",
            "4 - Classified Sea Ice", "5 - Misclassified water",
            "6 - Omitted Snow/Ice", "7 - Misclassified Snow/Ice",
            "8 - Backfilled Label", "9 - Forest Type Changed",
            "10 - No Data"]
        
    ### Method for accessing file variables
    def get(self, varname):
        return self.hfile.select(varname)[:]
        
    ### Method to clean up after object destruction
    def __del__(self):
        try:
            self.hfile.end()
        except:
            pass
        
    