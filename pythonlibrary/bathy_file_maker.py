import os
import requests
import numpy as np
from osgeo import gdal
from osgeo import osr


class BathyFileMaker(object):
    """Used to create a bathymetry file for SI3D."""
    def __init__(self, name='', dem=None, shoreline_shp=None, wse=None, out_dir=os.getcwd(), **kwargs):
        """
        :param name: Name for the domain bathymetry, copied to SI3D bathy file header.
        :type name: str
        :param dem: Path to input DEM/bathymetry raster. Must be in a projected GCS with horizontal and vertical units
                    of meters. If values are absolute elevations (not depths), need to provide "wse" argument as well.
        :type dem: str
        :param shoreline_shp: (optional) Path to input shoreline polygon shapefile.
                              DEM will be clipped within this polygon area to set model domain.
        :type shoreline_shp: str
        :param wse: (optional) Elevation of the water surface (used if the input DEM values are elevations). If given,
                    depths calculated from DEM elevations as: depth = WSE - DEM. Otherwise, assumed that DEM values
                    are depths.
        :type wse: float
        :param out_dir: output directory, default = working directory
        :type out_dir: str
        :param kwargs:

        TODO:
            - output processed bathy array as georaster for reference
        """
        self.name = name
        self.dem = dem
        self.shoreline_shp = shoreline_shp
        self.wse = wse
        self.out_dir = out_dir
        self.kwargs = kwargs
        # use name of DEM raster if no name is given
        self.dem_name = os.path.basename(dem).split('.')[0]
        if not len(self.name):
            self.name = self.dem_name
        # get projection, horizontal units and cell size of DEM
        self.proj, self.h_unit, self.cell_size = self.get_projection()
        # get DEM raster as array
        self.dem_array = self.get_dem_array()
        self.num_rows, self.num_cols = np.shape(self.dem_array)
        # generate bathy file
        self.make_bathy_file()

    @property
    def dem(self):
        return self._dem

    @staticmethod
    def valid_input(dem):
        valid_path = os.path.exists(dem)
        valid_url = False
        try:
            valid_url = requests.head(dem).status_code < 400
        except:
            pass
        return valid_path or valid_url

    @dem.setter
    def dem(self, dem):
        # make sure the given DEM exists
        if not self.valid_input(dem):
            raise FileNotFoundError(f"Cannot find input DEM: {dem}")
        self._dem = dem
        # if given DEM is in ASCII, convert to GeoTIFF
        if dem.lower().endswith('.asc'):
            print("Given DEM was in ASCII format, converting to GeoTIFF.")
            self._dem = self._asc_to_tif(dem)
        # if given DEM doesn't have NoData value, assume it is zero
        self._check_dem_nodata()
        return

    @property
    def out_dir(self):
        return self._out_dir

    @out_dir.setter
    def out_dir(self, out_dir):
        if not os.path.exists(out_dir):
            raise IOError(f"The output directory does not exist: {out_dir}")
        self._out_dir = out_dir
        return

    def _asc_to_tif(self, asc):
        """Convert ascii grid to geotiff"""
        tif_name = os.path.basename(asc).split('.')[0] + ".tif"
        gdal.Translate(f"{tif_name}", f"{asc}", format="GTiff")
        return tif_name

    def _check_dem_nodata(self):
        """Check that input DEM has a NoData value set. If not, set to zero."""
        r = gdal.Open(self._dem, gdal.GA_Update)
        band = r.GetRasterBand(1)
        if band.GetNoDataValue() is None:
            raise IOError("WARNING: NoData value not found for input DEM.")
        return

    def get_projection(self):
        """Get EPSG code for DEM raster projection."""
        ras = gdal.Open(self.dem, gdal.GA_ReadOnly)
        proj = osr.SpatialReference(wkt=ras.GetProjection())
        epsg_code = proj.GetAttrValue('AUTHORITY', 1)
        epsg_code = "EPSG:" + epsg_code if epsg_code is not None else None
        h_unit = proj.GetAttrValue('UNIT')
        if epsg_code is None or h_unit is None:
            raise IOError("ERROR: CRS metadata is missing for input DEM.")
        if h_unit == "degree":
            raise IOError("ERROR: Input DEM must be in a projected CRS with units of meters.")
        if h_unit != "metre":
            raise IOError("ERROR: Input DEM must have units of meters.")
        upper_left_x, x_size, x_rot, upper_left_y, y_rot, y_size = ras.GetGeoTransform()
        if abs(x_size) != abs(y_size):
            raise IOError("ERROR: Raster cells must be square (i.e. cell width = cell height) for SI3D input.")
        return epsg_code, h_unit, x_size

    def get_dem_array(self):
        """Get DEM as array"""
        if self.shoreline_shp:
            print('Cropping DEM to shoreline polygon...')
            gdal.Warp(f"{self.name}_crop.tif", self.dem, cutlineDSName=self.shoreline_shp)
            dem = f"{self.name}_crop.tif"
        else:
            dem = self.dem
        ras = gdal.Open(dem, gdal.GA_ReadOnly)
        arr = ras.GetRasterBand(1).ReadAsArray()
        # convert elevation values to depth values (if applicable)
        if self.wse is not None:
            print("Converting elevation values to depths...")
            arr = self.wse - arr
        # set NoData values to -99
        print("Setting NoData values to -99")
        nodata_val = ras.GetRasterBand(1).GetNoDataValue()
        arr = np.where(np.isnan(arr) | (arr == nodata_val), -99, arr)
        # set shallow areas with depths <= 0.5 dm (i.e. not rounding to at least 1 dm depth) to nodata/-99
        arr = np.where(arr <= 0.05, -99, arr)
        # depth must be int < 5 digits
        if np.max(arr) > 999.9:
            raise IOError(f"ERORR: Input DEM depth of {np.max(arr):.1f} m exceeds max value of 999.9 meters.")
        # convert m to dm
        print("Converting depths to dm...")
        arr = np.where(arr != -99, arr * 10, -99)
        return arr

    def get_header(self):
        header = f"{self.name} (dx= {self.cell_size}m),   " \
                 f"imx =  {self.num_cols:d},jmx =  {self.num_rows:d},ncols = {self.num_cols:d}"
        return header

    def make_bathy_file(self):
        print('Writing SI3D bathy file...')
        header = self.get_header()
        filename = f"h{self.cell_size:.0f}m_Lake"
        with open(os.path.join(self.out_dir, filename), 'w+') as f:
            lines = []
            # write descriptive header
            lines.append(header)
            # write H/V line for tracking row/cols
            lines.append("HV    " + ("   V" * self.num_cols))
            # write line numbering each col
            lines.append("     " + "".join([f"{i + 2:d}".rjust(5) for i in range(self.num_cols)]))
            for j in range(self.num_rows):
                dem_row = self.dem_array[j]
                # number rows decreasing
                row_num = f"{self.num_rows - j + 1:d}".rjust(5)
                # 4.0f float precision for bathy values
                dem_row_str = "".join([f"{val:4.0f}".rjust(5) for val in dem_row])
                lines.append(row_num + dem_row_str)
            f.writelines([line + '\n' for line in lines])
        print(f'Saved SI3D bathy file: {filename}')
        return filename


if __name__ == "__main__":
    dem = 'tahoe_bathy_400m.tif'
    shoreline_shp = 'tahoe_shoreline.shp'
    #shoreline_shp = None
    bfm = BathyFileMaker(name='Tahoe', dem=dem, shoreline_shp=shoreline_shp)
