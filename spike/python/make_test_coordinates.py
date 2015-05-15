from ocgis import SpatialDimension, CoordinateReferenceSystem, SpatialGridDimension, VectorDimension, RequestDataset
import itertools


def get_vector():
    rowv = [39.4, 39.8]
    colv = [-105.4, -105.2]
    proj_src = "+proj=longlat +a=6370997 +b=6370997 +towgs84=0,0,0,0,0,0,0 +no_defs "
    proj_dst = "+proj=lcc +lat_1=30 +lat_2=60 +lat_0=47.5 +lon_0=-97 +x_0=3325000 +y_0=2700000 +ellps=WGS84 +units=m +no_defs "

    row = VectorDimension(value=rowv)
    col = VectorDimension(value=colv)
    grid = SpatialGridDimension(row=row, col=col)
    crs = CoordinateReferenceSystem(proj4=proj_src)
    sdim = SpatialDimension(grid=grid, crs=crs)
    sdim.write_fiona('/tmp/original.shp', target='point')
    print '-- original --'
    print 'row: {0}'.format(sdim.grid.value[0].flatten().tolist())
    print 'col: {0}'.format(sdim.grid.value[1].flatten().tolist())
    print ''

    crs_dst = CoordinateReferenceSystem(proj4=proj_dst)
    sdim.write_fiona('/tmp/transformed.shp', target='point')
    sdim.update_crs(crs_dst)
    print '-- transformed --'
    print 'row: {0}'.format(sdim.grid.value[0].flatten().tolist())
    print 'col: {0}'.format(sdim.grid.value[1].flatten().tolist())


def get_large_vector():
    uri = '/home/ben.koziol/ocgis_test_data/nc/CanCM4/tas_day_CanCM4_decadal2000_r2i1p1_20010101-20101231.nc'
    rd = RequestDataset(uri=uri)
    field = rd.get()
    print field.shape
    with open('/tmp/my_coords.txt', 'w') as f:
        f.write(str(field.spatial.grid.value[1].flatten().tolist()))
        f.write('/n')
        f.write(str(field.spatial.grid.value[0].flatten().tolist()))


def get_edge_effects():
    x = [0, 180, 360]
    y = [-90, 0, 90]
    for c in itertools.product(x, y):
        print c


if __name__ == '__main__':
    # get_vector()
    # get_large_vector()
    get_edge_effects()