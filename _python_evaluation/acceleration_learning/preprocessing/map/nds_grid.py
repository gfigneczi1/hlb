S32_MIN = int(-2147483647 - 1)
S32_MAX = int(2147483647)
ATL = pow(2, 32) / 360.0
from shapely.geometry import Polygon
import os
os.environ['USE_PYGEOS'] = '0'
import geopandas as gpd

def create_nds_grid(north: float, south: float, east: float, west: float, ndsLevel: int=13):
    pow2k = pow(2, ndsLevel)
    dx = 180.0 / pow2k
    dy = 90.0 / (pow2k / 2.0)
    ulTileSizeInLev0MAPU = int(S32_MAX + 1)
    sizeOfTilesInLevelMAPU = int(ulTileSizeInLev0MAPU >> ndsLevel)

    sxmin = _s32(west)
    symin = _s32(south)
    sxmax = _s32(east)
    symax = _s32(north)

    minCol = _getColumn(sxmin, sizeOfTilesInLevelMAPU) - 1
    maxCol = _getColumn(sxmax, sizeOfTilesInLevelMAPU) + 2
    minRow = _getRow(symin, sizeOfTilesInLevelMAPU) - 1
    maxRow = _getRow(symax, sizeOfTilesInLevelMAPU) + 2

    geoms = []
    ids = []
    for col in range(minCol, maxCol):
        for row in range(minRow, maxRow):
            geo, id = _createPolygon(row, col, dx, dy, ndsLevel)
            geoms.append(geo)
            ids.append(id)

    return gpd.GeoDataFrame({'NDS Tile ID':ids, 'geometry':geoms})

def _createPolygon(rownum, colnum, dx, dy, ndsLevel):
    x = dx * colnum
    y = dy * rownum

    tileId = _tileId(rownum, colnum, ndsLevel)

    polygonPoints = [(x, y), (x, y+dy), (x+dx, y+dy), (x+dx, y)]

    return Polygon(polygonPoints), tileId

def get_tileids_for_extent(north: float, south: float, east: float, west: float, nds_level: int=13):
    ulTileSizeInLev0MAPU = int(S32_MAX + 1)
    sizeOfTilesInLevelMAPU = int(ulTileSizeInLev0MAPU >> nds_level)
    sxmin = _s32(west)
    symin = _s32(south)
    sxmax = _s32(east)
    symax = _s32(north)
    minCol = _getColumn(sxmin, sizeOfTilesInLevelMAPU)
    maxCol = _getColumn(sxmax, sizeOfTilesInLevelMAPU)
    minRow = _getRow(symin, sizeOfTilesInLevelMAPU)
    maxRow = _getRow(symax, sizeOfTilesInLevelMAPU)
    # create all cells of the grid
    # print(minCol, maxCol, minRow, maxRow)
    tiles=set()

    col = minCol

    while col <= maxCol:
        row = minRow
        while row <= maxRow:
            tileId = _tileId(row, col, nds_level)
            tiles.add(tileId)
            row +=1
        col += 1

    return list(tiles)

def _getRow(yInt, sizeOfTilesInLevelMAPU):
    if 0 > yInt:
        # minus 1 because the first negative row is -1 and not -0
        # and because all digits after the decimal point are truncated
        return int((yInt / sizeOfTilesInLevelMAPU) - 1)
    return int(yInt / sizeOfTilesInLevelMAPU)

def _getColumn(xInt, sizeOfTilesInLevelMAPU):
    if 0 > xInt:
        # minus 1 because the first negative column is -1 and not -0
        # and because all digits after the decimal point are truncated
        return int((xInt / sizeOfTilesInLevelMAPU) - 1)
    return int(xInt / sizeOfTilesInLevelMAPU)

def _s32(flDegree):
    if +180.0 == flDegree:
        return S32_MAX
    if -180.0 == flDegree:
        return S32_MIN
    return int(flDegree * ATL)

def _tileId(slRow, slColumn, ndsLevel):
    # Loop over the level+1 last Bits found in slRow and slColumn.
    # Always copy the row-Bit to the same location in tileNumber.
    # Then shift tileNumber one to the Left and copy the column-Bit.
    tileNumber = 0
    mask = (1 << ndsLevel)

    # Level 0 gets only one bit.
    # It is the inverse of the bit found in slColumn[level]
    tileNumber |= (slColumn & mask)
    mask >>= 1

    while mask > 0:
        # The remaining Levels get the Bits as they are found in slRow and slColumn
        tileNumber |= (slRow & mask)
        tileNumber <<= 1
        tileNumber |= (slColumn & mask)
        mask >>= 1

    # Now all necessary Bits of the TileNumber are set.
    # The Tile-Id can get derived from Level and TileNumber :
    compactId = 1
    compactId <<= (16 + ndsLevel)
    return int(compactId | tileNumber)