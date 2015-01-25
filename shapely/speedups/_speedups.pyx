# geos_linestring_from_py was transcribed from shapely.geometry.linestring
# geos_linearring_from_py was transcribed from shapely.geometry.polygon
# coordseq_ctypes was transcribed from shapely.coords.CoordinateSequence.ctypes
#
# Copyright (c) 2007, Sean C. Gillies
# Transcription to cython: Copyright (c) 2011, Oliver Tonnhofer

import ctypes

from shapely.coords import required, require_numeric
from shapely.geos import lgeos
from shapely.geometry import Point, LineString, LinearRing
from shapely.geometry.base import geom_factory

include "../_geos.pxi"
    

cdef inline GEOSGeometry *cast_geom(unsigned long geom_addr):
    return <GEOSGeometry *>geom_addr


cdef inline GEOSContextHandle_t cast_handle(unsigned long handle_addr):
    return <GEOSContextHandle_t>handle_addr


cdef inline GEOSCoordSequence *cast_seq(unsigned long handle_addr):
    return <GEOSCoordSequence *>handle_addr


def destroy(geom):
    GEOSGeom_destroy_r(cast_handle(lgeos.geos_handle), cast_geom(geom))


def fill_coordseq_array(self, seq, int n_points, int seq_offset, int array_offset):
    cdef GEOSContextHandle_t handle = cast_handle(lgeos.geos_handle)
    cdef GEOSCoordSequence *cs = <GEOSCoordSequence *><unsigned long>seq

    cdef double dx, dy, dz
    cdef double *cp
    cdef int i, sm, sn, seq_i, array_i
    cdef int n = self.shape[1]

    # Make pointer to the coordinate array
    if isinstance(self.array['data'], ctypes.Array):
        cp = <double *><unsigned long>ctypes.addressof(self.array['data'])
    else:
        cp = <double *><unsigned long>self.array['data'][0]

    # Use strides to properly index into cp
    # ob[i, j] == cp[sm*i + sn*j]
    # Just to avoid a referenced before assignment warning.
    dx = 0
    if self.array.get('strides', None):
        sm = self.array['strides'][0]/sizeof(dx)
        sn = self.array['strides'][1]/sizeof(dx)
    else:
        sm = n
        sn = 1

    for i in range(n_points):
        array_i = array_offset + i
        dx = cp[sm*array_i]
        dy = cp[sm*array_i+sn]
        dz = 0
        if n == 3:
            dz = cp[sm*array_i+2*sn]
            
        seq_i = seq_offset + i
        # Because of a bug in the GEOS C API, 
        # always set X before Y
        GEOSCoordSeq_setX_r(handle, cs, seq_i, dx)
        GEOSCoordSeq_setY_r(handle, cs, seq_i, dy)
        if n == 3:
            GEOSCoordSeq_setZ_r(handle, cs, seq_i, dz)


def fill_coordseq_object(self, seq, int n_points, int seq_offset, int array_offset):
    cdef GEOSContextHandle_t handle = cast_handle(lgeos.geos_handle)
    cdef GEOSCoordSequence *cs = <GEOSCoordSequence *><unsigned long>seq

    cdef double dx, dy, dz
    cdef int i, seq_i, array_i
    cdef int n = self.shape[1]

    for i in range(n_points):
        array_i = array_offset + i
        o = self.ob[array_i]

        if hasattr(o, 'coords'):
            coords = o.coords[0]
        else:
            coords = o

        dx = coords[0]
        dy = coords[1]
        dz = 0
        if n == 3:
            if len(coords) != 3:
                raise ValueError("Inconsistent coordinate dimensionality")
            dz = coords[2]

        seq_i = seq_offset + i
        # Because of a bug in the GEOS C API, 
        # always set X before Y
        GEOSCoordSeq_setX_r(handle, cs, seq_i, dx)
        GEOSCoordSeq_setY_r(handle, cs, seq_i, dy)
        if n == 3:
            GEOSCoordSeq_setZ_r(handle, cs, seq_i, dz)


def geos_linearring_from_py(ob, update_geom=None, update_ndim=0):
    cdef double *cp
    cdef GEOSContextHandle_t handle = cast_handle(lgeos.geos_handle)
    cdef GEOSGeometry *g
    cdef GEOSCoordSequence *cs
    cdef double dx, dy, dz
    cdef int i, n, m, M, sm, sn

    # If a LinearRing is passed in, just clone it and return
    # If a LineString is passed in, clone the coord seq and return a LinearRing
    if isinstance(ob, LineString):
        g = cast_geom(ob._geom)
        if GEOSHasZ_r(handle, g):
            n = 3
        else:
            n = 2

        if type(ob) == LinearRing:
            return <unsigned long>GEOSGeom_clone_r(handle, g), n
        else:
            cs = GEOSGeom_getCoordSeq_r(handle, g)
            GEOSCoordSeq_getSize_r(handle, cs, &m)
            if GEOSisClosed_r(handle, g) and m >= 4:
                cs = GEOSCoordSeq_clone_r(handle, cs)
                return <unsigned long>GEOSGeom_createLinearRing_r(handle, cs), n

    # If numpy is present, we use numpy.require to ensure that we have a
    # C-continguous array that owns its data. View data will be copied.
    ob = required(ob)
    try:
        # From array protocol
        array = ob.__array_interface__
        assert len(array['shape']) == 2
        m = array['shape'][0]
        n = array['shape'][1]
        if m < 3:
            raise ValueError(
                "A LinearRing must have at least 3 coordinate tuples")
        assert n == 2 or n == 3

        # Make pointer to the coordinate array
        if isinstance(array['data'], ctypes.Array):
            cp = <double *><unsigned long>ctypes.addressof(array['data'])
        else:
            cp = <double *><unsigned long>array['data'][0]

        # Use strides to properly index into cp
        # ob[i, j] == cp[sm*i + sn*j]
        dx = 0  # Just to avoid a referenced before assignment warning.
        if array.get('strides', None):
            sm = array['strides'][0]/sizeof(dx)
            sn = array['strides'][1]/sizeof(dx)
        else:
            sm = n
            sn = 1

        # Add closing coordinates to sequence?
        # Check whether the first set of coordinates matches the last.
        # If not, we'll have to close the ring later
        if (cp[0] != cp[sm*(m-1)] or cp[sn] != cp[sm*(m-1)+sn] or
            (n == 3 and cp[2*sn] != cp[sm*(m-1)+2*sn])):
            M = m + 1
        else:
            M = m

        # Create a coordinate sequence
        if update_geom is not None:
            cs = GEOSGeom_getCoordSeq_r(handle, cast_geom(update_geom))
            if n != update_ndim:
                raise ValueError(
                "Wrong coordinate dimensions; this geometry has dimensions: %d" \
                % update_ndim)
        else:
            cs = GEOSCoordSeq_create_r(handle, M, n)

        # add to coordinate sequence
        for i in xrange(m):
            dx = cp[sm*i]
            dy = cp[sm*i+sn]
            dz = 0
            if n == 3:
                dz = cp[sm*i+2*sn]
        
            # Because of a bug in the GEOS C API, 
            # always set X before Y
            GEOSCoordSeq_setX_r(handle, cs, i, dx)
            GEOSCoordSeq_setY_r(handle, cs, i, dy)
            if n == 3:
                GEOSCoordSeq_setZ_r(handle, cs, i, dz)

        # Add closing coordinates to sequence?
        if M > m:
            dx = cp[0]
            dy = cp[sn]
            dz = 0
            if n == 3:
                dz = cp[2*sn]
        
            # Because of a bug in the GEOS C API, 
            # always set X before Y
            GEOSCoordSeq_setX_r(handle, cs, M-1, dx)
            GEOSCoordSeq_setY_r(handle, cs, M-1, dy)
            if n == 3:
                GEOSCoordSeq_setZ_r(handle, cs, M-1, dz)
            
    except AttributeError:
        # Fall back on list
        try:
            m = len(ob)
        except TypeError:  # Iterators, e.g. Python 3 zip
            ob = list(ob)
            m = len(ob)
        n = len(ob[0])
        if m < 3:
            raise ValueError(
                "A LinearRing must have at least 3 coordinate tuples")
        assert (n == 2 or n == 3)

        # Add closing coordinates if not provided
        if m == 3 or ob[0][0] != ob[-1][0] or ob[0][1] != ob[-1][1]:
            M = m + 1
        else:
            M = m

        # Create a coordinate sequence
        if update_geom is not None:
            cs = GEOSGeom_getCoordSeq_r(handle, cast_geom(update_geom))
            if n != update_ndim:
                raise ValueError(
                "Wrong coordinate dimensions; this geometry has dimensions: %d" \
                % update_ndim)
        else:
            cs = GEOSCoordSeq_create_r(handle, M, n)
        
        # add to coordinate sequence
        for i in xrange(m):
            coords = ob[i]
            dx = coords[0]
            dy = coords[1]
            dz = 0
            if n == 3:
                dz = coords[2]
        
            # Because of a bug in the GEOS C API, 
            # always set X before Y
            GEOSCoordSeq_setX_r(handle, cs, i, dx)
            GEOSCoordSeq_setY_r(handle, cs, i, dy)
            if n == 3:
                GEOSCoordSeq_setZ_r(handle, cs, i, dz)

        # Add closing coordinates to sequence?
        if M > m:
            coords = ob[0]
            dx = coords[0]
            dy = coords[1]
            dz = 0
            if n == 3:
                dz = coords[2]
        
            # Because of a bug in the GEOS C API, 
            # always set X before Y
            GEOSCoordSeq_setX_r(handle, cs, M-1, dx)
            GEOSCoordSeq_setY_r(handle, cs, M-1, dy)
            if n == 3:
                GEOSCoordSeq_setZ_r(handle, cs, M-1, dz)

    if update_geom is not None:
        return None
    else:
        return <unsigned long>GEOSGeom_createLinearRing_r(handle, cs), n


def coordseq_ctypes(self):
    cdef int i, n, m
    cdef double temp = 0
    cdef GEOSContextHandle_t handle = cast_handle(lgeos.geos_handle)
    cdef GEOSCoordSequence *cs
    cdef double *data_p
    self._update()
    n = self._ndim
    m = self.__len__()
    array_type = ctypes.c_double * (m * n)
    data = array_type()
    
    cs = cast_seq(self._cseq)
    data_p = <double *><unsigned long>ctypes.addressof(data)
    
    for i in xrange(m):
        GEOSCoordSeq_getX_r(handle, cs, i, &temp)
        data_p[n*i] = temp
        GEOSCoordSeq_getY_r(handle, cs, i, &temp)
        data_p[n*i+1] = temp
        if n == 3: # TODO: use hasz
            GEOSCoordSeq_getZ_r(handle, cs, i, &temp)
            data_p[n*i+2] = temp
    return data

def coordseq_iter(self):
    cdef int i
    cdef double dx
    cdef double dy
    cdef double dz
    cdef int has_z

    self._update()

    cdef GEOSContextHandle_t handle = cast_handle(lgeos.geos_handle)
    cdef GEOSCoordSequence *cs
    cs = cast_seq(self._cseq)

    has_z = self._ndim == 3
    for i in range(self.__len__()):
        GEOSCoordSeq_getX_r(handle, cs, i, &dx)
        GEOSCoordSeq_getY_r(handle, cs, i, &dy)
        if has_z == 1:
            GEOSCoordSeq_getZ_r(handle, cs, i, &dz)
            yield (dx, dy, dz)
        else:
            yield (dx, dy)

cdef GEOSCoordSequence* transform(GEOSCoordSequence* cs,
                                  int ndim,
                                  double a,
                                  double b,
                                  double c,
                                  double d,
                                  double e,
                                  double f,
                                  double g,
                                  double h,
                                  double i,
                                  double xoff,
                                  double yoff,
                                  double zoff):
    """Performs an affine transformation on a GEOSCoordSequence
    
    Returns the transformed coordinate sequence
    """
    cdef GEOSContextHandle_t handle = cast_handle(lgeos.geos_handle)
    cdef int m
    cdef GEOSCoordSequence *cs_t
    cdef double x, y, z
    cdef double x_t, y_t, z_t
    
    # create a new coordinate sequence with the same size and dimensions
    GEOSCoordSeq_getSize_r(handle, cs, &m)
    cs_t = GEOSCoordSeq_create_r(handle, m, ndim)
    
    # perform the transform
    for n in range(0, m):
        GEOSCoordSeq_getX_r(handle, cs, n, &x)
        GEOSCoordSeq_getY_r(handle, cs, n, &y)
        x_t = a * x + b * y + xoff
        y_t = d * x + e * y + yoff
        GEOSCoordSeq_setX_r(handle, cs_t, n, x_t)
        GEOSCoordSeq_setY_r(handle, cs_t, n, y_t)
    if ndim == 3:
        for n in range(0, m):
            GEOSCoordSeq_getZ_r(handle, cs, n, &z)
            z_t = g * x + h * y + i * z + zoff
            GEOSCoordSeq_setZ_r(handle, cs_t, n, z_t)
    
    return cs_t

cpdef affine_transform(geom, matrix):
    cdef double a, b, c, d, e, f, g, h, i, xoff, yoff, zoff
    if geom.is_empty:
        return geom
    if len(matrix) == 6:
        ndim = 2
        a, b, d, e, xoff, yoff = matrix
        if geom.has_z:
            ndim = 3
            i = 1.0
            c = f = g = h = zoff = 0.0
            matrix = a, b, c, d, e, f, g, h, i, xoff, yoff, zoff
    elif len(matrix) == 12:
        ndim = 3
        a, b, c, d, e, f, g, h, i, xoff, yoff, zoff = matrix
        if not geom.has_z:
            ndim = 2
            matrix = a, b, d, e, xoff, yoff
    else:
        raise ValueError("'matrix' expects either 6 or 12 coefficients")
    
    cdef GEOSContextHandle_t handle = cast_handle(lgeos.geos_handle)
    cdef GEOSCoordSequence *cs
    cdef GEOSCoordSequence *cs_t
    cdef GEOSGeometry *the_geom
    cdef GEOSGeometry *the_geom_t
    cdef int m, n
    cdef double x, y, z
    cdef double x_t, y_t, z_t
    
    # Process coordinates from each supported geometry type
    if geom.type in ('Point', 'LineString', 'LinearRing'):
        the_geom = cast_geom(geom._geom)
        cs = GEOSGeom_getCoordSeq_r(handle, the_geom)
        
        # perform the transformation
        cs_t = transform(cs, ndim, a, b, c, d, e, f, g, h, i, xoff, yoff, zoff)
        
        # create a new geometry from the coordinate sequence
        if geom.type == 'Point':
            the_geom_t = GEOSGeom_createPoint_r(handle, cs_t)
        elif geom.type == 'LineString':
            the_geom_t = GEOSGeom_createLineString_r(handle, cs_t)
        elif geom.type == 'LinearRing':
            the_geom_t = GEOSGeom_createLinearRing_r(handle, cs_t)
        
        # return the geometry as a Python object
        return geom_factory(<unsigned long>the_geom_t)
    elif geom.type == 'Polygon':
        ring = geom.exterior
        shell = affine_transform(ring, matrix)
        holes = list(geom.interiors)
        for pos, ring in enumerate(holes):
            holes[pos] = affine_transform(ring, matrix)
        return type(geom)(shell, holes)
    elif geom.type.startswith('Multi') or geom.type == 'GeometryCollection':
        return type(geom)([affine_transform(part, matrix) for part in geom.geoms])
    else:
        raise ValueError('Type %r not recognized' % geom.type)
