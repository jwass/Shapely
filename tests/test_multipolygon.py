from . import unittest, numpy
from shapely.geometry import (
    Point, LineString, Polygon, MultiPolygon, asMultiPolygon)
from shapely.geometry.base import dump_coords


class MultiPolygonTestCase(unittest.TestCase):

    def test_multipolygon(self):

        # From coordinate tuples
        geom = MultiPolygon(
            [(((0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)),
              [((0.25, 0.25), (0.25, 0.5), (0.5, 0.5), (0.5, 0.25))])])
        self.assertIsInstance(geom, MultiPolygon)
        self.assertEqual(len(geom.geoms), 1)
        self.assertEqual(
            dump_coords(geom),
            [[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0),
              [(0.25, 0.25), (0.25, 0.5), (0.5, 0.5), (0.5, 0.25),
               (0.25, 0.25)]]])

        # Or from polygons
        p = Polygon(((0, 0), (0, 1), (1, 1), (1, 0)),
                    [((0.25, 0.25), (0.25, 0.5), (0.5, 0.5), (0.5, 0.25))])
        geom = MultiPolygon([p])
        self.assertEqual(len(geom.geoms), 1)
        self.assertEqual(
            dump_coords(geom),
            [[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0),
              [(0.25, 0.25), (0.25, 0.5), (0.5, 0.5), (0.5, 0.25),
               (0.25, 0.25)]]])

        # Or from another multi-polygon
        geom2 = MultiPolygon(geom)
        self.assertEqual(len(geom2.geoms), 1)
        self.assertEqual(
            dump_coords(geom2),
            [[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0),
              [(0.25, 0.25), (0.25, 0.5), (0.5, 0.5), (0.5, 0.25),
               (0.25, 0.25)]]])

        # Sub-geometry Access
        self.assertIsInstance(geom.geoms[0], Polygon)
        self.assertEqual(
            dump_coords(geom.geoms[0]),
            [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0),
             [(0.25, 0.25), (0.25, 0.5), (0.5, 0.5), (0.5, 0.25),
              (0.25, 0.25)]])
        with self.assertRaises(IndexError):  # index out of range
            geom.geoms[1]

        # Geo interface
        self.assertEqual(
            geom.__geo_interface__,
            {'type': 'MultiPolygon',
             'coordinates': [(((0.0, 0.0), (0.0, 1.0), (1.0, 1.0),
                               (1.0, 0.0), (0.0, 0.0)),
                              ((0.25, 0.25), (0.25, 0.5), (0.5, 0.5),
                               (0.5, 0.25), (0.25, 0.25)))]})

        # Adapter
        coords = ((0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0))
        holes_coords = [((0.25, 0.25), (0.25, 0.5), (0.5, 0.5), (0.5, 0.25))]
        mpa = asMultiPolygon([(coords, holes_coords)])
        self.assertEqual(len(mpa.geoms), 1)
        self.assertEqual(len(mpa.geoms[0].exterior.coords), 5)
        self.assertEqual(len(mpa.geoms[0].interiors), 1)
        self.assertEqual(len(mpa.geoms[0].interiors[0].coords), 5)


@unittest.skipIf(not numpy, 'Numpy required')
class NumpyMultiPolygonStringTestCase(unittest.TestCase):
    def setUp(self):
        self.ex1 = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        self.in1 = [(0.25, 0.25), (0.5, 0.25), (0.5, 0.5), (0.25, 0.5)]
        self.in2 = [(0.6, 0.25), (0.9, 0.25), (0.9, 0.5), (0.6, 0.6)]

        self.ex2 = [(5.0, 0.0), (6.0, 0.0), (6.0, 1.0), (6.0, 1.0)]
        self.in3 = [(5.25, 0.25), (5.5, 0.25), (5.5, 0.5), (5.25, 0.5)]
        self.in4 = [(5.6, 0.25), (5.9, 0.25), (5.9, 0.5), (5.6, 0.6)]

        self.mp = MultiPolygon([(self.ex1, [self.in1, self.in2]),
                                (self.ex2, [self.in3, self.in4])])

    def test_from_array(self):
        arr = numpy.array
        mp = MultiPolygon([(arr(self.ex1), [arr(self.in1), arr(self.in2)]),
                           (arr(self.ex2), [arr(self.in3), arr(self.in4)])])

        self.assertEqual(mp, self.mp)

    def test_from_polygon_array(self):
        p1 = Polygon(self.ex1, [self.in1, self.in2])
        p2 = Polygon(self.ex2, [self.in3, self.in4])

        mp = MultiPolygon(numpy.array([p1, p2], dtype='object'))

        self.assertEqual(mp, self.mp)

    def test_from_array_mix(self):
        arr = numpy.array

        in2 = [self.in2[0], Point(self.in2[1]), self.in2[2], Point(self.in2[3])]
        # p1 defined as (list of coords, [array([coords]), mixed array])
        p1 = (self.ex1, [arr(self.in1), arr(in2, dtype='object')])

        in4 = arr([Point(c) for c in self.in4], dtype='object')
        # p2 defined as (linestring, [list of coords, mixed array])
        p2 = (LineString(self.ex2), [self.in3, in4])

        mp = MultiPolygon([p1, p2])
        self.assertEqual(mp, self.mp)


def test_suite():
    return unittest.TestLoader().loadTestsFromTestCase(MultiPolygonTestCase)
