from __future__ import print_function
from collections import defaultdict
from sage.all import *
from sage.misc.cachefunc import cached_method

from torus_triangulation.remove_one import remove_one_iter
from torus_triangulation.translate_point import TranslatePoint
from torus_triangulation.exceptions import VertexOrderError
from torus_triangulation.cell_complex import CellCollection, SimplexRelation


class Builder(object):

    def __init__(self, dimension):
        self.dimension = dimension
        self._sorted_simplices = dict()
        self.simplices = defaultdict(set)
        self.boundary = dict()

    def _make_simplex(self, *vertices):
        simplex = tuple(map(tuple, vertices))
        sorted_simplex = tuple(sorted(simplex))
        try:
            previous = self._sorted_simplices[sorted_simplex]
        except KeyError:
            self._sorted_simplices[sorted_simplex] = simplex
        else:
            if previous != simplex:
                raise VertexOrderError(
                    'Previously-constructed simplex has conflicting vertex order: {0} != {1}'
                    .format(previous, simplex))
        return tuple(map(tuple, vertices))
    
    @property
    def cells(self):
        return CellCollection(self.dimension, self.simplices, self.boundary)

    def delta_complex(self):
        return self.cells.delta_complex()

    def test(self):
        cells = self.cells
        cells.test()
        assert self.simplices[self.dimension + 1]
        assert cells.dimension == self.dimension


class BoxBuilder(Builder):

    def __init__(self, *side_length):
        self.side_length = tuple(side_length)
        super(BoxBuilder, self).__init__(len(self.side_length))
        ranges = []
        for r in self.side_length:
            ranges.append(range(r + 1))
        self.points = tuple(cartesian_product(ranges))

    def add_simplex(self, *vertices):
        """
        Add the simplex defined by the vertices

        Returns:
            Normal form of the simplex
        """
        simplex = self._make_simplex(*vertices)
        assert all(v in self.points for v in simplex)
        dim = len(simplex)
        d_simplices = self.simplices[dim]
        d_simplices.add(simplex)
        simplex_boundary = []
        for face in remove_one_iter(simplex):
            simplex_boundary.append(self.add_simplex(*face))
        self.boundary[simplex] = tuple(simplex_boundary)
        return simplex

    def facets_at_boundary(self, coordinate_index, inner=True):
        value = 0 if inner else self.side_length[coordinate_index]
        return tuple(facet for facet in self.simplices[self.dimension]
                     if all(p[coordinate_index] == value for p in facet))

    def test(self):
        super(BoxBuilder, self).test()
        assert len(self.side_length) == self.dimension
        HH = self.delta_complex().homology()
        for dim in range(self.dimension + 1):
            assert HH[dim].ngens() == 0, 'Homology of box should vanish, got {0}'.format(HH)

    
class TorusBuilder(BoxBuilder):

    def torus_relation(self):
        relation = SimplexRelation(self)
        for index in range(self.dimension):
            facets0 = self.facets_at_boundary(index, True)
            translation = TranslatePoint(index, self.side_length[index])
            for facet0 in facets0:
                facet1 = self._make_simplex(*map(translation, facet0))
                assert facet1 in self.simplices[len(facet0)]
                relation.identify(facet1, facet0)
        return relation
    
    @property
    def torus_cells(self):
        relation = self.torus_relation()
        return self.cells.quotient(relation)

    def test(self):
        super(TorusBuilder, self).test()
        self.torus_relation().test()
        torus_cells = self.torus_cells
        torus_cells.test()
        HH = torus_cells.delta_complex().homology(reduced=False)
        for dim in range(1, self.dimension + 1):
            expected_rank = binomial(self.dimension, dim)
            expected_invariants = (0,) * expected_rank
            assert HH[dim].invariants() == expected_invariants, \
                'Homology of torus should have rank {0}, got {1}'.format(expected_rank, HH[dim])

    

class TorusQuotientBuilder(TorusBuilder):

    def __init__(self, *circumference, **kwds):
        super(TorusQuotientBuilder, self).__init__(*circumference)
        self.group_gens = tuple(kwds['group_gens'])

    def _make_fundamental_region_simplex(self, *vertices):
        """
        Translate simplex given by vertices to fundamental region (box)
        """
        if len(vertices) == 0:
            return self._make_simplex()
        for i, r in enumerate(self.side_length):
            d = min(v[i] for v in vertices) // r
            t = TranslatePoint(i, -d * r)
            if t:
                vertices = map(t, vertices)
        return self._make_simplex(*vertices)
    
    def orbit(self, *vertices):
        """
        Return the orbit of simplex
        """
        simplex = self._make_fundamental_region_simplex(*vertices)
        orbit = set([simplex])
        todo = set([simplex])
        while todo:
            simplex = todo.pop()
            for gen in self.group_gens:
                img = self._make_fundamental_region_simplex(*map(gen, simplex))
                if img not in orbit:
                    orbit.add(img)
                    todo.add(img)
        return tuple(orbit)

    def top_orbit_iter(self):
        """
        Iterate over the orbits of the top-dimensional cells
        """
        simplices = set(self.simplices[self.dimension + 1])
        while simplices:
            simplex = next(iter(simplices))
            orbit = self.orbit(*simplex)
            # print(simplex, orbit)
            for img in orbit:
                assert img in simplices, 'invalid group action'
            yield orbit
            simplices.difference_update(orbit)

    def add_simplex_orbit(self, *vertices):
        """
        Add simplex and its images under the group action
        """
        orbit = self.orbit(*vertices)
        for simplex in orbit:
            self.add_simplex(*simplex)
    
    def quotient_relation(self):
        relation = self.torus_relation()
        for orbit in self.top_orbit_iter():
            relation.identify(*orbit)
        return relation

    @property
    def quotient_cells(self):
        relation = self.quotient_relation()
        return self.cells.quotient(relation)

    def test(self):
        super(TorusQuotientBuilder, self).test()
        self.quotient_relation().test()

