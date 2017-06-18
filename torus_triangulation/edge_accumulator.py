from __future__ import print_function, absolute_import

from collections import defaultdict
from sage.all import *
from sage.misc.cachefunc import cached_method


from torus_triangulation.remove_one import remove_one_iter
from torus_triangulation.translate_point import TranslatePoint
from torus_triangulation.exceptions import VertexOrderError


    
class EdgeBuilder(object):

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

    def test(self):
        cells = self.cells
        cells.test()
        assert self.simplices[self.dimension + 1]
        assert cells.dimension == self.dimension

    def edge_plot(self):
        return sum(
            arrow3d(*edge) for edge in self.simplices[2]
        ) + sum(
            Polyhedron(tetrahedron).plot() for tetrahedron in self.simplices[4]
        )

    def edge_graph(self):
        dg = DiGraph()
        for vertex in self.simplices[1]:
            dg.add_vertex(*vertex)
        for edge in self.simplices[2]:
            dg.add_edge(*edge)
        return dg

    
class BoxEdgeBuilder(EdgeBuilder):

    def __init__(self, *side_length):
        self.side_length = tuple(side_length)
        super(BoxEdgeBuilder, self).__init__(len(self.side_length))
        ranges = []
        for r in self.side_length:
            ranges.append(range(r + 1))
        self.points = tuple(cartesian_product(ranges))

    def one_simplices(self):
        ranges = []
        for r in self.side_length:
            ranges.append(range(r))
        result = []
        for vertex in cartesian_product(ranges):
            for i in range(len(self.side_length)):
                tip = list(vertex)
                tip[i] += 1
                result.append((tuple(vertex), tuple(tip)))
        return result
        
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
        super(BoxEdgeBuilder, self).test()
        assert len(self.side_length) == self.dimension


class GroupBoxEdgeBuilder(BoxEdgeBuilder):

    def __init__(self, *side_length, **kwds):
        super(GroupBoxEdgeBuilder, self).__init__(*side_length)
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

    def _make_torus_images(self, *vertices):
        images = [self._make_simplex(*vertices)]
        for i, r in enumerate(self.side_length):
            if all(v[i] == 0 for v in vertices):
                t = TranslatePoint(i, r)
                images.extend([self._make_simplex(*map(t, s)) for s in images])
        return images
    
    def orbit(self, *vertices):
        """
        Return the orbit of simplex
        """
        simplex = self._make_fundamental_region_simplex(*vertices)
        orbit = set(self._make_torus_images(*simplex))
        todo = set([simplex])
        while todo:
            simplex = todo.pop()
            for gen in self.group_gens:
                img = self._make_fundamental_region_simplex(*map(gen, simplex))
                if img not in orbit:
                    orbit.add(img)
                    orbit.update(self._make_torus_images(*img))
                    todo.add(img)
        return tuple(orbit)

    def add_simplex_orbit(self, *vertices):
        """
        Add simplex and its images under the group action
        """
        orbit = self.orbit(*vertices)
        for simplex in orbit:
            self.add_simplex(*simplex)

    def __add__(self, other):
        assert self.dimension == other.dimension
        assert self.side_length == other.side_length
        assert self.group_gens == other.group_gens
        total = self.__class__(*self.side_length, group_gens=self.group_gens)
        for dim in range(0, self.dimension + 1):
            for simplex in self.simplices[dim]:
                total.add_simplex(*simplex)
            for simplex in other.simplices[dim]:
                total.add_simplex(*simplex)
        return total
    
    def reversed(self):
        rev = self.__class__(*self.side_length, group_gens=self.group_gens)
        for dim, simplices in self.simplices.items():
            for simplex in simplices:
                rev.add_simplex(*reversed(simplex))
        return rev
    
    @classmethod
    def edge_builders(cls, *args, **kwds):
        builder = cls(*args, **kwds)
        edges = set(tuple(sorted(edge)) for edge in builder.one_simplices())
        result = []
        while edges:
            edge = edges.pop()
            builder = cls(*args, **kwds)
            builder.add_simplex_orbit(*edge)
            result.append(builder)
            orbit = builder.orbit(*edge)
            edges.difference_update(tuple(sorted(edge)) for edge in orbit)
        return result

    @classmethod
    def acyclic_edge_iter(cls, *args, **kwds):
        builders = cls.edge_builders(*args, **kwds)
        builder_and_reversed = [(b, b.reversed()) for b in builders]
        for b_test in cartesian_product(builder_and_reversed):
            b_sum = b_test[0]
            for b in b_test[1:]:
                b_sum = b_sum + b
            if b_sum.edge_graph().is_directed_acyclic():
                # print('-' * 79)
                # print(b_sum, [b_i is builders_i for b_i, builders_i in zip(b_test, builders)])
                # print(map(len, Poset(b_sum.edge_graph()).level_sets()))
                # break
                yield b_sum
    
    def ordered(self):
        acyclic = self.edge_graph().is_directed_acyclic(certificate=True)
        assert acyclic[0]
        return acyclic[1]

    def level_sets(self):
        return Poset(self.edge_graph()).level_sets()


