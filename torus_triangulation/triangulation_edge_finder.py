from __future__ import print_function, absolute_import

from copy import deepcopy

from sage.all import *
from torus_triangulation.exceptions import VertexOrderError
from torus_triangulation.cube_triangulation import CubeTriangulation
from torus_triangulation.builder import TorusQuotientBuilder
from torus_triangulation.edge_accumulator import GroupBoxEdgeBuilder
from torus_triangulation.cube_filter import CubeFilter
from torus_triangulation.translate_point import TranslatePoint



class TriangulationEdgeFinder(object):
    
    def __init__(self, builder):
        self.builder = builder
        self.unit_filter = CubeFilter(*([[0, 1]] * self.builder.dimension))
        edge_graph = builder.edge_graph()
        vertices = self.unit_filter(edge_graph.vertices())
        self.unit_graph = DiGraph(edge_graph.subgraph(vertices))
        
    def unit_edges(self):
        """
        All edges that are going to be in the triangulation 
        """
        vertices = cartesian_product([[0,1]] * self.builder.dimension)
        triangulation = CubeTriangulation(*vertices)
        edges = set()
        for simplex in triangulation.simplex_iter():
            for i in range(len(simplex)):
                for j in range(i+1, len(simplex)):
                    edges.add((simplex[i], simplex[j]))
        return frozenset(edges)

    def triangulation_edges(self):
        """
        Currently missing edges compared to :meth:`unit_edges`
        """
        given = []
        for edge in self.unit_graph.edges():
            given.extend([(edge[0], edge[1]), (edge[1], edge[0])])
        return self.unit_edges().difference(given)

    def with_all_edges_to_origin(self):
        """
        Return new triangulation finder where all triangulation edges that can end at the origin do so.
        """
        builder = deepcopy(self.builder)
        origin = (0,) * self.builder.dimension
        for edge in self.triangulation_edges():
            if edge[0] == origin:
                builder.add_simplex_orbit(edge[1], origin)
            elif edge[1] == origin:
                builder.add_simplex_orbit(edge[0], origin)
            else:
                pass
        return TriangulationEdgeFinder(builder)
    
    def obvious_edge_iter(self):
        """
        Find triangulation edges that only fit in one orientation
        """
        for edge in self.triangulation_edges():
            b = deepcopy(self.builder)
            try:
                b.add_simplex_orbit(edge[0], edge[1])
                fit_01 = b.edge_graph().is_directed_acyclic()
            except VertexOrderError:
                fit_01 = False
            b = deepcopy(self.builder)
            try:
                b.add_simplex_orbit(edge[1], edge[0])
                fit_10 = b.edge_graph().is_directed_acyclic()
            except VertexOrderError:
                fit_10 = False
            if fit_01 and not fit_10:
                yield (edge[0], edge[1])
            elif fit_10 and not fit_01:
                yield (edge[1], edge[0])
            else:
                # direction not obvious
                pass
        
    def with_obvious_edges_added(self):
        builder = deepcopy(self.builder)
        for edge in self.obvious_edge_iter():
            builder.add_simplex_orbit(edge[0], edge[1])
        return TriangulationEdgeFinder(builder)

    def with_edge_added(self, v0, v1):
        builder = deepcopy(self.builder)
        builder.add_simplex_orbit(v0, v1)
        return TriangulationEdgeFinder(builder)
    
    def order_iter(self):
        for order in self.order_iter_for(self.unit_graph):
            yield order

    def order_iter_for(self, unit_graph):
        poset = Poset(unit_graph)
        translation = TranslatePoint(0, 1)
        for order1 in poset.linear_extensions():
            try:
                builder = deepcopy(self.builder)
                cube_triangulation = CubeTriangulation(*order1)
                for simplex in cube_triangulation:
                    builder.add_simplex_orbit(*simplex)
                order2 = map(translation, order1)
                cube_triangulation = CubeTriangulation(*order2)
                for simplex in cube_triangulation:
                    builder.add_simplex_orbit(*simplex)
                print('success', order1, order2)
                yield order1
            except VertexOrderError:
                print('failed')
    
    def edge_orientations_iter(self):
        """
        Iterate over all orientations of the remaining :meth:`triangulation_edges`
        """
        builder = deepcopy(self.builder)
        try:
            v0, v1 = next(iter(self.triangulation_edges()))
        except StopIteration:
            yield self   # terminate recursion
            return
        unit_graph = DiGraph(self.unit_graph)
        unit_graph.add_edge(v0, v1)
        if unit_graph.is_directed_acyclic():
            # print('trying', v0, v1)
            builder = deepcopy(self.builder)
            try:
                builder.add_simplex_orbit(v0, v1)
            except VertexOrderError:
                # print('vertex order')
                pass
            else:
                for recursion in TriangulationEdgeFinder(builder).edge_orientations_iter():
                    yield recursion
        unit_graph = DiGraph(self.unit_graph)
        unit_graph.add_edge(v1, v0)
        if unit_graph.is_directed_acyclic():
            # print('trying', v1, v0)
            builder = deepcopy(self.builder)
            try:
                builder.add_simplex_orbit(v1, v0)
            except VertexOrderError:
                # print('vertex order')
                pass
            else:
                for recursion in TriangulationEdgeFinder(builder).edge_orientations_iter():
                    yield recursion
        # print('done')
                    
    
