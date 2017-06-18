from __future__ import print_function, absolute_import

from collections import defaultdict
from sage.all import *
from sage.misc.cachefunc import cached_method

from torus_triangulation.remove_one import remove_one_iter
from torus_triangulation.translate_point import TranslatePoint
from torus_triangulation.exceptions import VertexOrderError




class CellCollection(object):
    
    def __init__(self, dimension, simplices, boundary):
        """
        Auxiliary class to build delta complex from simplicial cells
        """
        self.dimension = dimension
        self.simplices = simplices
        self.boundary = boundary
        
    def delta_complex(self):
        delta = dict()
        facet_enumeration = dict()
        for i in range(self.dimension + 2):
            simplices = sorted(self.simplices[i])
            delta[i - 1] = tuple(
                tuple(facet_enumeration[facet] for facet in self.boundary[simplex])
                for simplex in simplices
            )
            facet_enumeration = dict([simplex, pos] for pos, simplex in enumerate(simplices))
        return DeltaComplex(delta)

    def test_boundaries(self):
        for i in range(self.dimension + 2):
            simplices = self.simplices[i]
            for simplex in simplices:
                print(simplex, self.boundary[simplex])
    
    def test_d_d(self, simplex):
        for j in range(len(simplex)):
            for i in range(j):
                didjminus1 = builder.boundary[builder.boundary[s][i]][j-1]
                djdi = builder.boundary[builder.boundary[s][j]][i]
                print(i, j, didjminus1, djdi)

    def test(self):
        for dim, simplices in self.simplices.items():
            for simplex in simplices:
                simplex_bdry = self.boundary[simplex]
                for j in range(dim):
                    for i in range(j):
                        didjminus1 = self.boundary[simplex_bdry[i]][j-1]
                        djdi = self.boundary[simplex_bdry[j]][i]
                        assert djdi == didjminus1, \
                            'Invalid boundary: i={0}, j={1}, simplex={2}'.format(i,j,simplex)
                
    def quotient(self, relation):
        # The quotient "simplices" are the equivalence classes
        quotient_simplices = defaultdict(set)
        for dim in range(self.dimension + 2):
            quotient_simplices[dim] = set(relation(simplex) for simplex in self.simplices[dim])
        # The quotient "boundary" are the induced boundary maps
        quotient_boundary = dict()
        for key, value in self.boundary.items():
            quotient_boundary[relation(key)] = [relation(facet) for facet in value]
        return CellCollection(self.dimension, quotient_simplices, quotient_boundary)

    def homology(self, **kwds):
        return self.delta_complex().homology(**kwds)

    def plot(self):
        return sum(Polyhedron(simplex).plot() for simplex in self.simplices[self.dimension + 1])
    
    

        
        
class SimplexRelation(object):

    def __init__(self, builder):
        self.builder = builder
        self.class_map = dict()

    def _unsafe_identify(self, *simplices):
        """
        Identify simplices ignoring boundary data
        
        Use :meth:`indentify` instead.
        """
        union = set()
        for simplex in simplices:
            union.update(self(simplex))
        union = frozenset(union)
        for simplex in simplices:
            self.class_map[simplex] = union
        
    def identify(self, *simplices):
        self._unsafe_identify(*simplices)
        boundaries = [self.builder.boundary[simplex] for simplex in simplices]
        for boundaries_i in zip(*boundaries):
            self.identify(*boundaries_i)

    def __call__(self, simplex):
        """
        Return the equivalence class of the simplex
        """
        try:
            return self.class_map[simplex]
        except KeyError:
            return frozenset([simplex])

    def test(self):
        """
        Check consistency

        This method tests that the induced vertex order given by the replacement commutes with
        boundary map
        """
        for src, cls in self.class_map.items():
            for dst in cls:
                vertex_map = dict([src_vertex, dst_vertex] for src_vertex, dst_vertex in zip(src, dst))
                src_boundary = self.builder.boundary[src]
                dst_boundary = self.builder.boundary[dst]
                for src_facet, dst_facet in zip(src_boundary, dst_boundary):
                    for src_vertex, dst_vertex in zip(src_facet, dst_facet):
                        assert vertex_map[src_vertex] == dst_vertex, \
                            'Vertex order mismatch in identification: {0} does not map to {1}'.format(
                                src_vertex, dst_vertex)
        
        

        

        
