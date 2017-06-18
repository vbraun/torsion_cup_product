from __future__ import print_function, absolute_import
from sage.all import (
    ZZ,
    Permutations,
    cartesian_product_iterator,
    cached_method,
    Polyhedron
)



class TriangulationBase(object):

    def plot(self):
        """
        EXAMPLES::

            sage: UnitCubeTriangulation(3).plot()
        """
        return sum(Polyhedron(simplex).plot() for simplex in self.simplex_iter())
    

    
class UnitCubeTriangulation(TriangulationBase):

    def __init__(self, dimension):
        """
        Reference triangulation of the `[0,1]^n` hypecube

        There is probably not too much sense in the choice of the reference triangulation; We just
        have to pick one. The merit of this one is that its is rather simple to implement.

        INPUT:

        - ``dimension`` -- integer. The dimension

        EXAMPLES::

            sage: square = UnitCubeTriangulation(2)
            sage: list(square.simplex_iter())
            [((0, 0), (1, 0), (1, 1)), 
             ((0, 0), (0, 1), (1, 1))]
        """
        self.dimension = dimension
        self.permutations = Permutations(self.dimension)

    @cached_method
    def points(self):
        """
        The vertices in a fixed order (lexicographic sorted)

        EXAMPLES::

            sage: square = UnitCubeTriangulation(2)
            sage: square.points()
            ((0, 0), (0, 1), (1, 0), (1, 1))
        """
        return tuple(sorted(
            tuple(p) for p in cartesian_product_iterator([[0,1]] * self.dimension)
        ))
        
    def simplex_iter(self):
        for p in self.permutations:
            vertices = []
            for i in range(self.dimension + 1):
                v = [1] * self.dimension
                for j in p[i:]:
                    v[j-1] = 0
                vertices.append(tuple(v))
            yield(tuple(vertices))

    __iter__ = simplex_iter

    def len(self):
        return self.permutations.cardinality()



    

class CubeTriangulation(TriangulationBase):

    def __init__(self, *ordered_vertices):
        """
        Triangulation of cube with prescribed vertex order

        That is, the vertices of each simplex will be in the given order.

        EXAMPLES::

            sage: square = CubeTriangulation((1, 2), (4, 3), (1, 3), (4, 2))
            sage: square.dimension
            2
            sage: square.vertices
            ((1, 2), (4, 3), (1, 3), (4, 2))
            sage: square.extent_min
            (1, 2)
            sage: square.extent_max
            (4, 3)
            sage: square.unit_vertex_map
            {(0, 0): (1, 2), (0, 1): (1, 3), (1, 0): (4, 2), (1, 1): (4, 3)}
        """
        self.vertices = tuple(tuple(v) for v in ordered_vertices)
        self.sort_key = dict(zip(self.vertices, range(len(self.vertices))))
        self.dimension = len(self.vertices[0])
        self.extent_min = tuple(map(min, zip(*self.vertices)))
        self.extent_max = tuple(map(max, zip(*self.vertices)))
        self.unit = UnitCubeTriangulation(self.dimension)
        self.unit_vertex_map = dict(zip(self.unit.points(), self.points()))

    @cached_method
    def points(self):
        """
        The vertices in a fixed order (lexicographic sorted)

        EXAMPLES::

            sage: square = CubeTriangulation((1,2), (3,2), (1,4), (3,4))
            sage: square.points()
            ((1, 2), (1, 4), (3, 2), (3, 4))
        """
        return tuple(sorted(
            tuple(p) for p in cartesian_product_iterator(zip(self.extent_min, self.extent_max))
        ))
        
    def simplex_iter(self):
        """
        Iterate over the simplices of a cube triangulation

        EXAMPLES::
    
            sage: square = CubeTriangulation((1,2), (4,3), (1,3), (4,2))
            sage: simplex_iter = iter(square)
            sage: next(simplex_iter)
            ((1, 2), (4, 3), (4, 2))
            sage: next(simplex_iter)
            ((1, 2), (4, 3), (1, 3))
            sage: next(simplex_iter)
            Traceback (most recent call last):
            ...
            StopIteration:
        """
        for unit_simplex in self.unit.simplex_iter():
            simplex = [self.unit_vertex_map[v] for v in unit_simplex]
            simplex.sort(key=lambda p: self.sort_key[p])
            yield tuple(simplex)
        
    __iter__ = simplex_iter

    def len(self):
        return len(self.unit)

    def plot(self):
        """
        EXAMPLES::

            sage: UnitCubeTriangulation(3).plot()
        """
        return sum(Polyhedron(simplex).plot() for simplex in self.simplex_iter())
    


    



            
