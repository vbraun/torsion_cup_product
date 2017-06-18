from __future__ import print_function, absolute_import



def remove_one_iter(sequence):
    """
    Iterate over sequence while removing one element

    INPUT:

    - ``sequence`` -- list, tuple, or anything that supports slice operator. The sequence to use.

    EXAMPLES::

        sage: from torus_triangulation.remove_one import remove_one_iter
        sage: roi = remove_one_iter(('a', 'b', 'c'))
        sage: next(roi)
        ('b', 'c')
        sage: next(roi)
        ('a', 'c')
        sage: next(roi)
        ('a', 'b')
        sage: next(roi)
        Traceback (most recent call last):
        ...
        StopIteration
        sage: list(remove_one_iter(range(5)))
        [[1, 2, 3, 4], [0, 2, 3, 4], [0, 1, 3, 4], [0, 1, 2, 4], [0, 1, 2, 3]]
    """
    for i in range(len(sequence)):
        yield sequence[0:i] + sequence[i+1:]

