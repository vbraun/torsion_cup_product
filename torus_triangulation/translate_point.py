


class TranslatePoint(object):
    
    def __init__(self, coordinate_index, amount):
        """
        Translate point (tuple of numbers) in a coordinate direction
        """
        self._i = coordinate_index
        self._d = amount
        
    def __call__(self, point):
        result = list(point)
        result[self._i] += self._d
        return tuple(result)

    def __nonzero__(self):
        return self._d != 0

    
