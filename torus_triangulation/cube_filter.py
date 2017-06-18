

class CubeFilter(object):

    def __init__(self, *sets):
        self.sets = tuple(map(tuple, sets))

    def contains(self, vertex):
        for i, values in enumerate(self.sets):
            if vertex[i] not in values:
                return False
        return True
        
    def __call__(self, vertex_list):
        return [v for v in vertex_list if self.contains(v)]

