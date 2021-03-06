"""

    QuadTree data structure

"""

def normalize_rect(rect):
    x1, y1, x2, y2 = rect
    if x1 > x2:
        x1, x2 = x2, x1
    if y1 > y2:
        y1, y2 = y2, y1
    return (x1, y1, x2, y2)


class QuadNode:
    def __init__(self, item, rect):
        self.item = item
        self.rect = rect


class QuadTree:
    MAX = 10
    MAX_DEPTH = 20

    def __init__(self, x, y, size, depth = 0):
        self.nodes = []
        self.children = []
        self.center = [x, y]
        self.size = size
        self.depth = depth

    def insert(self, item, rect):
        rect = normalize_rect(rect)

        if len(self.children) == 0:
            node = QuadNode(item, rect)
            self.nodes.append(node)

            if len(self.nodes) > self.MAX and self.depth < self.MAX_DEPTH:
                self.split()
                return node
        else:
            return self.insert_into_children(item, rect)

    def insert_into_children(self, item, rect):

        # if rect spans center then insert here
        if ((rect[0] <= self.center[0] and rect[2] > self.center[0]) and
            (rect[1] <= self.center[1] and rect[3] > self.center[1])):
            node = QuadNode(item, rect)
            self.nodes.append(node)
            return node
        else:

            # try to insert into children
            if rect[0] <= self.center[0]:
                if rect[1] <= self.center[1]:
                    return self.children[0].insert(item, rect)
                if rect[3] > self.center[1]:
                    return self.children[1].insert(item, rect)
            if rect[2] > self.center[0]:
                if rect[1] <= self.center[1]:
                    return self.children[2].insert(item, rect)
                if rect[3] > self.center[1]:
                    return self.children[3].insert(item, rect)


    def split(self):
        self.children = [QuadTree(self.center[0] - self.size/2,
                                  self.center[1] - self.size/2,
                                  self.size/2, self.depth + 1),
                         QuadTree(self.center[0] - self.size/2,
                                  self.center[1] + self.size/2,
                                  self.size/2, self.depth + 1),
                         QuadTree(self.center[0] + self.size/2,
                                  self.center[1] - self.size/2,
                                  self.size/2, self.depth + 1),
                         QuadTree(self.center[0] + self.size/2,
                                  self.center[1] + self.size/2,
                                  self.size/2, self.depth + 1)]

        nodes = self.nodes
        self.nodes = []
        for node in nodes:
            self.insert_into_children(node.item, node.rect)


    def query(self, rect, results=None):

        if results is None:
            rect = normalize_rect(rect)
            results = set()


        # search children
        if len(self.children) > 0:
            if rect[0] <= self.center[0]:
                if rect[1] <= self.center[1]:
                    self.children[0].query(rect, results)
                if rect[3] > self.center[1]:
                    self.children[1].query(rect, results)
            if rect[2] > self.center[0]:
                if rect[1] <= self.center[1]:
                    self.children[2].query(rect, results)
                if rect[3] > self.center[1]:
                    self.children[3].query(rect, results)

        # search node at this level
        for node in self.nodes:
            if (node.rect[2] > rect[0] and node.rect[0] <= rect[2] and
                node.rect[3] > rect[1] and node.rect[1] <= rect[3]):
                results.add(node.item)

        return results


    def get_size(self):
        size = 0
        for child in self.children:
            size += child.get_size()
        size += len(self.nodes)
        return size


