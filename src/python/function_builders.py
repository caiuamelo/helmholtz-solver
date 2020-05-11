class RegionFunctionBuilder():
    class _RectangleRegion():
        def __init__(self, start_point, end_point, value):
            self.start_point = start_point
            self.end_point = end_point
            self.value = value
        
        def is_point_inside(self, point):
            xini, yini = self.start_point
            xend, yend = self.end_point
            x, y = point
            return (xini < x < xend and yini < y < yend)

    def __init__(self, overall_value):
        self._overall_value = overall_value
        self._region_list = []

    def set_rectangle_interval_value(self, start_point, end_point, value):
        self._region_list.append(
            RegionFunctionBuilder._RectangleRegion(
                start_point,
                end_point,
                value
            )
        )
        return self

    def build(self):
        def f(x, y):
            for region in self._region_list:
                if region.is_point_inside((x, y)):
                    return region.value

            return self._overall_value
        return f
