from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
...
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)

