import sys
from PyQt5 import QtWidgets, QtGui, QtCore
 
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
 
import random
 
class Window(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
 
        # a figure instance to plot on
        self.figure = Figure()
 
        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)
 
        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)
 
        # Just some button connected to `plot` method
        self.button = QtWidgets.QPushButton('Plot')
        self.button.clicked.connect(self.plot)
 
        # set the layout
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(self.button)
        self.setLayout(layout)
 
    def plot(self):
        ''' plot some random stuff '''
        # random data
        data = [random.random() for i in range(10)]
 
        # create an axis
        ax = self.figure.add_subplot(111)
 
        # discards the old graph
        ax.clear()
 
        # plot data
        ax.plot(data, '-')

        # plot decoration
        ax.set_title()

        self.figure.tight_layout()
        # refresh canvas
        self.canvas.draw()
 
class Principal(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
 
        self.resize(1560, 960)
 
        self.window = Window(self)
        self.window.setGeometry(30, 30, 1500, 900)
 
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
 
    main = Principal()
    main.show()
 
    sys.exit(app.exec_())