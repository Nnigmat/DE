import sys
from math import exp, isnan, isinf, inf
from numpy import linspace, concatenate, subtract
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout, QLabel, QLineEdit, QComboBox

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt


class NumericalProcedures():

    def __init__(self, x0, X, n, y0):
        self.x0 = x0 # Initial Values
        self.y0 = y0

        self.X = X # X end value
        self.n = n # Number of steps
        self.constant = 1/(self.y0 + exp(self.x0)) + self.x0 # Constant of exact solution
        self.h = (self.X - self.x0) / self.n # Height of one step
        self.set_x_axis(self.n)

    def set_x_axis(self, n):
        self.h = (self.X - self.x0) / n
        steps = int(abs(self.x0 - self.constant) / abs(self.x0 - self.X) * n)
        if steps != 0:
            self.x = concatenate((linspace(self.x0, self.constant - 0.1, steps), linspace(self.constant + 0.1, self.X, n - steps)))
        else:
            self.x = linspace(self.x0, self.X, n)

        # X axis with points
        for i in range(1, len(self.x[1:])):
            if self.x[i-1] < self.constant and self.x[i] > self.constant:
                self.breakpoint = i
                break

    def __repr__(self):
        return 'x0: {}, X: {}, n: {}, y0: {}, h: {}'.format(self.x0, self.X, self.n, self.y0, self.h)

    def formula(self, x, y):
        return (1-2*y)*exp(x) + y * y + exp(2*x)

    def eulers_method(self):
        prev_y = self.y0
        data = [prev_y]
        for i in range(1, len(self.x)):
            if i == self.breakpoint:
                prev_y = 1 / (self.constant - self.x[i]) + exp(self.x[i])  # We have jumped over breakpoint and now we use new y0 which we get from exact solution
            else:
                prev_y = prev_y + self.h * self.formula(self.x[i-1], prev_y)
            data.append(prev_y)
        return self.x, data

    def improved_method(self):
        prev_y = self.y0
        data = [self.y0]
        for i in range(1, len(self.x)):
            if i == self.breakpoint:
                prev_y = 1 / (self.constant - self.x[i]) + exp(self.x[i])  # We have jumped over breakpoint and now we use new y0 which we get from exact solution
            else:
                temp_x = self.x[i-1] + self.h / 2
                temp_y = prev_y + self.h / 2 * prev_y
                delta_y = self.h * self.formula(temp_x, temp_y)
                prev_y = prev_y + delta_y
            data.append(prev_y)
        return self.x, data

    def runge_kutta_method(self):
        prev_y = self.y0
        data = [self.y0]
        for i in range(1, len(self.x)):
            if i == self.breakpoint:
                prev_y = 1 / (self.constant - self.x[i]) + exp(self.x[i])  # We have jumped over breakpoint and now we use new y0 which we get from exact solution
            else:
                prev_x = self.x[i-1]
                k1 = self.formula(prev_x, prev_y)
                k2 = self.formula(prev_x + self.h / 2, prev_y + self.h * k1 / 2)
                k3 = self.formula(prev_x + self.h / 2, prev_y + self.h * k2 / 2)
                k4 = self.formula(prev_x + self.h, prev_y + self.h * k3)
                prev_y = prev_y + self.h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            data.append(prev_y)
        return self.x, data

    def exact_solution(self):
        exact = []
        for el in self.x:
            exact.append(1 / (self.constant - el) + exp(el))
        return exact

    def local_error(self, data):
        return list(range(len(data))), [abs(data[i] - self.exact_solution()[i]) for i in range(len(data))]

    def global_error(self, method='Euler'):
        if method == 'Euler':
            function = self.eulers_method
        elif method == 'Improved':
            function = self.improved_method
        else:
            function = self.runge_kutta_method

        error = []
        for i in range(30, 520, 10):
            self.set_x_axis(i)
            _, data = function()
            _, local_error = self.local_error(data)
            error.append(max(local_error))

        return list(range(30, 520, 10)), error


class Window(QDialog):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)

        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Just some button connected to `plot` method
        self.button = QPushButton('Plot')
        self.button.clicked.connect(self.plot)

        # set the layout
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Formula: (1-2y)e^x + y^2 + e^(2x)"))

        self.X0 = QLineEdit()
        self.Y0 = QLineEdit()
        self.N = QLineEdit()
        self.X0.setText('-5')
        self.Y0.setText('2')
        self.N.setText('100')

        layout.addWidget(QLabel("X0:"))
        layout.addWidget(self.X0)
        layout.addWidget(QLabel("Y0:"))
        layout.addWidget(self.Y0)
        layout.addWidget(QLabel("Step:"))
        layout.addWidget(self.N)

        self.comboBox = QComboBox(self)
        self.comboBox.addItems(["Euler's method", "Improved Euler's method", "Runge-Kutta method"])

        self.comboBox.setMinimumWidth(150)

        self.comboBox.move(500, 30)

        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(self.button)
        self.setLayout(layout)

    def plot(self):
        X = 0
        error = False
        try:
            x0 = int(self.X0.text())
        except:
            error = True
            self.X0.setText("Invalid X0 value (used default: -5)")
            x0 = -5

        try:
            y0 = int(self.Y0.text())
        except:
            error = True
            self.X0.setText("Invalid Y0 value (used default: 2)")
            y0 = 2

        try:
            n = int(self.N.text())
        except:
            error = True
            self.X0.setText("Invalid number of steps (used default: 100)")
            n = 100

        if x0 >= X:
            error = True
            self.X0.setText("X0 must be less than 5 (used default: -5)")
            x0 = -5

        if error:
            return

        procedures = NumericalProcedures(x0, X, n, y0)
        method = ''

        if self.comboBox.currentText() == "Euler's method":
            x, data = procedures.eulers_method()
            method = 'Euler'
        elif self.comboBox.currentText() == "Improved Euler's method":
            x, data = procedures.improved_method()
            method = 'Improved'
        elif self.comboBox.currentText() == "Runge-Kutta method":
            x, data = procedures.runge_kutta_method()

        error_x, error = procedures.local_error(data=data)

        # clear from previous values
        self.figure.clear()

        # method plot
        a1 = self.figure.add_subplot(311)
        data, = a1.plot(x, data, label=self.comboBox.currentText(), linewidth=2.5)
        a1.set_xlabel('x', fontsize=15)
        a1.set_ylabel('y', fontsize=15)
        plt.legend(handles=[data])

        # local error plot
        a2 = self.figure.add_subplot(312)
        total, = a2.plot(error_x, error, label='Local error', color='r', linewidth=2.5)
        a2.set_xlabel('step', fontsize=15)
        a2.set_ylabel('local error', fontsize=15)
        a2.text(n / 3.5 , max(error) / 3, 'Average approximation error: {}'.format(round(sum(error) / len(error), 3)), style='italic', bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
        plt.legend(handles=[total])

        # global
        a3 = self.figure.add_subplot(313)
        glob_x, global_error = procedures.global_error(method=method)
        global_err, = a3.plot(glob_x, global_error, label='Global error', color='y', linewidth=2.5)
        a3.set_xlabel('step', fontsize=15)
        a3.set_ylabel('global error', fontsize=15)
        plt.legend(handles=[global_err])

        self.canvas.draw()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main = Window()
    main.show()
    sys.exit(app.exec_())
