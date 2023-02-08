# Kivy Imports

from kivy.app import App
from kivy.uix.widget import Widget
from kivy.vector import Vector
from kivy.properties import NumericProperty, \
    ReferenceListProperty, ObjectProperty

# Other Imports

from sys import exit

# Widgets

class SimulationBall(Widget):

    # Properties

    x_pos = NumericProperty(0)
    y_pos = NumericProperty(0)
    pos = ReferenceListProperty(x_pos, y_pos)
    x_vel = NumericProperty(0)
    y_vel = NumericProperty(0)
    vel = ReferenceListProperty(x_vel, y_vel)
    mass = NumericProperty(0)
    radius = NumericProperty(25)

class SimulationBlock(Widget):

    # Properties

    x_pos = NumericProperty(0)
    y_pos = NumericProperty(0)
    pos = ReferenceListProperty(x_pos, y_pos)
    x_size = NumericProperty(150)
    y_size = NumericProperty(50)
    size = ReferenceListProperty(x_size, y_size)
    theta = NumericProperty(0)

class SimulationWindow(Widget):

    debug_ball = ObjectProperty(None)
    debug_block = ObjectProperty(None)

# Applications

class SimulationApp(App):

    def build(self):

        window = SimulationWindow()

        return window

# Main

if __name__ == "__main__":
    exit(SimulationApp().run())