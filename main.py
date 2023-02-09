# Kivy Imports

from kivy.app import App
from kivy.uix.widget import Widget
from kivy.vector import Vector
from kivy.properties import NumericProperty, \
    ReferenceListProperty, ObjectProperty
from kivy.clock import Clock

# Other Imports

from sys import exit

# Widgets

class SimulationBall(Widget):

    # Class Constants

    GRAVITY = -10

    # Properties

    x_pos = NumericProperty(0)
    y_pos = NumericProperty(0)
    pos = ReferenceListProperty(x_pos, y_pos)
    x_vel = NumericProperty(0)
    y_vel = NumericProperty(0)
    vel = ReferenceListProperty(x_vel, y_vel)
    mass = NumericProperty(0)
    radius = NumericProperty(25)

    # Methods

    def update(self, delta):

        self.vel[1] += SimulationBall.GRAVITY * delta

        self.pos[0] += self.vel[0]
        self.pos[1] += self.vel[1]

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
    
    # DEBUG PROPERTIES

    debug_ball = ObjectProperty(None)
    debug_block = ObjectProperty(None)

    # Methods

    def update(self, delta):
        
        self.debug_ball.update(delta)

# Applications

class SimulationApp(App):

    def build(self):

        window = SimulationWindow()

        Clock.schedule_interval(window.update, 1/60)

        return window

# Main

if __name__ == "__main__":
    exit(SimulationApp().run())