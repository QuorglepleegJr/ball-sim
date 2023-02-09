# Kivy Imports

from kivy.app import App
from kivy.uix.widget import Widget
from kivy.vector import Vector
from kivy.properties import NumericProperty, \
    ReferenceListProperty, ObjectProperty, ListProperty
from kivy.clock import Clock

# Other Imports

from sys import exit
from math import sin, cos

# Widgets

class SimulationBall(Widget):

    # Class Constants

    GRAVITY = -500

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

    def collision_check(self, other):

        if isinstance(other, SimulationBall):

            self.ball_collision_check(other)

        elif isinstance(other, SimulationBlock):

            self.block_collision_check(other)
    
    def ball_collision_check(self, ball):

        # Manual method because self.collide_widget(ball) appeared to not work
        if (Vector(self.pos) - Vector(ball.pos)).length() <  \
            self.radius + ball.radius:

            print(f"{self} and {ball} have collided!")

    def block_collision_check(self, block):

        # Manual method because self.collide_widget(ball) appeared to not work

        centers_vec = Vector(self.pos) - Vector(block.pos)

        if centers_vec.length() < self.radius + \
           block.calculate_radius(centers_vec.angle((1,0))):

            print(f"{self} and {block} have collided!")

    def update_before_collision(self, delta):

        self.vel[1] += SimulationBall.GRAVITY * delta

    def update_after_collision(self, delta):

        self.pos[0] += self.vel[0] * delta
        self.pos[1] += self.vel[1] * delta

class SimulationBlock(Widget):

    # Properties

    x_pos = NumericProperty(0)
    y_pos = NumericProperty(0)
    pos = ReferenceListProperty(x_pos, y_pos)
    x_size = NumericProperty(150)
    y_size = NumericProperty(50)
    size = ReferenceListProperty(x_size, y_size)
    theta = NumericProperty(0)

    # Methods

    def calculate_radius(self, phi):

        phi -= self.theta # Account for rotated blocks

        # Mathematical formula for "radius" of 1x1 square

        if sin(phi) == 0:

            r = 0.5/cos(phi)

        elif cos(phi) == 0:

            r = 0.5/sin(phi)
        
        else:

            r = 0.5 * min(abs(1/sin(phi)), abs(1/cos(phi)))

        # Find vector equivalent to required radius, origin center of block

        radial_vec = Vector(r, 0).rotate(phi)

        radial_vec.x *= self.x_size
        radial_vec.y *= self.y_size

        return radial_vec.length()

class SimulationManager(Widget):

    # Properties

    balls = ListProperty()
    blocks = ListProperty()

    # Methods

    def initialise(self, balls=None, blocks=None):

        if balls is not None:

            for ball in balls:

                ball[0].pos = ball[1:]
                self.add_widget(ball[0])
                self.balls.append(ball[0])
        
        if blocks is not None:

            for block in blocks:

                block[0].pos = block[1:]
                self.add_widget(block[0])
                self.blocks.append(block[0])

    def update(self, delta):

        # Three stage update - handle gravity, collisions, then finally move

        for ball in self.balls:

            ball.update_before_collision(delta)

        for ball in self.balls:
            
            for obj in self.balls + self.blocks:

                # Working around ListProperty's shallow copies
                if id(obj) != id(ball):

                    ball.collision_check(obj)
        
        for ball in self.balls:

            ball.update_after_collision(delta)

        # DEBUG LINE - TESTING BALL COLLISION

        #self.balls[1].vel[1] -= SimulationBall.GRAVITY * delta

# Applications

class SimulationApp(App):

    def build(self):

        window = SimulationManager()

        window.initialise(balls=((SimulationBall(), 500, 500),  \
            (SimulationBall(), 700, 200)),  \
            blocks=((SimulationBlock(), 500, 200),  \
            (SimulationBlock(), 300, 200)))

        Clock.schedule_interval(window.update, 1/60)

        return window

# Main

if __name__ == "__main__":
    exit(SimulationApp().run())