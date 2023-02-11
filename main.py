# Kivy Imports

from kivy.app import App
from kivy.uix.widget import Widget
from kivy.vector import Vector
from kivy.properties import NumericProperty, \
    ReferenceListProperty, ListProperty
from kivy.clock import Clock

# Other Imports

from sys import exit
from math import sin, cos

# Widgets

class SimulationBall(Widget):

    # Class Constants

    GRAVITY = -500
    COLLISION_PRECISION = 5

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

    def collision_check(self, other, delta):

        if isinstance(other, SimulationBall):

            self.ball_collision_check(other, delta)

        elif isinstance(other, SimulationBlock):

            self.block_collision_check(other, delta)
    
    def ball_collision_check(self, ball, delta):

        # Manual method because self.collide_widget(ball) appeared to not work
        # And it makes the actual collisions easier

        p1 = Vector(self.pos)
        p2 = Vector(ball.pos)
        v1 = Vector(self.vel) * delta
        v2 = Vector(ball.vel) * delta

        # Calculate scaler along velocities until balls centers overlap

        print((p1-p2).length()/(v2-v1).length())
        
        t = (p1-p2).length()/(v2-v1).length()

        if t < 0 or t > 1:
            return

        ''' Finding the t value that makes the balls just touch:
         This has an algebraical approach leading to a quadratic
         in t with horrible coefficients. This would be awkward 
         and prone to floating point error. Thus, the quadratic
         solution required (range (0,t)) is approximated using
         the Newton-Raphson method, starting with t as x0'''

        p_vec = p1-p2
        v_vec = v1-v2

        a = p_vec.x
        b = v_vec.x
        c = p_vec.y
        d = v_vec.y

        a2 = a**2
        b2 = b**2
        c2 = c**2
        d2 = d**2

        ab = a*b
        cd = c*d

        r2 = (self.radius + ball.radius)**2

        xn = t

        for iteration in range(SimulationBall.COLLISION_PRECISION):

            fxn = (b2+d2)*xn**2 + 2*(ab+cd)*xn + r2-a2-c2
            fdashxn = 2*(b2+d2)*xn + 2*(ab+cd)

            print(p1, p2, p_vec, v1, v2, v_vec, a2, b2, c2, d2, ab, cd, r2, xn, fxn, fdashxn, 2*(b2+d2)*xn, 2*(ab+cd))

            xn -= fxn/fdashxn

        # xn should now be an approximation to the required value

        self.vel.x *= xn
        self.vel.y *= xn
        ball.vel.x *= xn
        ball.vel.y *= xn

        # SOMEHOW fdashxn IS TURNING OUT ZERO, FIX


    def block_collision_check(self, block, delta):

        # Manual method because self.collide_widget(ball) appeared to not work

        centers_vec = Vector(self.pos) + Vector(self.vel) * delta - Vector(block.pos)

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

                ball[0].pos = ball[1:3]
                ball[0].vel = ball[3:]
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

                    ball.collision_check(obj, delta)
        
        for ball in self.balls:

            ball.update_after_collision(delta)

# Applications

class SimulationApp(App):

    def build(self):

        window = SimulationManager()

        window.initialise(balls=((SimulationBall(), 400, 300, 100, 500),  \
            (SimulationBall(), 700, 300, -100, 500)),  \
            blocks=((SimulationBlock(), 500, 200), \
            (SimulationBlock(), 500, 100)))

        Clock.schedule_interval(window.update, 1/60)

        return window

# Main

if __name__ == "__main__":
    exit(SimulationApp().run())