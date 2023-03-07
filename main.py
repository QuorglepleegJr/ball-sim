# Kivy Imports

from kivy.app import App
from kivy.uix.widget import Widget
from kivy.vector import Vector
from kivy.properties import NumericProperty, ObjectProperty, \
    ReferenceListProperty, ListProperty, BooleanProperty
from kivy.clock import Clock
from kivy.lang import Builder

# Other Imports

from sys import exit
from math import sin, cos, sqrt

# Widgets

class SimulationBall(Widget):

    '''
    Modelled as a smooth circle, this class represents an object
    to be simulated under gravity and physics.
    '''

    # Class Constants

    GRAVITY = 0 #-500 - Debug removal for the time being
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

        '''
        Acts as a factory function of sorts, calling the relevant
        collision handle method for all relevant objects.
        '''

        if isinstance(other, SimulationBall):

            self.ball_collision_check(other, delta)

        elif isinstance(other, SimulationBlock):

            self.block_collision_check(other, delta)
    
    def ball_collision_check(self, ball, delta):

        '''
        Handles the collision between two balls.
        '''

        # Manual method because self.collide_widget(ball) appeared to not work
        # And it makes the actual collisions easier

        p1 = Vector(self.pos)
        p2 = Vector(ball.pos)
        v1 = Vector(self.vel) * delta
        v2 = Vector(ball.vel) * delta

        if (v2-v1).length() != 0:

            t = None

            v = v2-v1
            p = p2-p1
            r = self.radius + ball.radius

            a = v.dot(v)
            b = 2*(p.dot(v))
            c = p.dot(p)-r*r

            disc = b**2 - 4*a*c

            print(a,b,c,disc)

            if disc >= 0:

                temp_sol = (b+sqrt(disc))/(2*a)
                solution_1 = -temp_sol # Tends to lower floatin point error
                solution_2 = c/solution_1 

                print(a,b,c,disc,solution_1,solution_2)

                if solution_1 >= 0 and solution_1 < 1:

                    t = solution_1
                
                elif solution_2 >= 0 and solution_2 < 1:

                    t = solution_2

            if t is not None:

                # Finding the vectors of the velocities component to positions
                # At the instant of bouncing to ensure the right bit is scaled

                v = (p2 + v2 * t) - (p1 + v1 * t)

                # Parallel components scaled by delta

                v1_component_parallel = v * (v1.dot(v)/v.dot(v))
                v2_component_parallel = v * (v2.dot(v)/v.dot(v))
                
                # Perpendicular components unscaled

                self_vel_parallel = Vector(v * \
                        (Vector(self.vel).dot(v)/v.dot(v)))
                ball_vel_parallel = Vector(v * \
                        (Vector(ball.vel).dot(v)/v.dot(v)))

                self_vel_perp = Vector(self.vel) - self_vel_parallel
                ball_vel_perp = Vector(ball.vel) - ball_vel_parallel

                print(v1, v2, p1, p2, v, v1_component_parallel, self_vel_perp, v2_component_parallel, ball_vel_perp, t) # Debug

                # Move to position to touch

                self.pos[0] += v1_component_parallel[0] * t
                self.pos[1] += v1_component_parallel[1] * t
                ball.pos[0] += v2_component_parallel[0] * t
                ball.pos[1] += v2_component_parallel[1] * t

                # Remove parallel components

                self.vel = self_vel_perp
                ball.vel = ball_vel_perp

    def block_collision_check(self, block, delta):

        '''
        Handles collisions between balls and blocks, UNFINISHED.
        '''

        # Manual method because self.collide_widget(ball) appeared to not work

        centers_vec = Vector(self.pos) + Vector(self.vel) *\
              delta - Vector(block.pos)

        # New Idea for method - model ball as a point on graph origin block.pos
        # with axes along rotation of block, then compute y and x in model
        # and check if crossing boundaries

        if centers_vec.length() < self.radius + \
            block.calculate_radius(centers_vec.angle((1,0))):

            print(f"{self} and {block} have collided!")

    def update_before_collision(self, delta):

        '''
        Performs any updates required before collision detection applied,
        notably (and as of now only) accelerations.
        '''

        self.vel[1] += SimulationBall.GRAVITY * delta

    def update_after_collision(self, delta):

        self.pos[0] += self.vel[0] * delta
        self.pos[1] += self.vel[1] * delta

class SimulationBlock(Widget):

    '''
    Modelled as a smooth rectangle, this class represents an
    immovable (and non-momentum-conserving) surface.
    '''

    # Properties

    x_pos = NumericProperty(0)
    y_pos = NumericProperty(0)
    pos = ReferenceListProperty(x_pos, y_pos)
    x_size = NumericProperty(150)
    y_size = NumericProperty(50)
    size = ReferenceListProperty(x_size, y_size)
    theta = NumericProperty(0) # Angle from positive x-axis, between 0 and 2 Pi

    # Methods

    def calculate_radius(self, phi) -> float: # Phi used to avoid confusion with theta

        '''
        Calculates the effective radius at a given angle from the positive
        x-axis, phi.
        '''

        phi -= self.theta # Account for rotated blocks

        # Mathematical formula for "radius" of 1x1 square
        # Simply written as max(sec(theta), csc(theta))
        # 0.5 works in terms of radius rather than diameter

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

    '''
    The core class of the simulation, this acts as the window
    and parent of all the other objects.
    '''

    # Properties

    balls = ListProperty()
    blocks = ListProperty()

    # Methods

    def initialise(self, balls=None, blocks=None):

        '''
        Add lists of balls and blocks as provided, in format 
        (object, pos_x, pos_y, vel_x, vel_y).
        '''

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

        '''
        Runs through every object, updating it as required.
        Also triggers objects to check collisions.
        Ran every frame.
        '''

        # Three stage update - handle gravity, collisions, then finally move

        # Stage 1

        for ball in self.balls:

            ball.update_before_collision(delta)

        # Stage 2
        
        ball_pairs = [{a, b} for a in self.balls for b in self.balls]

        while len(ball_pairs) > 0:

            balls = ball_pairs.pop()

            if balls not in ball_pairs:
            
                balls = list(balls)
                if len(balls) == 2:
                    balls[0].collision_check(balls[1], delta)

        # Stage 3
        
        for ball in self.balls:

            ball.update_after_collision(delta)

    def on_touch_down(self, touch):

        self.update(1/60)

# Applications

class SimulationApp(App):

    '''
    Functions as the wrapper for the SimulationManager.
    '''

    def build(self):

        Builder.load_file("simulation.kv")

        self.window = SimulationManager()

        if hasattr(self, "balls") and hasattr(self, "blocks"):

            init_balls = [(SimulationBall(), *b) for b in self.balls]
            init_blocks = [(SimulationBlock(), *b) for b in self.blocks]

            self.window.initialise(balls=init_balls, blocks=init_blocks)

            delattr(self, "balls")
            delattr(self, "blocks")
        
        if not hasattr(self, "frame_advance") or not self.frame_advance:

            Clock.schedule_interval(self.window.update, 1/60)

        return self.window

    def initialise(self, balls=None, blocks=None, frame_advance=False):

        '''
        Sets the initial state of the simulation.
        '''
        
        self.balls = balls
        self.blocks = blocks

        self.frame_advance = frame_advance

# Main

if __name__ == "__main__":

    sim = SimulationApp()

    sim.initialise(
        balls=((200, 500, 0, -100),
               (500, 500, 0, -100)),
        blocks=((200, 200),
                (500, 200)),
        frame_advance=False)
    
    sim.run()
    