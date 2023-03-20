# Kivy Imports

from kivy.app import App
from kivy.uix.widget import Widget
from kivy.vector import Vector
from kivy.properties import NumericProperty, ObjectProperty, \
    ReferenceListProperty, ListProperty, BooleanProperty
from kivy.clock import Clock
from kivy.lang import Builder

# Other Imports

from math import sin, cos, sqrt, degrees, floor

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
    adjust_x = NumericProperty(0)
    adjust_y = NumericProperty(0)
    adjustments = ReferenceListProperty(adjust_x, adjust_y)
    x_vel = NumericProperty(0)
    y_vel = NumericProperty(0)
    vel = ReferenceListProperty(x_vel, y_vel)
    mass = NumericProperty(0)
    radius = NumericProperty(25)

    # Methods

    def collision_check(self, other, delta):

        '''
        Acts as a factory function of sorts, calling the relevant
        collision detection method for all relevant objects.
        '''

        if isinstance(other, SimulationBall):

            return self.ball_collision_check(other, delta)

        elif isinstance(other, SimulationBlock):

            return self.block_collision_check(other, delta)
    
    def handle_collision(self, other, t, delta):

        '''
        Acts as a factory function of sorts, calling the relevant
        collision handle method for all relevant objects.
        '''

        if isinstance(other, SimulationBall):

            self.handle_ball_collision(other, t, delta)

        elif isinstance(other, SimulationBlock):

            print("COLLIDING")

            self.handle_block_collision(other, t, delta)
    
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

        t = None

        if (v2-v1).length() != 0:

            v = v2-v1
            p = p2-p1
            r = self.radius + ball.radius

            a = v.dot(v)
            b = 2*(p.dot(v))
            c = p.dot(p) - r * r

            disc = b**2 - 4*a*c

            if disc >= 0:

                temp_sol = (b+sqrt(disc))/(2*a)
                solution_1 = -temp_sol # Tends to lower floating point error
                solution_2 = 0
                if solution_1 != 0:
                    solution_2 = c/solution_1 

                if solution_1 >= 0 and solution_1 < 1:

                    t = solution_1
                
                elif solution_2 >= 0 and solution_2 < 1:

                    t = solution_2

        return t

    def handle_ball_collision(self, ball, t, delta):

        '''
        Handles the collision of balls and blocks
        '''
            
        p1 = Vector(self.pos)
        p2 = Vector(ball.pos)
        v1 = Vector(self.vel) * delta
        v2 = Vector(ball.vel) * delta

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
        Handles detecting collisions between two balls
        '''

        # Manual method because self.collide_widget(ball) appeared to not work

        vel = Vector(self.vel) * delta

        start_pos = Vector(self.pos) - Vector(block.pos)
        end_pos = start_pos + vel

        x_axis_vec = Vector(1, 0).rotate(degrees(block.theta))
        y_axis_vec = x_axis_vec.rotate(90)

        # Find the equivalent positions in the vector space with
        # origin at center of block and axes parallel and perpendicular
        # to the rotation of the block

        aligned_start = Vector(start_pos.dot(x_axis_vec), \
                               start_pos.dot(y_axis_vec))
        aligned_end = Vector(end_pos.dot(x_axis_vec), \
                             end_pos.dot(y_axis_vec))
        aligned_vel = aligned_end-aligned_start

        t1 = 1
        t2 = 1

        # a ^ b < 0 returns whether signs are opposite via bitwise xor
        # Flooring casts to int whilst preserving sign

        if aligned_vel.x != 0 and \
            floor(aligned_start.x) ^ floor(aligned_vel.x) < 0:

            t1 = (abs(aligned_start.x) - (self.radius+block.size[0]/2))\
                / abs(aligned_vel.x)
            
        print(abs(aligned_start.x), (self.radius + block.size[0]/2), abs(aligned_vel.x), self.radius, block.size[0])
        
        if aligned_vel.y != 0 and \
            floor(aligned_start.y) ^ floor(aligned_vel.y) < 0:

            t2 = (abs(aligned_start.y) - (self.radius+block.size[1]/2))\
                / abs(aligned_vel.y)

        #print(t1, t2, x_axis_vec, y_axis_vec, start_pos, end_pos, vel, aligned_start, aligned_end, aligned_vel)

        valid = [t for t in (t1, t2) if t < 1 and t >= 0]

        if len(valid) > 0:

            return min(valid)
        
        return None
    
    def handle_block_collision(self, block, t, delta):

        '''
        Handles collisions between balls and blocks
        '''

        print("Colliding at", t)

        vel = Vector(self.vel) * delta

        start_pos = Vector(self.pos) - Vector(block.pos)
        end_pos = start_pos + vel

        x_axis_vec = Vector(1, 0).rotate(degrees(block.theta))
        y_axis_vec = x_axis_vec.rotate(90)

        aligned_start = Vector(start_pos.dot(x_axis_vec), \
            start_pos.dot(y_axis_vec))
        aligned_end = Vector(end_pos.dot(x_axis_vec), \
            end_pos.dot(y_axis_vec))
        aligned_vel = aligned_end - aligned_start

        required_index = 0 # Stores whether x or y should be modified, x=0, y=1
        
        if aligned_vel.x != 0 and \
            floor(aligned_start.x) ^ floor(aligned_vel.x) < 0:

            if t == (abs(aligned_start.x) - (self.radius+block.size[0]/2))\
                / abs(aligned_vel.x):

                required_index = 1 # Passing via x, so y axis
        
        axis = [x_axis_vec, y_axis_vec][required_index]

        perpendicular = (Vector(self.vel).dot(axis)\
            / axis.dot(axis)) * axis
        unscaled_parallel = Vector(self.vel) - perpendicular
        scaled_parallel = unscaled_parallel * delta

        print(axis, self.vel, scaled_parallel, unscaled_parallel, perpendicular)

        self.adjustments = scaled_parallel * t * 2

        self.vel = perpendicular - unscaled_parallel

        print(t, x_axis_vec, y_axis_vec, start_pos, end_pos)

    def update_before_collision(self, delta):

        '''
        Performs any updates required before collision detection applied,
        notably (and as of now only) accelerations.
        '''

        self.vel[1] += SimulationBall.GRAVITY * delta

    def update_after_collision(self, delta):

        '''
        Performs updates such as applying velocity and position adjustments
        after all the calculations to modify them have been completed.
        '''

        self.pos[0] += self.vel[0] * delta + self.adjustments[0]
        self.pos[1] += self.vel[1] * delta + self.adjustments[1]

        self.adjustments = [0,0]

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

    def collision_check(self, other, delta):

        if isinstance(other, SimulationBall):

            other.collision_check(self, delta)

        raise ValueError("Must collide with ball")
    
    def handle_collision(self, other, t, delta):

        if isinstance(other, SimulationBall):

            other.handle_collision(self, t, delta)

        raise ValueError("Must collide with ball")

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
        (pos_x, pos_y, vel_x, vel_y).
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

        # Four stage update - handle gravity, obtain collisions, 
        # handle collisions in appropriate order, then finally move

        # Stage 1

        for ball in self.balls:

            ball.update_before_collision(delta)

        # Stage 2 & 3 until resolved all collisions

        collisions = self.get_all_collisions(delta)

        while len(collisions) > 0:

            current = collisions.pop()

            #print(current, collisions)

            current[0].handle_collision(current[1], current[2], delta)
            
            for index, collision in enumerate(collisions):

                if (current[0] in collision and isinstance(current[0], \
                    SimulationBall)) or (current[1] in collision and \
                    isinstance(current[1], SimulationBall)):

                    collisions.pop(index)

            collisions += self.get_all_collisions(delta, \
                *[b for b in [current[0], current[1]] \
                if isinstance(b, SimulationBall)])

        # Stage 3
        
        for ball in self.balls:

            ball.update_after_collision(delta)

    def get_all_collisions(self, delta, *collidants):

        collisions = []

        if len(collidants) == 0:

            collidants = self.balls
        
        all_pairs = [{a, b} for a in collidants for b in self.balls+self.blocks]

        while len(all_pairs) > 0:

            pair = all_pairs.pop()

            pair = list(pair)

            if pair not in all_pairs and len(pair) == 2:

                #print(pair)
                    
                collision = pair[0].collision_check(pair[1], delta)

                #print(collision)

                if collision is not None:

                    if len(collisions) == 0:

                        collisions.append((pair[0], \
                                pair[1], collision))

                    else:

                        for index, collide in enumerate(collisions):

                            if collide[2] > collision:

                                collisions = collisions[:index] + [(pair[0], \
                                    pair[1], collision)] + collisions[index:]
                                
                                break
        print(collisions)
        
        return collisions


    def on_touch_down(self, touch):

        self.update(1/60)

# Applications

class SimulationApp(App):

    '''
    Functions as the wrapper for the SimulationManager.
    '''

    def build(self):

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

    # sim.initialise(
    #     balls=((200, 300, 0, -600),
    #            (500, 300, -600, -3600)),
    #     blocks=((200, 200),
    #             (500, 200),
    #             (200, 500),
    #             (500, 500)),
    #     frame_advance=True)

    sim.initialise(
        balls=((200,300,600,0),),
        blocks=((350,300),),
        frame_advance=True)
    
    sim.run()
    