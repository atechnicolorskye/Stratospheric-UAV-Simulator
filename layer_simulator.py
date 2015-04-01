"""
layer_simulator.py
UAV Flight Planner

DESCRIPTION
-----------

Layer Simulator Module:
This module simulates the descent of a UAV through multiple atmopsheric layers. The UAV has a certain lift/drag ratio (= glide ratio).
Using code from the base_simulator module, displacement, heading and flight time are found layer-by-layer. It is assumed that wind blowing
steady in a certain direction and speed, air density is constant and the UAV is assumed to have been carried up by a weather balloon (hence
acquiring its speed and direction).

GLIDING FLIGHT PHYSICS
----------------------

When an aircraft is gliding, its lift/drag ratio is equivalent to its glide ratio (= forward speed/sink rate) in steady wind condtions.
As the UAV assumes a predefined speed in a direction predefinited by a wind field, the speed vector of the UAV can be said to be
equivalent to the forward speed compoent of the glide ratio which gives us a corresponding sink rate. In this simulator, we assume a 3D
wind field extending through the atmosphere that includes thermals to give rise to a component of wind speed in the z-direction. According to
aerodynamic convetion, z is postive downwards. The sink rate is variable in this case unlike the base_simulator module. Like in the
base_simulator module, we also consider the effects of constant headwinds and tailwinds in this simulation: headwinds and tailwinds reduces
and increases respectively the forward speed of the UAV. This leads to changes in the distance covered by the UAV, but not the time taken
for the UAV to reach the ground.

USAGE
-----

As this simulates descent through layers of atmospheric, a basic user interface is used to define the number of layers, the thickness of
the layers for the code of the base_simulator to work.  The factors considered will be wind speed (wind_xy (X & Y directions), wind_z (Z direction)),
wind direction (wind_dir), aircraft speed (uav_xy), aircraft direction (uav_dir) altitude (h), lift/drag ratio (LD). This simulator returns
displacement, heading, flight time layer-by-layer and the final displacement, heading and flight time.

University of Southampton
Si Kai Lee, skl2g11@soton.ac.uk
"""

__author__ = "Si Kai Lee, University of Southampton, skl2g11@soton.ac.uk"

import numpy
import pylab
import logging
import sys


class Layer(object):
    """
    Simulator of gliding UAV flight from a provided altitude, assuming wind field is constant.
    Returns [distance covered, flight time].
    """
    def __init__(self, wind_xy, wind_z, wind_dir, uav_xy, uav_dir, h, LD):

        # Defines and processes parameters
        # Directions are rounded to nearest integer and converted to between 0 and 359 degrees

        if wind_dir < 0:
            wind_dir = 360 - (abs(wind_dir) % 360)
        elif wind_dir > 360:
            wind_dir = wind_dir % 360
        else:
            wind_dir = wind_dir

        if uav_dir < 0:
            uav_dir = 360 - (abs(uav_dir) % 360)
        elif uav_dir > 360:
            uav_dir = uav_dir % 360
        else:
            uav_dir = uav_dir

         # Altitude must be positive
        if h <= 0:
            print "Layer Thickness must be greater than 0. Please try again.\n\n"
            LayerSimulator()

        # Wind speed must be postive
        if wind_xy < 0:
            print "2D Wind Speed must be positive. Please try again.\n\n"
            LayerSimulator()

        # UAV Speed must be positive
        if uav_xy <= 0:
            print "UAV Speed must be greater than 0. Please try again.\n\n"
            LayerSimulator()

        self.wind_xy = abs(float(wind_xy))
        self.wind_z = float(wind_z)
        self.wind_direction = round(wind_dir) * (2 * numpy.pi) / 360
        self.uav_speed = abs(float(uav_xy))
        self.uav_direction = round(uav_dir) * (2 * numpy.pi) / 360
        self.altitude = float(h)
        self.lift_drag_ratio = abs(float(LD))
        self.wind_U = self.wind_xy * numpy.cos(self.wind_direction)
        self.wind_V = self.wind_xy * numpy.sin(self.wind_direction)
        self.uav_U = self.uav_speed * numpy.cos(self.uav_direction)
        self.uav_V = self.uav_speed * numpy.sin(self.uav_direction)

    def LayerFly(self):
        """
        The basic flght simulator, takes wind speed in u and v components, UAV speed in u and v components, altitude, lift/drag ratio.
        """

        # Calculate sink rate (sink_rate) and relative groundspeed (UV) in u (U) and v (V) components to find new glide ratio (glide_ratio).
        sink_rate = self.uav_speed / self.lift_drag_ratio - self.wind_z
        U = self.uav_U + self.wind_U
        V = self.uav_V + self.wind_V
        UV = (U ** 2 + V ** 2) ** 0.5
        glide_ratio = UV / sink_rate

        # From final lift/drag ratio (= glide ratio), find displacement (displ), heading (head), flight time (time)
        displ = glide_ratio * self.altitude
        # Use arctan2 to ensure no errors with angles
        if abs(U) < 1e-4:
            U = 0
        if abs(V) < 1e-4:
            V = 0
        head = numpy.arctan2(V, U)
        time = self.altitude / sink_rate

        return [displ, head, time]


def LayerSimulator():
    print "This is the simulator for UAV flights across severals layers of atmosphere. Please follow the instructions below."
    print "Only the number of atmospheric layers need to be an integer; the rest of your inputs can be floating point numbers. "
    layers = raw_input("Enter the number of atmospheric layers: ")
    # Checks if layers is an integer
    try:
        layers = int(layers)
    except ValueError:
        print "This is not an integer. Please try again.\n\n"
        LayerSimulator()
    # Creates two empty lists to store inputs and outputs
    # Layer Input
    l_i = []
    # Layer Output
    l_o = []

    for x in range(1, layers + 1, 1):
        if x == 1:
            print "\n"
            print "Layer {}".format(x)
            print "Please input a list with the following variables in order separated by commas."
            print "2D Wind Speed, Wind Speed in Z direction, Wind Direction, 2D UAV Speed,"
            print "UAV Heading, Layer Thickness, Lift to Drag Ratio."
            var = raw_input(": ")
            try:
                l_i.append([float(y) for y in var.split(',')])
                if l_i[0][1] > float(l_i[0][3]) / l_i[0][6]:
                    print 'Wind Speed in Z Direction must be less than UAV Sink Rate. Please try again.\n\n'
                    LayerSimulator()
                else:
                    y = Layer(l_i[0][0], l_i[0][1], l_i[0][2], l_i[0][3], l_i[0][4], l_i[0][5], l_i[0][6])
                    l_o.append(y.LayerFly())
                    print "\n"
            except:
                print "There was an input error. Please try again.\n\n"
                LayerSimulator()
        else:
            print "Layer {}".format(x)
            print "Please input a list with the following variables in order separated by commas:"
            print "2D Wind Speed, Wind Speed in Z direction, Wind Direction, Layer Thickness"
            var = raw_input(": ")
            try:
                l_i.append([float(y) for y in var.split(',')])
                x_ = int(x) - 1
                if l_i[x_][1] > float(l_i[0][3]) / l_i[0][6]:
                    print 'Wind Speed in Z Direction must be less than UAV Sink Rate. Please try again.\n\n'
                    LayerSimulator()
                else:
                    y = Layer(l_i[x_][0], l_i[x_][1], l_i[x_][2], l_i[0][3], l_i[0][4], l_i[x_][3], l_i[0][6])
                    l_o.append(y.LayerFly())
                    print "\n"
            except:
                print "There was an input error. Please try again.\n\n"
                LayerSimulator()

    # Creates nested lists to store x and y coordinates of the UAV
    coords = [[0, 0]]
    descent_time = 0

    # Creates counter for moving x and y coordinates in coords
    count = 0

    for x in l_o:
        if count == 0:
            coords.append([x[0] * numpy.cos(x[1]), x[0] * numpy.sin(x[1])])
            count += 1
        else:
            coords.append([coords[count][0] + x[0] * numpy.cos(x[1]), coords[count][1] + x[0] * numpy.sin(x[1])])
            count += 1
        descent_time = descent_time + x[2]

    # arctan2(y, x)
    if abs(coords[-1][1]) < 1e-4:
        coords[-1][1] = 0
    if abs(coords[-1][0]) < 1e-4:
        coords[-1][0] = 0
    total_head = numpy.arctan2(coords[-1][1], coords[-1][0]) * 360 / (2 * numpy.pi)
    if total_head < 0:
        total_head = 360 + total_head
    else:
        total_head = total_head

    xy = (coords[-1][0] ** 2 + coords[-1][1] ** 2) ** 0.5

    fig = pylab.figure()
    # Adjusts margins
    fig.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)

    # Plots points
    for x in range(0, layers + 1, 1):
        # Plots start points
        if x < 1:
            pylab.plot(coords[x][0], coords[x][1], 'bo', markersize=10, label='Start')
        # Plots intermediate points
        elif x < layers:
            pylab.plot(coords[x][0], coords[x][1], 'ro', markersize=10, label='Displacement Point {}'.format(x))
        # Plots final points
        else:
            pylab.plot(coords[x][0], coords[x][1], 'ko', markersize=10, label='Final Displacement')

    TextL = 'Total Displacement (m):\nFinal Heading (Degrees):\nFlight Time (s):'
    TextR = str(int(xy)) + '\n' + str(total_head) + '\n' + str(descent_time)
    fig.text(0.25, 0.01, TextL, multialignment='left')
    fig.text(0.73, 0.01, TextR, multialignment='right')

    pylab.title('UAV Displacement in Y against UAV Displacement in X')
    pylab.xlabel('UAV Displacement in X')
    pylab.ylabel('UAV Displacement in Y')
    #pylab.axis('equal')
    #pylab.axis([-xy-100, xy+100, -xy-100, xy+100])
    pylab.xlim(-xy - 100, xy + 100)
    pylab.ylim(-xy - 100, xy + 100)
    pylab.grid(True)
    pylab.axes().set_aspect('equal')
    pylab.legend()

    pylab.show()

    return [xy, total_head, descent_time]
