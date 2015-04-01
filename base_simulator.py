"""
base_simulator.py
UAV Flight Planner

DESCRIPTION
-----------

Base Simulator Module:
This module simulates the descent of a UAV with a certain lift/drag ratio (= glide ratio), finds displacement, heading and flight time from a
certain predefined altitude with 2D wind blowing steady in a certain direction and speed. The UAV is assumed to have been carried up by a
weather balloon. The UAV has only speed in the direction dictated by the wind field. The wind field is assumed to be constant
throughout all layers of the atmosphere and air density is constant.

GLIDING FLIGHT PHYSICS
----------------------

When an aircraft is gliding, its lift/drag ratio is equivalent to its glide ratio (= forward speed/sink rate) in steady wind condtions.
As the UAV assumes a predefined speed in a direction predefinited by a wind field, the speed vector of the UAV can be said to be
equivalent to the forward speed compoent of the glide ratio. Hence this gives us a corresponding sink rate. In this basic simulator, we
assume a 2D wind field extending through the atmosphere, there are no thermals to give rise to an upward component of speed. However, we
do consider the effects of constant headwinds and tailwinds in this simulation: headwinds and tailwinds reduces and increases respectively
the forward speed of the UAV. This leads to changes in the distance covered by the UAV, but not the time taken for the UAV to reach the
ground. The sink rate is assumed constant as the presence of vertical winds (updrafts and downdrafts) is ignored.


USAGE
-----

As this is the most basic simulator module of the project, the only factors considered will be wind speed (wind_v), wind direction (wind_dir),
aircraft speed (uav_v), aircraft direction (uav_dir) altitude (h), lift/drag ratio (LD). This simulator returns [displacement, heading, flight time].

University of Southampton
Si Kai Lee, skl2g11@soton.ac.uk
"""

__author__ = "Si Kai Lee, University of Southampton, skl2g11@soton.ac.uk"

import numpy
import pylab


class BaseSimulator(object):
    """
    Simulator of gliding UAV flight from a provided altitude, assuming wind field is constant.
    Returns [distance covered, flight time].
    """
    def __init__(self, wind_v, wind_dir, uav_v, uav_dir, h, LD, marker1='bo', marker2='ro', marker3='r.'):

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
            raise ValueError("Altitude must be greater than 0.")

        # Wind speed must be postive
        if wind_v < 0:
            raise ValueError("Wind speed must be positive.")

        # UAV Speed must be positive
        if uav_v <= 0:
            raise ValueError("UAV Speed must be greater than 0.")

        self.wind_speed = abs(float(wind_v))
        self.wind_direction = numpy.radians(round(wind_dir))
        self.uav_speed = abs(float(uav_v))
        self.uav_direction = round(uav_dir) * (2 * numpy.pi) / 360
        self.altitude = float(h)
        self.lift_drag_ratio = abs(float(LD))
        self.wind_U = self.wind_speed * numpy.cos(self.wind_direction)
        self.wind_V = self.wind_speed * numpy.sin(self.wind_direction)
        self.uav_U = self.uav_speed * numpy.cos(self.uav_direction)
        self.uav_V = self.uav_speed * numpy.sin(self.uav_direction)
        self.marker1 = marker1
        self.marker2 = marker2
        self.marker3 = marker3

    def BaseFly(self):
        """
        The basic flght simulator, takes wind speed in u and v components, UAV speed in u and v components, altitude, lift/drag ratio.
        """

        # Calculate sink rate (sink_rate) and relative groundspeed (UV) in u (U) and v (V) components to find new glide ratio (glide_ratio).
        sink_rate = self.uav_speed / self.lift_drag_ratio

        if abs(self.wind_U) > 1e-4:
            U = self.uav_U + self.wind_U
        else:
            U = self.uav_U

        if abs(self.wind_V) > 1e-4:
            V = self.uav_V + self.wind_V
        else:
            V = self.uav_V

        UV = (U ** 2 + V ** 2) ** 0.5

        glide_ratio = UV / sink_rate

        # From final lift/drag ratio (= glide ratio), find displacement (displ), heading (head), flight time (time)
        displ = glide_ratio * self.altitude
        # Use arctan2 to ensure no errors with angles
        if abs(U) < 1e-4:
            U = 0
        if abs(V) < 1e-4:
            V = 0
        head = numpy.arctan2(V, U) * 360 / (2 * numpy.pi)
        if head < 0:
            head = 360 + head
        else:
            head = head
        time = self.altitude / sink_rate

        # Sets variables to allow plotting
        self.displ = displ
        self.head = head
        self.time = time

        return [displ, head, time]

    def PlotBaseFly(self):
        self.BaseFly()
        # Creates lists of x and y coordinates
        x_coords = [0, (self.displ * numpy.cos(self.head * (2 * numpy.pi) / 360))]
        y_coords = [0, (self.displ * numpy.sin(self.head * (2 * numpy.pi) / 360))]

        fig = pylab.figure()
        # Adjusts margins
        fig.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)

        # Plots start point
        pylab.plot(x_coords[0], y_coords[0], self.marker1, markersize=10, label='Start')
        # Plots end point
        pylab.plot(x_coords[1], y_coords[1], self.marker2, markersize=10, label='End')

        # Generates 2D wind vector field
        if self.displ < 10:
            X, Y = numpy.meshgrid(numpy.arange(-30, 30, 6), numpy.arange(-30, 30, 6))
            Q = pylab.quiver(X, Y, self.wind_U, self.wind_V, label='Wind Direction')
            QK = pylab.quiverkey(Q, 0.1, -0.08, 10, 'Wind Direction', labelpos='W',
                                 fontproperties={'weight': 'bold'})
        else:
            X, Y = numpy.meshgrid(numpy.arange(-self.displ - 20, self.displ + 20, self.displ / 5), numpy.arange(-self.displ - 20, self.displ + 20, self.displ / 5))
            Q = pylab.quiver(X, Y, self.wind_U, self.wind_V, label='Wind Direction')
            QK = pylab.quiverkey(Q, 0.1, -0.08, self.displ / 5, 'Wind Direction', labelpos='W',
                                 fontproperties={'weight': 'bold'})

        # Shows displ, head, time
        # Refer to http://docs.python.org/2/library/string.html for formatting
        # and http://stackoverflow.com/questions/19261089/how-to-use-new-style-string-formatting-in-matplotlib-figure-text
        '''
        Text = '{0} {1:10}\n'.format('Displacement: ', int(self.displ)) + \
               '{0} {1:10}\n'.format('Heading:         ', int(self.head)) + \
               '{0} {1:10}\n'.format('Flight Time:  ', int(self.time))

        fig.text(0.25, -0.01, Text, family='monospace')
        '''

        TextL = 'Displacement (m):\nHeading (Degrees):\nDescent Time (s):'
        TextR = str(int(self.displ)) + '\n' + str(int(self.head)) + '\n' + str(int(self.time))
        fig.text(0.4, 0.01, TextL, multialignment='left')
        fig.text(0.6, 0.01, TextR, multialignment='right')

        pylab.title('UAV Displacement in Y against UAV Displacement in X')
        pylab.xlabel('UAV Displacement in X')
        pylab.ylabel('UAV Displacement in Y')
        pylab.xlim(-self.displ - 20, self.displ + 20)
        pylab.ylim(-self.displ - 20, self.displ + 20)

        pylab.legend()

        pylab.show()

    def PlotWindContour(self):
        """
        Plots the contour of the distance travelled with a wind blowing from 0 to 359 degrees. Wind speed, UAV speed and UAV direction is
        set using the BaseSimulator class. BaseFly is used to find displacement.
        """

        # Creates list for wind direction from 0 to 359
        wind_direction_0_359 = range(0, 360, 1)

        # Creates empty list to store distance travelled by UAV
        displacement = []

        # Creates empty lists to store x and y coordinates
        x_coords = []
        y_coords = []

        # For-loop to run through the possible wind directions
        for x in wind_direction_0_359:
            self.wind_U = self.wind_speed * numpy.cos(x * (2 * numpy.pi) / 360)
            self.wind_V = self.wind_speed * numpy.sin(x * (2 * numpy.pi) / 360)
            y = self.BaseFly()
            # Appends results of BaseFly into relevant lists
            displacement.append(y[0])
            x_coords.append(y[0] * numpy.cos(y[1] * (2 * numpy.pi) / 360))
            y_coords.append(y[0] * numpy.sin(y[1] * (2 * numpy.pi) / 360))

        pylab.plot(x_coords, y_coords, self.marker3, label='Displacement Contour')
        pylab.plot(0, 0, self.marker1, label="Start")

        pylab.title('UAV Displacement in Y against UAV Displacement in X\n(Varying Wind Direction from 0 to 359)')
        pylab.xlabel('UAV Displacement in X')
        pylab.ylabel('UAV Displacement in Y')
        pylab.legend()

        pylab.show()

    def PlotFlightContour(self):

        # Creates list for UAV directions from 0 to 359
        uav_direction_0_359 = numpy.arange(0, 360, 0.1)

        # Creates empty list to store distance travelled by UAV
        displacement = []

        # Creates empty lists to store x and y coordinates
        x_coords = []
        y_coords = []

        # For-loop to run through the possible UAV directions
        for x in uav_direction_0_359:
            self.uav_U = self.uav_speed * numpy.cos(x * (2 * numpy.pi) / 360)
            self.uav_V = self.uav_speed * numpy.sin(x * (2 * numpy.pi) / 360)
            y = self.BaseFly()
            # Appends results of BaseFly into relevant lists
            displacement.append(y[0])
            x_coords.append(y[0] * numpy.cos(y[1] * (2 * numpy.pi) / 360))
            y_coords.append(y[0] * numpy.sin(y[1] * (2 * numpy.pi) / 360))

        # Generates 2D wind vector field
        X, Y = numpy.meshgrid(numpy.arange(min(x_coords) - 80, max(x_coords) + 80, (max(x_coords) - min(x_coords)) / 5),
                              numpy.arange(min(y_coords) - 80, max(y_coords) + 80, (max(y_coords) - min(y_coords)) / 5))
        Q = pylab.quiver(X, Y, self.wind_U, self.wind_V, label='Wind Direction')
        QK = pylab.quiverkey(Q, 0.1, -0.08, 10, 'Wind Direction', labelpos='W',
                             fontproperties={'weight': 'bold'})

        pylab.plot(x_coords, y_coords, self.marker3, label='Displacement Contour')
        pylab.plot(0, 0, self.marker1, label="Start")
        pylab.title('UAV Displacement in Y against UAV Displacement in X\n(Varying UAV Direction from 0 to 359)')
        pylab.xlabel('UAV Displacement in X')
        pylab.ylabel('UAV Displacement in Y')

        pylab.legend()

        pylab.show()
