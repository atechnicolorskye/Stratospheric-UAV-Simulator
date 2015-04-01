"""
gfs_data_simulator_local.py
GFS Local Data Simulator

DESCRIPTION
-----------

GFS Data Local Simulator Module:
This module simulates the descent of a UAV through multiple atmopsheric layers. The properties of the layers are obtained from National
Oceanic and Atmospheric Administration's  (NOAA) Global Forecast System (GFS) using Niccolo' Zapponi's GFS and weather modules. The module
stores the initial Weather object to improve processing times by reducing the number of downloads. Using code from the base_simulator module,
displacement, heading and flight time are found layer-by-layer. The location and other relevant data of the UAV are updated every x seconds
of flight until it is above 100m where updates are carried out every y seconds. The UAV is assumed to have been carried up by a weather balloon
and there are no updrafts and downdrafts.

GLIDING FLIGHT PHYSICS
----------------------

When an aircraft is gliding, its lift/drag ratio is equivalent to its glide ratio (= forward velocity/sink rate) in steady wind condtions.
As the UAV assumes a predefined velocity in a direction predefinited by a wind field, the velocity vector of the UAV can be said to be
equivalent to the forward velocity compoent of the glide ratio which gives us a corresponding sink rate. In this simulator, we assume a 2D
wind field extending through the atmosphere. Like in the base_simulator module, we also consider the effects of constant headwinds and tailwinds
in this simulation: headwinds and tailwinds reduces and increases respectively the forward velocity of the UAV. This leads to changes in the
distance covered by the UAV, but not the time taken for the UAV to reach the ground.

USAGE
-----

The code requires the user to create an instance of the Flight() object. The user has the choice to call the methods Fly_1, Fly_2, Fly_3 and
Fly_4. Fly_1, Fly_2 and Fly_3 serves to simulate the descent of the UAV while keeping to a certain heading using either a distance-based
or a heading-based algorithm. Fly_4 demonstrates the UAV's capability to reach various different waypoints that are user determined using
the distance-based algorithm from the moethod Fly from gfs_data_simulator. Once one of the above methods hs been called, the user can call
the method PlotContour that plots relevant figures to visualise the descent process. Fly_Range is used to allow the user to plot the descent
of the UAV for various set headings and is used for Fly_1, Fly_2 and Fly_3 specifically. Fly_Range_2 is used to plot the results from Fly_4 for
various pre-determined headings.

Exmaple:
y = Flight(50.1, -5.00, 30000, (2014,03,01,19,31,01), 10, 2, 0.6, 0.4, 3, 0.5)
y.Fly_1(90) # Same for Fly_2 and Fly_3
y.PlotContour()

y = Flight(50.1, -5.00, 30000, (2014,03,01,19,31,01), 10, 2, 0.6, 0.4, 3, 0.5)
y.Fly_4(5000, 5000, 500, 45)

y = Flight(50.1, -5.00, 7500, (2014,03,05,19,31,01), 5, 2, 0.2, 0.15, 3, 0.5)
y.Fly_Range(3, [45,90], False, False, True, True, False)

y = Flight(50.1, -5.00, 5000, (2014,03,10,19,31,01), 5, 2, 0.2, 0.15, 3, 0.5)
y.Fly_Range_2(0, 359, 10000, 100000, 5000, 5)

University of Southampton
Si Kai Lee, skl2g11@soton.ac.uk
"""

__author__ = "Si Kai Lee, University of Southampton, skl2g11@soton.ac.uk"


from global_tools import m2deg, getUTCOffset
from datetime import datetime, timedelta
import weather

import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from scipy.integrate import odeint
from write_file import WriteToFile_Dist, WriteToFile_Head1, WriteToFile_Head2, WriteToFile_Local, WriteToFile_Simulator

import logging
from sys import exit


class Flight(object):
    """
    Class that contains the simulator of the simulation of UAV flight from a provided altitude.
    Instantiates the Flight_LatLon object when called.
    Inputs: Starting Latitutde, Starting Longitude, Starting Altitude, Starting Time, Average Lift-to-Drag Ratio,
            Mass of UAV, PlanformArea (Assuming UAV is a flying wind),Average Coefficient of Lift, Time Step used
            for ODEs at altitudes above 100m, Time Step used for ODEs at altitudes below 100m.
    """
    def __init__(self, Lat, Lon, Alt, StartTime, LD, Mass, PlanformArea, CL, TimeStep1, TimeStep2):
        # Initialises flight container object

        # Sets up logging configuration
        logging.basicConfig(filename='Fly.log', format='%(levelname)s:%(message)s', level=logging.DEBUG)

        # User defined variables
        self.Lat = Lat
        self.Lon = Lon
        self.Alt = Alt
        self._StartTime = StartTime
        self.StartTime = datetime(*StartTime)
        self.LD = LD
        self.Mass = Mass
        self.PlanformArea = PlanformArea
        self.CL = CL
        self.TimeStep1 = TimeStep1
        self.TimeStep2 = TimeStep1

        # Creates empty lists to store x and y coordinates and defined heading for UAV for Fly_4 and Fly_Range_2
        self.x_coords = []
        self.y_coords = []
        self.x_coords_max = []
        self.y_coords_max = []
        self.x_coords_end = []
        self.y_coords_end = []
        self.set_heading = []
        self.end_heading = []
        self.x_coords_fail = []
        self.y_coords_fail = []
        self.x_coords_end_fail = []
        self.y_coords_end_fail = []
        self.set_heading_fail = []
        self.end_heading_fail = []

        # Checks if StartTime is between 30 days before and 6 after current time
        if self.StartTime > datetime.now():
            if self.StartTime - datetime.now() > timedelta(days=6):
                raise Exception('Your selected date and time is greater than 6 days in the future. Please try again.')
                exit()
        else:
            if datetime.now() - self.StartTime > timedelta(days=30):
                raise Exception('Your selected date and time is greater than 30 days in the past. Please try again.')
                exit()

        # Dependent variables
        self._BaseLat = self.Lat
        self._BaseLon = self.Lon
        self._xyWindSpeed = 0
        self._WindHead = 0
        self._GlideAngle = numpy.arctan(float(1) / self.LD)
        self._XZCurrentSpeed = 0
        self._XZTempSpeed = 0
        self._CurrentDensity = 0
        self._CurrentLatLon = 0
        self._CurrentLat = self.Lat
        self._CurrentLon = self.Lon
        self._CurrentAlt = self.Alt
        self._CurrentDist = 0
        self._CurrentTime = self.StartTime
        self._TimeDump = [0]
        self._AltDump = [self.Alt]
        self._DistDump = [0]
        self._DensityDump = []
        self._XDirDump = [0]
        self._YDirDump = [0]
        self._HeadingDump = []
        self._XWindDump = []
        self._YWindDump = []

        # Intialises Weather Environment object
        self.WeatherData = [weather.forecastEnvironment()]

        # Obtains data from GFS for current flight conditions
        self.WeatherData[0].launchSiteLat = self.Lat
        self.WeatherData[0].launchSiteLon = self.Lon
        self.WeatherData[0].launchSiteElev = self.Alt
        self.WeatherData[0].dateAndTime = self.StartTime
        self.WeatherData[0].UTC_offset = getUTCOffset(self.Lat, self.Lon, self.StartTime)

        # Downloads GFS weather data
        print 'Downloading the forecast (might take a while)...'
        self.WeatherData[0].loadForecast()
        logging.info('GFS Data downloaded!')
        print "Forecast downloaded!"

        # Get relevant properties of current air mass
        self._xyWindSpeed = self.WeatherData[0].getWindSpeed(self.Lat, self.Lon, self.Alt, self.StartTime)
        logging.info('XY Wind Speed: %s', self._xyWindSpeed)

        self._WindHead = self.WeatherData[0].getWindDirection(self.Lat, self.Lon, self.Alt, self.StartTime)
        logging.info('Wind Heading: %s', self._WindHead)

        self._CurrentDensity = self.WeatherData[0].getDensity(self.Lat, self.Lon, self.Alt, self.StartTime)
        logging.info('Current Density: %s', self._CurrentDensity)

        # Appends Density in data dumps
        self._DensityDump.append(self._CurrentDensity)

    def Clear(self):
        # Resets object attributes

        self.WeatherData = [self.WeatherData[0]]

        # Dependent variables
        self._BaseLat = self.Lat
        self._BaseLon = self.Lon
        self._xyWindSpeed = self.WeatherData[0].getWindSpeed(self.Lat, self.Lon, self.Alt, self.StartTime)
        self._WindHead = self.WeatherData[0].getWindDirection(self.Lat, self.Lon, self.Alt, self.StartTime)
        self._GlideAngle = numpy.arctan(float(1) / self.LD)
        self._XZCurrentSpeed = 0
        self._XZTempSpeed = 0
        self._CurrentDensity = self.WeatherData[0].getDensity(self.Lat, self.Lon, self.Alt, self.StartTime)
        self._CurrentLatLon = 0
        self._CurrentLat = self.Lat
        self._CurrentLon = self.Lon
        self._CurrentAlt = self.Alt
        self._CurrentDist = 0
        self._CurrentTime = self.StartTime
        self._TimeDump = [0]
        self._AltDump = [self.Alt]
        self._DistDump = [0]
        self._DensityDump = [self._CurrentDensity]
        self._XDirDump = [0]
        self._YDirDump = [0]
        self._HeadingDump = []
        self._XWindDump = []
        self._YWindDump = []

    def Fly_1(self, FinalHead):
        """
        Simulates of gliding UAV flight from a provided altitude; weather data is provided by GFS.
        Uses distance-based algorithm to calculate UAV heading.
        Inputs: Desired Heading
        Prints and Returns: End Latitude, End Longtitude, Ideal End Latitude, Ideal End Longitude, End Heading, Desired Heading
        """
        # Fly for TimeStep1 seconds and recalculates dependent variables

        print 'Running Fly_1'
        # Converts FinalHead from deg to rad
        self.FinalHead = numpy.deg2rad(FinalHead)

        logging.info('User Input:\nStart Latitude: %s\nStart Longtitude: %s\nStart Altitude: %s\nDesired Heading: %s\nStart Time: %s\nLift to Drag Ratio: %s\n \
                     Mass: %s\nWing Planform Area: %s\nCoefficient of Lift: %s', self.Lat, self.Lon, self.Alt, self.FinalHead, self.StartTime, self.LD,
                     self.Mass, self.PlanformArea, self.CL)

        while self._CurrentAlt > 100:
            # Fly for TimeStep1 seconds and recalculates dependent variables

            # Calculates heading of UAV
            self._BaseXDiff = self._DistDump[-1] * numpy.cos(self.FinalHead)
            self._BaseYDiff = self._DistDump[-1] * numpy.sin(self.FinalHead)

            self._XDiff, self._YDiff = self._XDirDump[-1], self._YDirDump[-1]

            if (self._BaseXDiff - self._XDiff) == 0 and (self._BaseYDiff - self._YDiff) > 0:
                self._Heading = numpy.pi / 2
            elif (self._BaseXDiff - self._XDiff) < 0 and (self._BaseYDiff - self._YDiff) == 0:
                self._Heading = numpy.pi
            elif (self._BaseXDiff - self._XDiff) == 0 and (self._BaseYDiff - self._YDiff) < 0:
                self._Heading = 3 * numpy.pi / 2
            else:
                self._Heading = numpy.arctan2(self._BaseYDiff - self._YDiff, self._BaseXDiff - self._XDiff)
            self._HeadingDump.append(self._Heading)

            # Provides list of time values for ODE solver
            time = numpy.linspace(0, self.TimeStep1)

            # Variables required for ODE solver
            SinA = numpy.sin(self._GlideAngle)
            CosA = numpy.cos(self._GlideAngle)
            B = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass)
            C = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass * self.LD)

            # Sets up ODEs to be solved
            def z_derivatives(z, t):
                return numpy.array([-z[1] * SinA, - B * (z[1] ** 2) + 9.81 * SinA])

            def xy_derivatives(xy, t):
                return numpy.array([xy[1] * CosA, - C * (xy[1] ** 2) + 9.81 * CosA])

            # Solves for change in altitude
            z_initial_conditions = numpy.array([self._CurrentAlt, self._XZCurrentSpeed])
            z = odeint(z_derivatives, z_initial_conditions, time)
            self._CurrentAlt, self._XZCurrentSpeed = z[-1]

            # Solves for distance travelled
            xy_initial_conditions = numpy.array([self._CurrentDist, self._XZTempSpeed])
            xy = odeint(xy_derivatives, xy_initial_conditions, time)
            self._CurrentDist, self._XZTempSpeed = xy[-1]

            # Appends changes in X and Y directions to X and Y data dumps
            if self._Heading == 0:
                self._XDirDump.append(self._XDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            elif self._Heading == numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            elif self._Heading == numpy.pi:
                self._XDirDump.append(self._XDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            elif self._Heading == 3 * numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            else:
                self._XDirDump.append(self._XDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.cos(self._Heading)) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.sin(self._Heading)) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)

            # Appends the wind conditions to wind data dumps
            self._XWindDump.append(self._xyWindSpeed * numpy.cos(self._WindHead))
            self._YWindDump.append(self._xyWindSpeed * numpy.sin(self._WindHead))

            # Updates parameters used to obtain the next set of forecasts
            self._CurrentLatLon = m2deg(self._YDirDump[-1], self._XDirDump[-1], self.Lat)
            self._CurrentLat = self.Lat + self._CurrentLatLon[0]
            self._CurrentLon = self.Lon + self._CurrentLatLon[1]
            self._CurrentTime += timedelta(seconds=self.TimeStep1)

            # Appends useful data to data dumps
            self._TimeDump.append(self._TimeDump[-1] + self.TimeStep1)
            self._AltDump.append(self._CurrentAlt)
            self._DistDump.append(self._CurrentDist)

            logging.info('Current Latitude: %s\nCurrent Longitude: %s\nCurrent Altitude: %s\nCurrent Time: %s', self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)

            # Check if region travelled exceeds 3 degrees Lat and 3 degrees Lon from the launch site
            if abs(self._CurrentLat - self._BaseLat) < 1.5 and abs(self._CurrentLon - self._BaseLon) < 1.5:
                # Retrieves changing variables from forecast
                self._xyWindSpeed = self.WeatherData[-1].getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = self.WeatherData[-1].getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = self.WeatherData[-1].getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)
            else:
                print 'UAV has travelled out of forecast region.'
                print 'New forecast has to be obtained, please wait while forecast is downloaded and processed.'
                self._BaseLat = self._CurrentLat
                self._BaseLon = self._CurrentLon

                # Intialises Weather Environment object
                self.WeatherData.append(weather.forecastEnvironment())

                # Obtains data from GFS for current flight conditions
                self.WeatherData[-1].launchSiteLat = self._BaseLat
                self.WeatherData[-1].launchSiteLon = self._BaseLon
                self.WeatherData[-1].launchSiteElev = self._CurrentAlt
                self.WeatherData[-1].dateAndTime = self._CurrentTime
                self.WeatherData[-1].UTC_offset = getUTCOffset(self._BaseLat, self._BaseLon, self._CurrentTime)

                # Downloads GFS weather data
                self.WeatherData[-1].loadForecast()
                logging.info('GFS Data downloaded!')
                print "Forecast downloaded!"

                # Retrieves changing variables from forecast
                self._xyWindSpeed = self.WeatherData[-1].getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = self.WeatherData[-1].getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = self.WeatherData[-1].getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)

        while self._CurrentAlt > 0:
            # Fly for TimeStep2 seconds and recalculates dependent variables

            self._BaseXDiff = self._DistDump[-1] * numpy.cos(self.FinalHead)
            self._BaseYDiff = self._DistDump[-1] * numpy.sin(self.FinalHead)

            self._XDiff, self._YDiff = self._XDirDump[-1], self._YDirDump[-1]

            if (self._BaseXDiff - self._XDiff) == 0 and (self._BaseYDiff - self._YDiff) > 0:
                self._Heading = numpy.pi / 2
            elif (self._BaseXDiff - self._XDiff) < 0 and (self._BaseYDiff - self._YDiff) == 0:
                self._Heading = numpy.pi
            elif (self._BaseXDiff - self._XDiff) == 0 and (self._BaseYDiff - self._YDiff) < 0:
                self._Heading = 3 * numpy.pi / 2
            else:
                self._Heading = numpy.arctan2(self._BaseYDiff - self._YDiff, self._BaseXDiff - self._XDiff)
            self._HeadingDump.append(self._Heading)

            # Provides list of time values for ODE solver
            time = numpy.linspace(0, self.TimeStep2, num=100)

            # Variables required for ODE solver
            SinA = numpy.sin(self._GlideAngle)
            CosA = numpy.cos(self._GlideAngle)
            B = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass)
            C = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass * self.LD)

            # Sets up ODEs to be solved
            def z_derivatives(z, t):
                return numpy.array([-z[1] * SinA, - B * (z[1] ** 2) + 9.81 * SinA])

            def xy_derivatives(xy, t):
                return numpy.array([xy[1] * CosA, - C * (xy[1] ** 2) + 9.81 * CosA])

            # Solves for change in altitude
            z_initial_conditions = numpy.array([self._CurrentAlt, self._XZCurrentSpeed])
            z = odeint(z_derivatives, z_initial_conditions, time)
            self._CurrentAlt, self._XZCurrentSpeed = z[-1]

            # Solves for distance travelled
            xy_initial_conditions = numpy.array([self._CurrentDist, self._XZTempSpeed])
            xy = odeint(xy_derivatives, xy_initial_conditions, time)
            self._CurrentDist, self._XZTempSpeed = xy[-1]

            # Appends changes in X and Y directions to X and Y data dumps
            if self._Heading == 0:
                self._XDirDump.append(self._XDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)
            elif self._Heading == numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)
            elif self._Heading == numpy.pi:
                self._XDirDump.append(self._XDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)
            elif self._Heading == 3 * numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)
            else:
                self._XDirDump.append(self._XDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.cos(self._Heading)) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.sin(self._Heading)) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)

            # Appends the wind conditions to wind data dumps
            self._XWindDump.append(self._xyWindSpeed * numpy.cos(self._WindHead))
            self._YWindDump.append(self._xyWindSpeed * numpy.sin(self._WindHead))

            # Updates parameters used to obtain the next set of forecasts
            self._CurrentLatLon = m2deg(self._YDirDump[-1], self._XDirDump[-1], self.Lat)
            self._CurrentLat = self.Lat + self._CurrentLatLon[0]
            self._CurrentLon = self.Lon + self._CurrentLatLon[1]
            self._CurrentTime += timedelta(seconds=self.TimeStep2)

            # Appends useful data to data dumps
            self._TimeDump.append(self._TimeDump[-1] + self.TimeStep2)
            self._AltDump.append(self._CurrentAlt)
            self._DistDump.append(self._CurrentDist)

            logging.info('Current Latitude: %s\nCurrent Longitude: %s\nCurrent Altitude: %s\nCurrent Time: %s', self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)

            # Check if region travelled exceeds 3 degrees Lat and 3 degrees Lon from the launch site
            if abs(self._CurrentLat - self._BaseLat) < 1.5 and abs(self._CurrentLon - self._BaseLon) < 1.5:
                # Retrieves changing variables from forecast
                self._xyWindSpeed = self.WeatherData[-1].getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = self.WeatherData[-1].getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = self.WeatherData[-1].getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)

                # Checks if UAV is low enough for calculations to stop
                if self._CurrentAlt < 1:
                    break
            else:
                print 'UAV has travelled out of forecast region.'
                print 'New forecast has to be obtained, please wait while forecast is downloaded and processed.'
                self._BaseLat = self._CurrentLat
                self._BaseLon = self._CurrentLon

                # Intialises Weather Environment object
                self.WeatherData.append(weather.forecastEnvironment())

                # Obtains data from GFS for current flight conditions
                self.WeatherData[-1].launchSiteLat = self._BaseLat
                self.WeatherData[-1].launchSiteLon = self._BaseLon
                self.WeatherData[-1].launchSiteElev = self._CurrentAlt
                self.WeatherData[-1].dateAndTime = self._CurrentTime
                self.WeatherData[-1].UTC_offset = getUTCOffset(self._BaseLat, self._BaseLon, self._CurrentTime)

                # Downloads GFS weather data
                self.WeatherData[-1].loadForecast()
                logging.info('GFS Data downloaded!')
                print "Forecast downloaded!"

                # Retrieves changing variables from forecast
                self._xyWindSpeed = self.WeatherData[-1].getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = self.WeatherData[-1].getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = self.WeatherData[-1].getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)

                # Checks if UAV is low enough for calculations to stop
                if self._CurrentAlt < 1:
                    break

        # Prepares output
        self._FinalLatLon = m2deg(self._YDirDump[-1], self._XDirDump[-1], self.Lat)
        self._EndHead = numpy.rad2deg(numpy.arctan2(self._YDirDump[-1], self._XDirDump[-1]))
        if self._EndHead < 0 and FinalHead > 90:
            self._EndHead += 360

        logging.info('Final Latitude: %s\nFinal Longitude: %s\nDescent Time: %s', self._FinalLatLon[0] + self.Lat, self._FinalLatLon[1] + self.Lon, (self._CurrentTime - self.StartTime))
        print 'Final Latitude: {}\nFinal Longitude: {}\nDescent Time: {}'.format(self._FinalLatLon[0] + self.Lat, self._FinalLatLon[1] + self.Lon, (self._CurrentTime - self.StartTime))
        print 'Final Heading: {}\nDesired Heading: {}\nFinal Altitude: {}'.format(self._EndHead, FinalHead, self._CurrentAlt)

        # Sets variables for plotting

        # Outputs
        return [self._XDirDump[-1], self._YDirDump[-1], self._CurrentDist * numpy.cos(self.FinalHead), self._CurrentDist * numpy.sin(self.FinalHead), FinalHead, self._EndHead]

    def Fly_2(self, FinalHead):
        """
        Simulates of gliding UAV flight from a provided altitude; weather data is provided by GFS.
        Uses heading-based algorithm to calculate UAV heading.
        Inputs: Desired Heading
        Prints and Returns: End Latitude, End Longtitude, Ideal End Latitude, Ideal End Longitude, End Heading, Desired Heading
        """
        # Fly for TimeStep1 seconds and recalculates dependent variables

        print 'Running Fly_2'
        # Processes desired heading to avoid degree/radian related issues
        if FinalHead < 0:
            FinalHead = 360 - (abs(FinalHead) % 360)
        elif FinalHead > 360:
            FinalHead = FinalHead % 360
        else:
            FinalHead = FinalHead

        if FinalHead > 315:
            # Makes self.FinalHead -ve in radians
            self.FinalHead = numpy.unwrap([0, numpy.deg2rad(FinalHead)])[1]
        else:
            self.FinalHead = numpy.deg2rad(FinalHead)

        logging.info('User Input:\nStart Latitude: %s\nStart Longtitude: %s\nStart Altitude: %s\nDesired Heading: %s\nStart Time: %s\nLift to Drag Ratio: %s\n \
                     Mass: %s\nWing Planform Area: %s\nCoefficient of Lift: %s', self.Lat, self.Lon, self.Alt, self.FinalHead, self.StartTime, self.LD,
                     self.Mass, self.PlanformArea, self.CL)

        while self._CurrentAlt > 100:
            # Fly for TimeStep1 seconds and recalculates dependent variables

            # Calculates heading of UAV
            self._CurrentHead = numpy.arctan2(self._YDirDump[-1], self._XDirDump[-1])

            # Avoids issues at 0, 180, 360 degress
            if self.FinalHead < 0:
                self._CurrentHead = self._CurrentHead
            elif self.FinalHead >= 0 and self.FinalHead < numpy.pi / 8:
                self._CurrentHead = self._CurrentHead
            else:
                if self._CurrentHead < 0:
                    self._CurrentHead += 2 * numpy.pi
                else:
                    self._CurrentHead = self._CurrentHead

            if abs(self._CurrentHead - self.FinalHead) > 0.8:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead
            elif abs(self._CurrentHead - self.FinalHead) > 0.5:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + (self.FinalHead - self._CurrentHead)
            elif abs(self._CurrentHead - self.FinalHead) > 0.3:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - 1.6 * (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + 1.6 * (self.FinalHead - self._CurrentHead)
            elif abs(self._CurrentHead - self.FinalHead) > 0.1:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - 2.7 * (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + 2.7 * (self.FinalHead - self._CurrentHead)
            elif abs(self._CurrentHead - self.FinalHead) > 0.05:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - 8 * (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + 8 * (self.FinalHead - self._CurrentHead)
            else:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - 16 * (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + 16 * (self.FinalHead - self._CurrentHead)
            self._HeadingDump.append(self._Heading)

            # Provides list of time values for ODE solver
            time = numpy.linspace(0, self.TimeStep1)

            # Variables required for ODE solver
            SinA = numpy.sin(self._GlideAngle)
            CosA = numpy.cos(self._GlideAngle)
            B = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass)
            C = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass * self.LD)

            # Sets up ODEs to be solved
            def z_derivatives(z, t):
                return numpy.array([-z[1] * SinA, - B * (z[1] ** 2) + 9.81 * SinA])

            def xy_derivatives(xy, t):
                return numpy.array([xy[1] * CosA, - C * (xy[1] ** 2) + 9.81 * CosA])

            # Solves for change in altitude
            z_initial_conditions = numpy.array([self._CurrentAlt, self._XZCurrentSpeed])
            z = odeint(z_derivatives, z_initial_conditions, time)
            self._CurrentAlt, self._XZCurrentSpeed = z[-1]

            # Solves for distance travelled
            xy_initial_conditions = numpy.array([self._CurrentDist, self._XZTempSpeed])
            xy = odeint(xy_derivatives, xy_initial_conditions, time)
            self._CurrentDist, self._XZTempSpeed = xy[-1]

            # Appends changes in X and Y directions to X and Y data dumps
            if self._Heading == 0:
                self._XDirDump.append(self._XDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            elif self._Heading == numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            elif self._Heading == numpy.pi:
                self._XDirDump.append(self._XDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            elif self._Heading == 3 * numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            else:
                self._XDirDump.append(self._XDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.cos(self._Heading)) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.sin(self._Heading)) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)

            # Appends the wind conditions to wind data dumps
            self._XWindDump.append(self._xyWindSpeed * numpy.cos(self._WindHead))
            self._YWindDump.append(self._xyWindSpeed * numpy.sin(self._WindHead))

            # Updates parameters used to obtain the next set of forecasts
            self._CurrentLatLon = m2deg(self._YDirDump[-1], self._XDirDump[-1], self.Lat)
            self._CurrentLat = self.Lat + self._CurrentLatLon[0]
            self._CurrentLon = self.Lon + self._CurrentLatLon[1]
            self._CurrentTime += timedelta(seconds=self.TimeStep1)

            # Appends useful data to data dumps
            self._TimeDump.append(self._TimeDump[-1] + self.TimeStep1)
            self._AltDump.append(self._CurrentAlt)
            self._DistDump.append(self._CurrentDist)

            logging.info('Current Latitude: %s\nCurrent Longitude: %s\nCurrent Altitude: %s\nCurrent Time: %s', self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)

            # Check if region travelled exceeds 3 degrees Lat and 3 degrees Lon from the launch site
            if abs(self._CurrentLat - self._BaseLat) < 1.5 and abs(self._CurrentLon - self._BaseLon) < 1.5:
                # Retrieves changing variables from forecast
                self._xyWindSpeed = self.WeatherData[-1].getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = self.WeatherData[-1].getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = self.WeatherData[-1].getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)
            else:
                print 'UAV has travelled out of forecast region.'
                print 'New forecast has to be obtained, please wait while forecast is downloaded and processed.'
                self._BaseLat = self._CurrentLat
                self._BaseLon = self._CurrentLon

                # Intialises Weather Environment object
                self.WeatherData.append(weather.forecastEnvironment())

                # Obtains data from GFS for current flight conditions
                self.WeatherData[-1].launchSiteLat = self._BaseLat
                self.WeatherData[-1].launchSiteLon = self._BaseLon
                self.WeatherData[-1].launchSiteElev = self._CurrentAlt
                self.WeatherData[-1].dateAndTime = self._CurrentTime
                self.WeatherData[-1].UTC_offset = getUTCOffset(self._BaseLat, self._BaseLon, self._CurrentTime)

                # Downloads GFS weather data
                self.WeatherData[-1].loadForecast()
                logging.info('GFS Data downloaded!')
                print "Forecast downloaded!"

                # Retrieves changing variables from forecast
                self._xyWindSpeed = self.WeatherData[-1].getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = self.WeatherData[-1].getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = self.WeatherData[-1].getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)

        while self._CurrentAlt > 0:
            # Fly for TimeStep2 seconds and recalculates dependent variables

            # Calculates heading of UAV
            self._CurrentHead = numpy.arctan2(self._YDirDump[-1], self._XDirDump[-1])

            # Avoids issues at 0, 180, 360 degress
            if self.FinalHead < 0:
                self._CurrentHead = self._CurrentHead
            elif self.FinalHead >= 0 and self.FinalHead < numpy.pi / 8:
                self._CurrentHead = self._CurrentHead
            else:
                if self._CurrentHead < 0:
                    self._CurrentHead += 2 * numpy.pi
                else:
                    self._CurrentHead = self._CurrentHead

            if abs(self._CurrentHead - self.FinalHead) > 0.8:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead
            elif abs(self._CurrentHead - self.FinalHead) > 0.5:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + (self.FinalHead - self._CurrentHead)
            elif abs(self._CurrentHead - self.FinalHead) > 0.3:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - 1.6 * (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + 1.6 * (self.FinalHead - self._CurrentHead)
            elif abs(self._CurrentHead - self.FinalHead) > 0.1:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - 2.7 * (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + 2.7 * (self.FinalHead - self._CurrentHead)
            elif abs(self._CurrentHead - self.FinalHead) > 0.05:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - 8 * (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + 8 * (self.FinalHead - self._CurrentHead)
            else:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - 16 * (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + 16 * (self.FinalHead - self._CurrentHead)
            self._HeadingDump.append(self._Heading)

            # Provides list of time values for ODE solver
            time = numpy.linspace(0, self.TimeStep2, num=100)

            # Variables required for ODE solver
            SinA = numpy.sin(self._GlideAngle)
            CosA = numpy.cos(self._GlideAngle)
            B = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass)
            C = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass * self.LD)

            # Sets up ODEs to be solved
            def z_derivatives(z, t):
                return numpy.array([-z[1] * SinA, - B * (z[1] ** 2) + 9.81 * SinA])

            def xy_derivatives(xy, t):
                return numpy.array([xy[1] * CosA, - C * (xy[1] ** 2) + 9.81 * CosA])

            # Solves for change in altitude
            z_initial_conditions = numpy.array([self._CurrentAlt, self._XZCurrentSpeed])
            z = odeint(z_derivatives, z_initial_conditions, time)
            self._CurrentAlt, self._XZCurrentSpeed = z[-1]

            # Solves for distance travelled
            xy_initial_conditions = numpy.array([self._CurrentDist, self._XZTempSpeed])
            xy = odeint(xy_derivatives, xy_initial_conditions, time)
            self._CurrentDist, self._XZTempSpeed = xy[-1]

            # Appends changes in X and Y directions to X and Y data dumps
            if self._Heading == 0:
                self._XDirDump.append(self._XDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)
            elif self._Heading == numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)
            elif self._Heading == numpy.pi:
                self._XDirDump.append(self._XDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)
            elif self._Heading == 3 * numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)
            else:
                self._XDirDump.append(self._XDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.cos(self._Heading)) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.sin(self._Heading)) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)

            # Appends the wind conditions to wind data dumps
            self._XWindDump.append(self._xyWindSpeed * numpy.cos(self._WindHead))
            self._YWindDump.append(self._xyWindSpeed * numpy.sin(self._WindHead))

            # Updates parameters used to obtain the next set of forecasts
            self._CurrentLatLon = m2deg(self._YDirDump[-1], self._XDirDump[-1], self.Lat)
            self._CurrentLat = self.Lat + self._CurrentLatLon[0]
            self._CurrentLon = self.Lon + self._CurrentLatLon[1]
            self._CurrentTime += timedelta(seconds=self.TimeStep2)

            # Appends useful data to data dumps
            self._TimeDump.append(self._TimeDump[-1] + self.TimeStep2)
            self._AltDump.append(self._CurrentAlt)
            self._DistDump.append(self._CurrentDist)

            logging.info('Current Latitude: %s\nCurrent Longitude: %s\nCurrent Altitude: %s\nCurrent Time: %s', self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)

            # Check if region travelled exceeds 3 degrees Lat and 3 degrees Lon from the launch site
            if abs(self._CurrentLat - self._BaseLat) < 1.5 and abs(self._CurrentLon - self._BaseLon) < 1.5:
                # Retrieves changing variables from forecast
                self._xyWindSpeed = self.WeatherData[-1].getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = self.WeatherData[-1].getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = self.WeatherData[-1].getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)

                # Checks if UAV is low enough for calculations to stop
                if self._CurrentAlt < 1:
                    break
            else:
                print 'UAV has travelled out of forecast region.'
                print 'New forecast has to be obtained, please wait while forecast is downloaded and processed.'
                self._BaseLat = self._CurrentLat
                self._BaseLon = self._CurrentLon

                # Intialises Weather Environment object
                self.WeatherData.append(weather.forecastEnvironment())

                # Obtains data from GFS for current flight conditions
                self.WeatherData[-1].launchSiteLat = self._BaseLat
                self.WeatherData[-1].launchSiteLon = self._BaseLon
                self.WeatherData[-1].launchSiteElev = self._CurrentAlt
                self.WeatherData[-1].dateAndTime = self._CurrentTime
                self.WeatherData[-1].UTC_offset = getUTCOffset(self._BaseLat, self._BaseLon, self._CurrentTime)

                # Downloads GFS weather data
                self.WeatherData[-1].loadForecast()
                logging.info('GFS Data downloaded!')
                print "Forecast downloaded!"

                # Retrieves changing variables from forecast
                self._xyWindSpeed = self.WeatherData[-1].getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = self.WeatherData[-1].getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = self.WeatherData[-1].getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)

                # Checks if UAV is low enough for calculations to stop
                if self._CurrentAlt < 1:
                    break

        # Prepares output
        self._FinalLatLon = m2deg(self._YDirDump[-1], self._XDirDump[-1], self.Lat)
        self._EndHead = numpy.rad2deg(numpy.arctan2(self._YDirDump[-1], self._XDirDump[-1]))
        if self._EndHead < 0 and FinalHead > 90:
            self._EndHead += 360

        logging.info('Final Latitude: %s\nFinal Longitude: %s\nDescent Time: %s', self._FinalLatLon[0] + self.Lat, self._FinalLatLon[1] + self.Lon, (self._CurrentTime - self.StartTime))
        print 'Final Latitude: {}\nFinal Longitude: {}\nDescent Time: {}'.format(self._FinalLatLon[0] + self.Lat, self._FinalLatLon[1] + self.Lon, (self._CurrentTime - self.StartTime))
        print 'Final Heading: {}\nDesired Heading: {}\nFinal Altitude: {}'.format(self._EndHead, FinalHead, self._CurrentAlt)

        # Sets variables for plotting
        WriteToFile_Simulator(self._StartTime, self._XDirDump, self._YDirDump, self._DensityDump, self._AltDump, self._TimeDump, self.TimeStep1)

        # Outputs
        return [self._XDirDump[-1], self._YDirDump[-1], self._CurrentDist * numpy.cos(self.FinalHead), self._CurrentDist * numpy.sin(self.FinalHead), FinalHead, self._EndHead]

    def Fly_3(self, FinalHead):
        """
        Simulates of gliding UAV flight from a provided altitude; weather data is provided by GFS.
        Uses heading-based algorithm to calculate UAV heading.
        Inputs: Desired Heading
        Prints and Returns: End Latitude, End Longtitude, Ideal End Latitude, Ideal End Longitude, End Heading, Desired Heading
        """
        # Fly for TimeStep1 seconds and recalculates dependent variables

        print 'Running Fly_3'
        # Processes desired heading to avoid degree/radian related issues
        if FinalHead < 0:
            FinalHead = 360 - (abs(FinalHead) % 360)
        elif FinalHead > 360:
            FinalHead = FinalHead % 360
        else:
            FinalHead = FinalHead

        if FinalHead > 315:
            # Makes self.FinalHead -ve in radians
            self.FinalHead = numpy.unwrap([0, numpy.deg2rad(FinalHead)])[1]
        else:
            self.FinalHead = numpy.deg2rad(FinalHead)

        logging.info('User Input:\nStart Latitude: %s\nStart Longtitude: %s\nStart Altitude: %s\nDesired Heading: %s\nStart Time: %s\nLift to Drag Ratio: %s\n \
                     Mass: %s\nWing Planform Area: %s\nCoefficient of Lift: %s', self.Lat, self.Lon, self.Alt, self.FinalHead, self.StartTime, self.LD,
                     self.Mass, self.PlanformArea, self.CL)

        while self._CurrentAlt > 0.6 * self.Alt:
            # Fly for TimeStep1 seconds and recalculates dependent variables

            # Calculates heading of UAV
            self._CurrentHead = numpy.arctan2(self._YDirDump[-1], self._XDirDump[-1])

            # Avoids issues at 0, 180, 360 degress
            if self.FinalHead < 0:
                self._CurrentHead = self._CurrentHead
            elif self.FinalHead >= 0 and self.FinalHead < numpy.pi / 8:
                self._CurrentHead = self._CurrentHead
            else:
                if self._CurrentHead < 0:
                    self._CurrentHead += 2 * numpy.pi
                else:
                    self._CurrentHead = self._CurrentHead

            if abs(self._CurrentHead - self.FinalHead) > 0.8:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead
            elif abs(self._CurrentHead - self.FinalHead) > 0.5:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + (self.FinalHead - self._CurrentHead)
            elif abs(self._CurrentHead - self.FinalHead) > 0.3:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - 1.6 * (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + 1.6 * (self.FinalHead - self._CurrentHead)
            elif abs(self._CurrentHead - self.FinalHead) > 0.1:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - 2.7 * (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + 2.7 * (self.FinalHead - self._CurrentHead)
            elif abs(self._CurrentHead - self.FinalHead) > 0.05:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - 8 * (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + 8 * (self.FinalHead - self._CurrentHead)
            else:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - 16 * (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + 16 * (self.FinalHead - self._CurrentHead)

            self._HeadingDump.append(self._Heading)

            # Provides list of time values for ODE solver
            time = numpy.linspace(0, self.TimeStep1)

            # Variables required for ODE solver
            SinA = numpy.sin(self._GlideAngle)
            CosA = numpy.cos(self._GlideAngle)
            B = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass)
            C = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass * self.LD)

            # Sets up ODEs to be solved
            def z_derivatives(z, t):
                return numpy.array([-z[1] * SinA, - B * (z[1] ** 2) + 9.81 * SinA])

            def xy_derivatives(xy, t):
                return numpy.array([xy[1] * CosA, - C * (xy[1] ** 2) + 9.81 * CosA])

            # Solves for change in altitude
            z_initial_conditions = numpy.array([self._CurrentAlt, self._XZCurrentSpeed])
            z = odeint(z_derivatives, z_initial_conditions, time)
            self._CurrentAlt, self._XZCurrentSpeed = z[-1]

            # Solves for distance travelled
            xy_initial_conditions = numpy.array([self._CurrentDist, self._XZTempSpeed])
            xy = odeint(xy_derivatives, xy_initial_conditions, time)
            self._CurrentDist, self._XZTempSpeed = xy[-1]

            # Appends changes in X and Y directions to X and Y data dumps
            if self._Heading == 0:
                self._XDirDump.append(self._XDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            elif self._Heading == numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            elif self._Heading == numpy.pi:
                self._XDirDump.append(self._XDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            elif self._Heading == 3 * numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            else:
                self._XDirDump.append(self._XDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.cos(self._Heading)) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.sin(self._Heading)) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)

            # Appends the wind conditions to wind data dumps
            self._XWindDump.append(self._xyWindSpeed * numpy.cos(self._WindHead))
            self._YWindDump.append(self._xyWindSpeed * numpy.sin(self._WindHead))

            # Updates parameters used to obtain the next set of forecasts
            self._CurrentLatLon = m2deg(self._YDirDump[-1], self._XDirDump[-1], self.Lat)
            self._CurrentLat = self.Lat + self._CurrentLatLon[0]
            self._CurrentLon = self.Lon + self._CurrentLatLon[1]
            self._CurrentTime += timedelta(seconds=self.TimeStep1)

            # Appends useful data to data dumps
            self._TimeDump.append(self._TimeDump[-1] + self.TimeStep1)
            self._AltDump.append(self._CurrentAlt)
            self._DistDump.append(self._CurrentDist)

            logging.info('Current Latitude: %s\nCurrent Longitude: %s\nCurrent Altitude: %s\nCurrent Time: %s', self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)

            # Check if region travelled exceeds 3 degrees Lat and 3 degrees Lon from the launch site
            if abs(self._CurrentLat - self._BaseLat) < 1.5 and abs(self._CurrentLon - self._BaseLon) < 1.5:
                # Retrieves changing variables from forecast
                self._xyWindSpeed = self.WeatherData[-1].getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = self.WeatherData[-1].getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = self.WeatherData[-1].getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)
            else:
                print 'UAV has travelled out of forecast region.'
                print 'New forecast has to be obtained, please wait while forecast is downloaded and processed.'
                self._BaseLat = self._CurrentLat
                self._BaseLon = self._CurrentLon

                # Intialises Weather Environment object
                self.WeatherData.append(weather.forecastEnvironment())

                # Obtains data from GFS for current flight conditions
                self.WeatherData[-1].launchSiteLat = self._BaseLat
                self.WeatherData[-1].launchSiteLon = self._BaseLon
                self.WeatherData[-1].launchSiteElev = self._CurrentAlt
                self.WeatherData[-1].dateAndTime = self._CurrentTime
                self.WeatherData[-1].UTC_offset = getUTCOffset(self._BaseLat, self._BaseLon, self._CurrentTime)

                # Downloads GFS weather data
                self.WeatherData[-1].loadForecast()
                logging.info('GFS Data downloaded!')
                print "Forecast downloaded!"

                # Retrieves changing variables from forecast
                self._xyWindSpeed = self.WeatherData[-1].getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = self.WeatherData[-1].getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = self.WeatherData[-1].getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)

        while self._CurrentAlt < 0.6 * self.Alt and self._CurrentAlt > 0.25 * self.Alt:
            # Fly for 3 seconds and recalculate dependent variables

            self._Heading = self.FinalHead
            self._HeadingDump.append(self._Heading)

            # Provides list of time values for ODE solver
            time = numpy.linspace(0, self.TimeStep1)

            # Variables required for ODE solver
            SinA = numpy.sin(self._GlideAngle)
            CosA = numpy.cos(self._GlideAngle)
            B = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass)
            C = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass * self.LD)

            # Sets up ODEs to be solved
            def z_derivatives(z, t):
                return numpy.array([-z[1] * SinA, - B * (z[1] ** 2) + 9.81 * SinA])

            def xy_derivatives(xy, t):
                return numpy.array([xy[1] * CosA, - C * (xy[1] ** 2) + 9.81 * CosA])

            # Solves for change in altitude
            z_initial_conditions = numpy.array([self._CurrentAlt, self._XZCurrentSpeed])
            z = odeint(z_derivatives, z_initial_conditions, time)
            self._CurrentAlt, self._XZCurrentSpeed = z[-1]

            # Solves for distance travelled
            xy_initial_conditions = numpy.array([self._CurrentDist, self._XZTempSpeed])
            xy = odeint(xy_derivatives, xy_initial_conditions, time)
            self._CurrentDist, self._XZTempSpeed = xy[-1]

            # Appends changes in X and Y directions to X and Y data dumps
            if self._Heading == 0:
                self._XDirDump.append(self._XDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            elif self._Heading == numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            elif self._Heading == numpy.pi:
                self._XDirDump.append(self._XDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            elif self._Heading == 3 * numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            else:
                self._XDirDump.append(self._XDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.cos(self._Heading)) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.sin(self._Heading)) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)

            # Appends the wind conditions to wind data dumps
            self._XWindDump.append(self._xyWindSpeed * numpy.cos(self._WindHead))
            self._YWindDump.append(self._xyWindSpeed * numpy.sin(self._WindHead))

            # Updates parameters used to obtain the next set of forecasts
            self._CurrentLatLon = m2deg(self._YDirDump[-1], self._XDirDump[-1], self.Lat)
            self._CurrentLat = self.Lat + self._CurrentLatLon[0]
            self._CurrentLon = self.Lon + self._CurrentLatLon[1]
            self._CurrentTime += timedelta(seconds=self.TimeStep1)

            # Appends useful data to data dumps
            self._TimeDump.append(self._TimeDump[-1] + self.TimeStep1)
            self._AltDump.append(self._CurrentAlt)
            self._DistDump.append(self._CurrentDist)

            logging.info('Current Latitude: %s\nCurrent Longitude: %s\nCurrent Altitude: %s\nCurrent Time: %s', self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)

            # Check if region travelled exceeds 3 degrees Lat and 3 degrees Lon from the launch site
            if abs(self._CurrentLat - self._BaseLat) < 1.5 and abs(self._CurrentLon - self._BaseLon) < 1.5:
                # Retrieves changing variables from forecast
                self._xyWindSpeed = self.WeatherData[-1].getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = self.WeatherData[-1].getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = self.WeatherData[-1].getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)
            else:
                print 'UAV has travelled out of forecast region.'
                print 'New forecast has to be obtained, please wait while forecast is downloaded and processed.'
                self._BaseLat = self._CurrentLat
                self._BaseLon = self._CurrentLon

                # Intialises Weather Environment object
                self.WeatherData.append(weather.forecastEnvironment())

                # Obtains data from GFS for current flight conditions
                self.WeatherData[-1].launchSiteLat = self._BaseLat
                self.WeatherData[-1].launchSiteLon = self._BaseLon
                self.WeatherData[-1].launchSiteElev = self._CurrentAlt
                self.WeatherData[-1].dateAndTime = self._CurrentTime
                self.WeatherData[-1].UTC_offset = getUTCOffset(self._BaseLat, self._BaseLon, self._CurrentTime)

                # Downloads GFS weather data
                self.WeatherData[-1].loadForecast()
                logging.info('GFS Data downloaded!')
                print "Forecast downloaded!"

                # Retrieves changing variables from forecast
                self._xyWindSpeed = self.WeatherData[-1].getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = self.WeatherData[-1].getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = self.WeatherData[-1].getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)

        while self._CurrentAlt < 0.25 * self.Alt:
            # Fly for TimeStep2 seconds and recalculates dependent variables

            # Calculates heading of UAV
            self._CurrentHead = numpy.arctan2(self._YDirDump[-1], self._XDirDump[-1])

            # Avoids issues at 0, 180, 360 degress
            if self.FinalHead < 0:
                self._CurrentHead = self._CurrentHead
            elif self.FinalHead >= 0 and self.FinalHead < numpy.pi / 8:
                self._CurrentHead = self._CurrentHead
            else:
                if self._CurrentHead < 0:
                    self._CurrentHead += 2 * numpy.pi
                else:
                    self._CurrentHead = self._CurrentHead

            if abs(self._CurrentHead - self.FinalHead) > 0.8:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead
            elif abs(self._CurrentHead - self.FinalHead) > 0.5:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + (self.FinalHead - self._CurrentHead)
            elif abs(self._CurrentHead - self.FinalHead) > 0.3:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - 1.6 * (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + 1.6 * (self.FinalHead - self._CurrentHead)
            elif abs(self._CurrentHead - self.FinalHead) > 0.1:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - 2.7 * (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + 2.7 * (self.FinalHead - self._CurrentHead)
            elif abs(self._CurrentHead - self.FinalHead) > 0.05:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - 8 * (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + 8 * (self.FinalHead - self._CurrentHead)
            else:
                if self._CurrentHead > self.FinalHead:
                    self._Heading = self.FinalHead - 16 * (self._CurrentHead - self.FinalHead)
                else:  # self.FinalHead > self._CurrentHead
                    self._Heading = self.FinalHead + 16 * (self.FinalHead - self._CurrentHead)

            self._HeadingDump.append(self._Heading)

            # Provides list of time values for ODE solver
            time = numpy.linspace(0, self.TimeStep2, num=100)

            # Variables required for ODE solver
            SinA = numpy.sin(self._GlideAngle)
            CosA = numpy.cos(self._GlideAngle)
            B = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass)
            C = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass * self.LD)

            # Sets up ODEs to be solved
            def z_derivatives(z, t):
                return numpy.array([-z[1] * SinA, - B * (z[1] ** 2) + 9.81 * SinA])

            def xy_derivatives(xy, t):
                return numpy.array([xy[1] * CosA, - C * (xy[1] ** 2) + 9.81 * CosA])

            # Solves for change in altitude
            z_initial_conditions = numpy.array([self._CurrentAlt, self._XZCurrentSpeed])
            z = odeint(z_derivatives, z_initial_conditions, time)
            self._CurrentAlt, self._XZCurrentSpeed = z[-1]

            # Solves for distance travelled
            xy_initial_conditions = numpy.array([self._CurrentDist, self._XZTempSpeed])
            xy = odeint(xy_derivatives, xy_initial_conditions, time)
            self._CurrentDist, self._XZTempSpeed = xy[-1]

            # Appends changes in X and Y directions to X and Y data dumps
            if self._Heading == 0:
                self._XDirDump.append(self._XDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)
            elif self._Heading == numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)
            elif self._Heading == numpy.pi:
                self._XDirDump.append(self._XDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)
            elif self._Heading == 3 * numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)
            else:
                self._XDirDump.append(self._XDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.cos(self._Heading)) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.sin(self._Heading)) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)

            # Appends the wind conditions to wind data dumps
            self._XWindDump.append(self._xyWindSpeed * numpy.cos(self._WindHead))
            self._YWindDump.append(self._xyWindSpeed * numpy.sin(self._WindHead))

            # Updates parameters used to obtain the next set of forecasts
            self._CurrentLatLon = m2deg(self._YDirDump[-1], self._XDirDump[-1], self.Lat)
            self._CurrentLat = self.Lat + self._CurrentLatLon[0]
            self._CurrentLon = self.Lon + self._CurrentLatLon[1]
            self._CurrentTime += timedelta(seconds=self.TimeStep2)

            # Appends useful data to data dumps
            self._TimeDump.append(self._TimeDump[-1] + self.TimeStep2)
            self._AltDump.append(self._CurrentAlt)
            self._DistDump.append(self._CurrentDist)

            logging.info('Current Latitude: %s\nCurrent Longitude: %s\nCurrent Altitude: %s\nCurrent Time: %s', self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)

            # Check if region travelled exceeds 3 degrees Lat and 3 degrees Lon from the launch site
            if abs(self._CurrentLat - self._BaseLat) < 1.5 and abs(self._CurrentLon - self._BaseLon) < 1.5:
                # Retrieves changing variables from forecast
                self._xyWindSpeed = self.WeatherData[-1].getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = self.WeatherData[-1].getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = self.WeatherData[-1].getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)

                # Checks if UAV is low enough for calculations to stop
                if self._CurrentAlt < 1:
                    break
            else:
                print 'UAV has travelled out of forecast region.'
                print 'New forecast has to be obtained, please wait while forecast is downloaded and processed.'
                self._BaseLat = self._CurrentLat
                self._BaseLon = self._CurrentLon

                # Intialises Weather Environment object
                self.WeatherData.append(weather.forecastEnvironment())

                # Obtains data from GFS for current flight conditions
                self.WeatherData[-1].launchSiteLat = self._BaseLat
                self.WeatherData[-1].launchSiteLon = self._BaseLon
                self.WeatherData[-1].launchSiteElev = self._CurrentAlt
                self.WeatherData[-1].dateAndTime = self._CurrentTime
                self.WeatherData[-1].UTC_offset = getUTCOffset(self._BaseLat, self._BaseLon, self._CurrentTime)

                # Downloads GFS weather data
                self.WeatherData[-1].loadForecast()
                logging.info('GFS Data downloaded!')
                print "Forecast downloaded!"

                # Retrieves changing variables from forecast
                self._xyWindSpeed = self.WeatherData[-1].getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = self.WeatherData[-1].getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = self.WeatherData[-1].getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)

                # Checks if UAV is low enough for calculations to stop
                if self._CurrentAlt < 1:
                    break

        # Prepares output
        self._FinalLatLon = m2deg(self._YDirDump[-1], self._XDirDump[-1], self.Lat)
        self._EndHead = numpy.rad2deg(numpy.arctan2(self._YDirDump[-1], self._XDirDump[-1]))
        if self._EndHead < 0 and FinalHead > 90:
            self._EndHead += 360

        logging.info('Final Latitude: %s\nFinal Longitude: %s\nDescent Time: %s', self._FinalLatLon[0] + self.Lat, self._FinalLatLon[1] + self.Lon, (self._CurrentTime - self.StartTime))
        print 'Final Latitude: {}\nFinal Longitude: {}\nDescent Time: {}'.format(self._FinalLatLon[0] + self.Lat, self._FinalLatLon[1] + self.Lon, (self._CurrentTime - self.StartTime))
        print 'Final Heading: {}\nDesired Heading: {}\nFinal Altitude: {}'.format(self._EndHead, FinalHead, self._CurrentAlt)

        # Sets variables for plotting

        # Outputs
        return [self._XDirDump[-1], self._YDirDump[-1], self._CurrentDist * numpy.cos(self.FinalHead), self._CurrentDist * numpy.sin(self.FinalHead), FinalHead, self._EndHead]

    def Fly_4(self, Final_X, Final_Y, Error, FinalHead):
        """
        Simulates of gliding UAV flight from a provided altitude; weather data is provided by GFS.
        Uses heading-based algorithm to calculate UAV heading.
        Inputs: Desired X Distance, Desired Y Distance, Error in X and Y Direction, Desired Heading
        """
        # Fly for TimeStep1 seconds and recalculates dependent variables
        self.Trigger = False
        self.Error = Error
        self._BaseXDiff, self._BaseYDiff = Final_X, Final_Y

        logging.info('User Input:\nStart Latitude: %s\nStart Longtitude: %s\End X: %s\End Y: %s\nStart Altitude: %s\nStart Time: %s\nLift to Drag Ratio: %s\n \
                     Mass: %s\nWing Planform Area: %s\nCoefficient of Lift: %s', self.Lat, self.Lon, self._BaseXDiff, self._BaseYDiff, self.Alt, self.StartTime, self.LD,
                     self.Mass, self.PlanformArea, self.CL)

        while self._CurrentAlt > 100:
            # Fly for TimeStep1 seconds and recalculates dependent variables

            # Calculates heading of UAV
            self._XDiff, self._YDiff = self._XDirDump[-1], self._YDirDump[-1]

            if (self._BaseXDiff - self._XDiff) == 0 and (self._BaseYDiff - self._YDiff) > 0:
                self._Heading = numpy.pi / 2
            elif (self._BaseXDiff - self._XDiff) < 0 and (self._BaseYDiff - self._YDiff) == 0:
                self._Heading = numpy.pi
            elif (self._BaseXDiff - self._XDiff) == 0 and (self._BaseYDiff - self._YDiff) < 0:
                self._Heading = 3 * numpy.pi / 2
            else:
                self._Heading = numpy.arctan2(self._BaseYDiff - self._YDiff, self._BaseXDiff - self._XDiff)
            self._HeadingDump.append(self._Heading)

            # Provides list of time values for ODE solver
            time = numpy.linspace(0, self.TimeStep1)

            # Variables required for ODE solver
            SinA = numpy.sin(self._GlideAngle)
            CosA = numpy.cos(self._GlideAngle)
            B = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass)
            C = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass * self.LD)

            # Sets up ODEs to be solved
            def z_derivatives(z, t):
                return numpy.array([-z[1] * SinA, - B * (z[1] ** 2) + 9.81 * SinA])

            def xy_derivatives(xy, t):
                return numpy.array([xy[1] * CosA, - C * (xy[1] ** 2) + 9.81 * CosA])

            # Solves for change in altitude
            z_initial_conditions = numpy.array([self._CurrentAlt, self._XZCurrentSpeed])
            z = odeint(z_derivatives, z_initial_conditions, time)
            self._CurrentAlt, self._XZCurrentSpeed = z[-1]

            # Solves for distance travelled
            xy_initial_conditions = numpy.array([self._CurrentDist, self._XZTempSpeed])
            xy = odeint(xy_derivatives, xy_initial_conditions, time)
            self._CurrentDist, self._XZTempSpeed = xy[-1]

            # Appends changes in X and Y directions to X and Y data dumps
            if self._Heading == 0:
                self._XDirDump.append(self._XDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            elif self._Heading == numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            elif self._Heading == numpy.pi:
                self._XDirDump.append(self._XDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            elif self._Heading == 3 * numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)
            else:
                self._XDirDump.append(self._XDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.cos(self._Heading)) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep1)
                self._YDirDump.append(self._YDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.sin(self._Heading)) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep1)

            # Appends the wind conditions to wind data dumps
            self._XWindDump.append(self._xyWindSpeed * numpy.cos(self._WindHead))
            self._YWindDump.append(self._xyWindSpeed * numpy.sin(self._WindHead))

            # Updates parameters used to obtain the next set of forecasts
            self._CurrentLatLon = m2deg(self._YDirDump[-1], self._XDirDump[-1], self.Lat)
            self._CurrentLat = self.Lat + self._CurrentLatLon[0]
            self._CurrentLon = self.Lon + self._CurrentLatLon[1]
            self._CurrentTime += timedelta(seconds=self.TimeStep1)

            # Appends useful data to data dumps
            self._TimeDump.append(self._TimeDump[-1] + self.TimeStep1)
            self._AltDump.append(self._CurrentAlt)
            self._DistDump.append(self._CurrentDist)

            if abs(self._BaseXDiff - self._XDirDump[-1]) < self.Error and abs(self._BaseYDiff - self._YDirDump[-1]) < self.Error:
                self._EndHead = numpy.rad2deg(numpy.arctan2(self._YDirDump[-1], self._XDirDump[-1]))
                self.x_coords.append(self._XDirDump[-1])
                self.y_coords.append(self._YDirDump[-1])
                self.x_coords_end.append(Final_X)
                self.y_coords_end.append(Final_Y)
                self.set_heading.append(FinalHead)
                self.end_heading.append(self._EndHead)
                return

            logging.info('Current Latitude: %s\nCurrent Longitude: %s\nCurrent Altitude: %s\nCurrent Time: %s', self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)

            # Check if region travelled exceeds 3 degrees Lat and 3 degrees Lon from the launch site
            if self.Trigger is False:
                if abs(self._CurrentLat - self._BaseLat) < 1.5 and abs(self._CurrentLon - self._BaseLon) < 1.5:
                    # Retrieves changing variables from forecast
                    self._xyWindSpeed = self.WeatherData[-1].getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                    self._WindHead = self.WeatherData[-1].getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                    self._CurrentDensity = self.WeatherData[-1].getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                    self._DensityDump.append(self._CurrentDensity)
                else:
                    print 'UAV has travelled out of forecast region.'
                    print 'Assume UAV stays at border of forecast region, hold xyWindSpeed, WindHead, CurrentDensity'

                    self._xyWindSpeed = self._xyWindSpeed
                    self._WindHead = self._WindHead
                    self._CurrentDensity = self._CurrentDensity
                    self._DensityDump.append(self._CurrentDensity)
            else:
                pass

        while self._CurrentAlt > 0:
            # Fly for TimeStep2 seconds and recalculates dependent variables
            self._XDiff, self._YDiff = self._XDirDump[-1], self._YDirDump[-1]

            if (self._BaseXDiff - self._XDiff) == 0 and (self._BaseYDiff - self._YDiff) > 0:
                self._Heading = numpy.pi / 2
            elif (self._BaseXDiff - self._XDiff) < 0 and (self._BaseYDiff - self._YDiff) == 0:
                self._Heading = numpy.pi
            elif (self._BaseXDiff - self._XDiff) == 0 and (self._BaseYDiff - self._YDiff) < 0:
                self._Heading = 3 * numpy.pi / 2
            else:
                self._Heading = numpy.arctan2(self._BaseYDiff - self._YDiff, self._BaseXDiff - self._XDiff)
            self._HeadingDump.append(self._Heading)

            # Provides list of time values for ODE solver
            time = numpy.linspace(0, self.TimeStep2, num=100)

            # Variables required for ODE solver
            SinA = numpy.sin(self._GlideAngle)
            CosA = numpy.cos(self._GlideAngle)
            B = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass)
            C = 0.5 * self._CurrentDensity * self.PlanformArea * self.CL / (self.Mass * self.LD)

            # Sets up ODEs to be solved
            def z_derivatives(z, t):
                return numpy.array([-z[1] * SinA, - B * (z[1] ** 2) + 9.81 * SinA])

            def xy_derivatives(xy, t):
                return numpy.array([xy[1] * CosA, - C * (xy[1] ** 2) + 9.81 * CosA])

            # Solves for change in altitude
            z_initial_conditions = numpy.array([self._CurrentAlt, self._XZCurrentSpeed])
            z = odeint(z_derivatives, z_initial_conditions, time)
            self._CurrentAlt, self._XZCurrentSpeed = z[-1]

            # Solves for distance travelled
            xy_initial_conditions = numpy.array([self._CurrentDist, self._XZTempSpeed])
            xy = odeint(xy_derivatives, xy_initial_conditions, time)
            self._CurrentDist, self._XZTempSpeed = xy[-1]

            # Appends changes in X and Y directions to X and Y data dumps
            if self._Heading == 0:
                self._XDirDump.append(self._XDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)
            elif self._Heading == numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] + (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)
            elif self._Heading == numpy.pi:
                self._XDirDump.append(self._XDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)
            elif self._Heading == 3 * numpy.pi / 2:
                self._XDirDump.append(self._XDirDump[-1] + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] - (self._CurrentDist - self._DistDump[-1]) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)
            else:
                self._XDirDump.append(self._XDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.cos(self._Heading)) + self._xyWindSpeed * numpy.cos(self._WindHead) * self.TimeStep2)
                self._YDirDump.append(self._YDirDump[-1] + ((self._CurrentDist - self._DistDump[-1]) * numpy.sin(self._Heading)) + self._xyWindSpeed * numpy.sin(self._WindHead) * self.TimeStep2)

            # Appends the wind conditions to wind data dumps
            self._XWindDump.append(self._xyWindSpeed * numpy.cos(self._WindHead))
            self._YWindDump.append(self._xyWindSpeed * numpy.sin(self._WindHead))

            # Updates parameters used to obtain the next set of forecasts
            self._CurrentLatLon = m2deg(self._YDirDump[-1], self._XDirDump[-1], self.Lat)
            self._CurrentLat = self.Lat + self._CurrentLatLon[0]
            self._CurrentLon = self.Lon + self._CurrentLatLon[1]
            self._CurrentTime += timedelta(seconds=self.TimeStep2)

            # Appends useful data to data dumps
            self._TimeDump.append(self._TimeDump[-1] + self.TimeStep2)
            self._AltDump.append(self._CurrentAlt)
            self._DistDump.append(self._CurrentDist)

            if abs(self._BaseXDiff - self._XDirDump[-1]) < self.Error and abs(self._BaseYDiff - self._YDirDump[-1]) < self.Error:
                self._EndHead = numpy.rad2deg(numpy.arctan2(self._YDirDump[-1], self._XDirDump[-1]))
                self.x_coords.append(self._XDirDump[-1])
                self.y_coords.append(self._YDirDump[-1])
                self.x_coords_end.append(Final_X)
                self.y_coords_end.append(Final_Y)
                self.set_heading.append(FinalHead)
                self.end_heading.append(self._EndHead)
                return

            logging.info('Current Latitude: %s\nCurrent Longitude: %s\nCurrent Altitude: %s\nCurrent Time: %s', self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)

            # Check if region travelled exceeds 3 degrees Lat and 3 degrees Lon from the launch site
            if self.Trigger is False:
                if abs(self._CurrentLat - self._BaseLat) < 1.5 and abs(self._CurrentLon - self._BaseLon) < 1.5:
                    # Retrieves changing variables from forecast
                    self._xyWindSpeed = self.WeatherData[-1].getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                    self._WindHead = self.WeatherData[-1].getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                    self._CurrentDensity = self.WeatherData[-1].getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                    self._DensityDump.append(self._CurrentDensity)
                    # Checks if UAV is low enough for calculations to stop
                    if self._CurrentAlt < 1:
                        self._EndHead = numpy.rad2deg(numpy.arctan2(self._YDirDump[-1], self._XDirDump[-1]))
                        self.x_coords_fail.append(self._XDirDump[-1])
                        self.y_coords_fail.append(self._YDirDump[-1])
                        self.x_coords_end_fail.append(Final_X)
                        self.y_coords_end_fail.append(Final_Y)
                        self.set_heading_fail.append(FinalHead)
                        self.end_heading_fail.append(self._EndHead)
                        return
                else:
                    print 'UAV has travelled out of forecast region.'
                    print 'Assume UAV stays at border of forecast region, hold xyWindSpeed, WindHead, CurrentDensity'

                    self._xyWindSpeed = self._xyWindSpeed
                    self._WindHead = self._WindHead
                    self._CurrentDensity = self._CurrentDensity
                    self._DensityDump.append(self._CurrentDensity)

                    # Checks if UAV is low enough for calculations to stop
                    if self._CurrentAlt < 1:
                        self._EndHead = numpy.rad2deg(numpy.arctan2(self._YDirDump[-1], self._XDirDump[-1]))
                        self.x_coords_fail.append(self._XDirDump[-1])
                        self.y_coords_fail.append(self._YDirDump[-1])
                        self.x_coords_end_fail.append(Final_X)
                        self.y_coords_end_fail.append(Final_Y)
                        self.set_heading_fail.append(FinalHead)
                        self.end_heading_fail.append(self._EndHead)
                        return
            else:
                if self._CurrentAlt < 1:
                        self._EndHead = numpy.rad2deg(numpy.arctan2(self._YDirDump[-1], self._XDirDump[-1]))
                        self.x_coords_fail.append(self._XDirDump[-1])
                        self.y_coords_fail.append(self._YDirDump[-1])
                        self.x_coords_end_fail.append(Final_X)
                        self.y_coords_end_fail.append(Final_Y)
                        self.set_heading_fail.append(FinalHead)
                        self.end_heading_fail.append(self._EndHead)
                        return
                else:
                    pass

    def PlotContour(self, Fig1, Fig2, Fig3, Fig4, Fig5):
        '''
        Plots relevant figures to visualise the simulation process.
        Figure 1: Density of Air Against Time
        Figure 2: Altitude Against Time
        Figure 3: Density of Ait Against Altitude
        Figure 4: Predicted Flight Path of UAV in 2D
        Figure 5: Predicted Flight Path of UAV in 3D
        '''
        if Fig1 is True:
            plt.figure(1)
            plt.plot(self._TimeDump, self._DensityDump, 'r--', label='Density')
            plt.title('Density of Air Against Time')
            plt.xlabel('Time/s')
            plt.ylabel('Density of Air/kg/m^3')
            plt.legend()

        if Fig2 is True:
            plt.figure(2)
            plt.plot(self._TimeDump, self._AltDump, label='Altitude')
            plt.title('Altitude Against Time')
            plt.xlabel('Time/s')
            plt.ylabel('Altitude/m')
            plt.legend()

        if Fig3 is True:
            plt.figure(3)
            plt.plot(self._AltDump, self._DensityDump, 'g', label='Density')
            plt.title('Density of Air Against Altitude')
            plt.xlabel('Altitude/m')
            plt.ylabel('Density of Air/kg/m^3')
            plt.legend()

        if Fig4 is True:
            plt.figure(4)
            # plt.plot(self._CurrentDist * numpy.cos(self.FinalHead), self._CurrentDist * numpy.sin(self.FinalHead), 'bo', label='Desired End Point')
            plt.plot(0, 0, 'ro', label='Start Point')
            plt.plot(self._XDirDump, self._YDirDump, 'b-', label='Flight Path')
            plt.plot(self._XDirDump[-1], self._YDirDump[-1], 'go', label='End Point')
            plt.title('Predicted Flight Path of UAV')
            plt.xlabel('Distance in X Direction/m')
            plt.ylabel('Distance in Y Direction/m')
            plt.legend()

        if Fig5 is True:
            fig5 = plt.figure(5)
            axes5 = fig5.gca(projection='3d')
            axes5.plot([self._CurrentDist * numpy.cos(self.FinalHead)], [self._CurrentDist * numpy.sin(self.FinalHead)], [0], 'bo', label='Desired End Point')
            axes5.plot([0], [0], [self.Alt], 'ro', label='Start Point')
            axes5.plot(self._XDirDump, self._YDirDump, self._AltDump, label='Flight Path')
            axes5.plot([self._XDirDump[-1]], [self._YDirDump[-1]], [self._AltDump[-1]], 'go', label='End Point')
            axes5.set_title('Predicted Flight Path of UAV in 3D')
            axes5.set_xlabel('Distance in X Direction/m')
            axes5.set_ylabel('Distance in Y Direction/m')
            axes5.set_zlabel('Distance in Z Direction/m')
            axes5.legend()

        plt.show()

    def Fly_Range(self, Fly, FinalHead, Fig1, Fig2, Fig3, Fig4, Fig5):
        """
        Allows user to run Fly_1/Fly_2/Fly_3 multiple times while capturing the data and writing it to a Python file for further prcoessing.
        Inputs: Fly_1/Fly_2/Fly_3, [Desired Heading/s], Fig1, Fig2, Fig3, Fig4, Fig5
        Plots the results of Fly_1/Fly_2/Fly_3 and writes to a Python file named plot_contours_date.py
        """
        # Creates list for UAV directions
        uav_direction_range = FinalHead

        # Creates empty lists to store x and y coordinates and defined heading for UAV
        x_coords = []
        y_coords = []
        x_coords_end = []
        y_coords_end = []
        set_heading = []
        end_heading = []

        # Runs for-loop to collect information for contour plot
        for x in uav_direction_range:
            print 'Calculations started for the UAV headed towards {} degrees'.format(x)
            if Fly == 1:
                y = self.Fly_1(x)
                self.PlotContour(Fig1, Fig2, Fig3, Fig4, Fig5)
                x_coords.append(y[0])
                y_coords.append(y[1])
                x_coords_end.append(y[2])
                y_coords_end.append(y[3])
                set_heading.append(y[4])
                end_heading.append(y[5])
                WriteToFile_Dist(self._StartTime, y[0], y[1], y[4], y[5])
            if Fly == 2:
                y = self.Fly_2(x)
                self.PlotContour(Fig1, Fig2, Fig3, Fig4, Fig5)
                x_coords.append(y[0])
                y_coords.append(y[1])
                x_coords_end.append(y[2])
                y_coords_end.append(y[3])
                set_heading.append(y[4])
                end_heading.append(y[5])
                WriteToFile_Head1(self._StartTime, y[0], y[1], y[4], y[5])
            if Fly == 3:
                y = self.Fly_3(x)
                self.PlotContour(Fig1, Fig2, Fig3, Fig4, Fig5)
                x_coords.append(y[0])
                y_coords.append(y[1])
                x_coords_end.append(y[2])
                y_coords_end.append(y[3])
                set_heading.append(y[4])
                end_heading.append(y[5])
                WriteToFile_Head2(self._StartTime, y[0], y[1], y[4], y[5])
            self.Clear()
            print 'Calculations completed for the UAV headed towards {} degrees'.format(x)

        # Plots collected data
        if len(FinalHead) > 1:
            plt.plot(x_coords, y_coords, 'bo', label='Displacement Contour')
            plt.plot(x_coords_end, y_coords_end, 'go', label='Defined Contour')
            plt.plot(0, 0, 'ro', label="Start")

            plt.title('UAV Displacement in Y against UAV Displacement in X\n')
            plt.xlabel('UAV Displacement in X')
            plt.ylabel('UAV Displacement in Y')

            plt.legend()

            plt.show()

    def Fly_Range_2(self, StartAngle, EndAngle, StartDist, EndDist, DistDiff, Error):
        """
        Allows user to run Fly_4 multiple times while capturing the data and writing it to a Python file for further prcoessing.
        Inputs: Start Angle, End Angle, First Set of Distances to be Reached, Last Set of Distances to be Reached, Steps between Distances, Error allowed for both X and Y Directions
        Plots the results of Fly_4 and writes to a Python file named plot_contours_date_error.py
        """

        # Creates list for UAV directions
        uav_direction_range = numpy.arange(StartAngle, EndAngle + 1, 1)

        # Runs for-loop to collect information for contour plot
        for x in uav_direction_range:
            print 'Calculations started for the UAV headed towards {} degrees'.format(x)
            x_1 = numpy.deg2rad(x)
            # Generates waypoints for the UAV to fly to
            self.Generate_XY(x_1, StartDist, EndDist + 1, DistDiff)
            for z in self.xy:
                self.Fly_4(z[0], z[1], Error, x)
                # Saves the furthest successful flight
                if len(self.x_coords) > 1:
                    if abs(self.x_coords[-1] ** 2 + self.y_coords[-1] ** 2) > abs(self.x_coords[-2] ** 2 + self.y_coords[-2] ** 2):
                        self.x_coords_angle_max, self.y_coords_angle_max = self.x_coords[-1], self.y_coords[-1]
                self.Clear()
            # Appends the furthest successful flight of the set angle and clear variables for next angle
            try:
                self.x_coords_max.append(self.x_coords_angle_max)
                self.y_coords_max.append(self.y_coords_angle_max)
                self.x_coords_angle_max, self.y_coords_angle_max = None, None
            except:
                pass
            print 'Calculations completed for the UAV headed towards {} degrees'.format(x)

        WriteToFile_Local(self._StartTime, self.x_coords, self.y_coords, self.x_coords_max, self.y_coords_max, self.x_coords_end, self.y_coords_end, self.set_heading, self.end_heading, self.x_coords_fail, self.y_coords_fail, self.x_coords_end_fail, self.y_coords_end_fail, self.set_heading_fail, self.end_heading_fail, self.Error)

        # Plots collected data
        plt.plot(self.x_coords, self.y_coords, 'bo', label='Displacements')
        plt.plot(self.x_coords_max, self.y_coords_max, 'yo', linestyle='-', label='Flight Contour')
        if len(self.x_coords) > 1:
            plt.fill(self.x_coords_max, self.y_coords_max, 'y', alpha=0.5)
            plt.plot(self.x_coords_end, self.y_coords_end, 'go', label='Defined Contour')

        plt.plot(0, 0, 'ro', label="Start")

        plt.title('Flight Contour\n')
        plt.xlabel('UAV Displacement in X')
        plt.ylabel('UAV Displacement in Y')

        plt.legend()

        plt.show()

    def Generate_XY(self, FinalHead, StartDist, EndDist, DistDiff):
        '''
        Generates the distances for the UAV to fly to for Fly_Range_2
        '''
        x = numpy.arange(StartDist, EndDist, DistDiff)
        x = [z * numpy.cos(FinalHead) for z in x]
        y = numpy.arange(StartDist, EndDist, DistDiff)
        y = [z * numpy.sin(FinalHead) for z in y]
        self.xy = zip(x, y)
