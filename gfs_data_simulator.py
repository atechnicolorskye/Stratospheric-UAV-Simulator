"""
gfs_data_simulator.py
GFS Data Simulator

DESCRIPTION
-----------

GFS Data Simulator Module:
This module simulates the descent of a UAV through multiple atmopsheric layers. The properties of the layers are obtained from National
Oceanic and Atmospheric Administration's  (NOAA) Global Forecast System (GFS) using Niccolo' Zapponi's GFS and weather modules. Using
code from the base_simulator module, displacement, heading and flight time are found layer-by-layer. The location and other relevant data
of the UAV are updated every x seconds of flight until it is above 100m where updates are carried out every y seconds. The UAV is assumed
to have been carried up by a weather balloon and there are no updrafts and downdrafts.

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

The code requires the user to create an instance of the Flight_LatLon() object. The user then calls the method Fly that predict the
trajectory of the UAV on its descent. Fly takes variables in the order of starting latitude, starting longitude, desired latitude, desired
longitude, altitude, time, average lift-to-drag ratio, mass, planform area, average coefficient of lift, time step at altitudes above 100m
and time step at altitudes below 100m. Once Fly has been called, the user can call the method PlotContour that plots relevant figures to
visualise the descent process.

Example:
y = Flight_LatLon()
y.Fly_LatLon(50.1, -5.00, 52, -3, 30000, (2014,03,01,19,30,01), 10, 2, 0.6, 0.4, 3, 0.5)
y.PlotContour(True, True, False, True, True)


University of Southampton
Si Kai Lee, skl2g11@soton.ac.uk
"""

__author__ = "Si Kai Lee, University of Southampton, skl2g11@soton.ac.uk"


from global_tools import m2deg, deg2m, getUTCOffset
import weather

from datetime import datetime, timedelta
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from scipy.integrate import odeint
import time as t  # t.time()
from write_file import WriteToFile_Simulator

import logging
from sys import exit


class Flight_LatLon(object):
    """
    Class that contains the simulator of the simulation of UAV flight from a provided altitude.
    Instantiates the Flight_LatLon object when called.
    """
    def __init__(self):
        # Initialises flight container object

        # User defined variables
        self.Lat = 0
        self.Lon = 0
        self.ToLat = 0
        self.ToLon = 0
        self.Alt = 0
        self.StartTime = 0
        self.LD = 0
        self.Mass = 0
        self.PlanformArea = 0
        self.CL = 0
        self.TimeStep1 = 0
        self.TimeStep2 = 0

        # Dependent variables
        self._BaseLat = 0
        self._BaseLon = 0
        self._xyWindSpeed = 0
        self._WindHead = 0
        self._GlideAngle = 0
        self._XZCurrentSpeed = 0
        self._XZTempSpeed = 0
        self._CurrentDensity = 0
        self._CurrentLatLon = 0
        self._CurrentLat = 0
        self._CurrentLon = 0
        self._CurrentAlt = 0
        self._CurrentDist = 0
        self._CurrentTime = 0
        self._TimeDump = [0]
        self._AltDump = []
        self._DistDump = [0]
        self._DensityDump = []
        self._XDirDump = [0]
        self._YDirDump = [0]
        self._HeadingDump = []
        self._XZDump = [0]

        # Output
        self.output = []

    def Fly_LatLon(self, Lat, Lon, ToLat, ToLon, Alt, StartTime, LD, Mass, PlanformArea, CL, TimeStep1, TimeStep2):
        """
        Simulates of gliding UAV flight from a provided altitude; weather data is provided by GFS.
        Uses distance-based algorithm to calculate UAV heading.
        Inputs: Starting Latitutde, Starting Longitude, Starting Altitude, Desired Heading, Starting Time,
                Average Lift-to-Drag Ratio, Mass of UAV, PlanformArea (Assuming UAV is a flying wind),
                Average Coefficient of Lift, Time Step used for ODEs at altitudes above 100m,
                Time Step used for ODEs at altitudes below 100m.
        Prints End Latitude, End Longtitude, Descent Time
        """
        # Sets up logging configuration
        logging.basicConfig(filename='Fly.log', format='%(levelname)s:%(message)s', level=logging.DEBUG)

        # Intialises Weather Environment object
        WeatherData = weather.forecastEnvironment()

        # Populates user defined variables with user inputs
        self.Lat = Lat
        self.Lon = Lon
        self.ToLat = ToLat
        self.ToLon = ToLon
        self.Alt = Alt
        self._StartTime = StartTime
        self.StartTime = datetime(*StartTime)
        self.LD = LD
        self.Mass = Mass
        self.PlanformArea = PlanformArea
        self.CL = CL
        self.TimeStep1 = TimeStep1
        self.TimeStep2 = TimeStep2

        logging.info('User Input:\nStart Latitude: %s\nStart Longtitude: %s\End Latitude: %s\End Longtitude: %s\nStart Altitude: %s\nStart Time: %s\nLift to Drag Ratio: %s\n \
                     Mass: %s\nWing Planform Area: %s\nCoefficient of Lift: %s', self.Lat, self.Lon, self.ToLat, self.ToLon, self.Alt, self.StartTime, self.LD,
                     self.Mass, self.PlanformArea, self.CL)

        # Checks if StartTime is between 30 days before and 6 after current time
        if self.StartTime > datetime.now():
            if self.StartTime - datetime.now() > timedelta(days=6):
                raise Exception('Your selected date and time is greater than 6 days in the future. Please try again.')
                exit()
        else:
            if datetime.now() - self.StartTime > timedelta(days=30):
                raise Exception('Your selected date and time is greater than 30 days in the past. Please try again.')
                exit()

        # Obtains data from GFS for current flight conditions
        WeatherData.launchSiteLat = self.Lat
        WeatherData.launchSiteLon = self.Lon
        WeatherData.launchSiteElev = self.Alt
        WeatherData.dateAndTime = self.StartTime
        WeatherData.UTC_offset = getUTCOffset(self.Lat, self.Lon, self.StartTime)

        # Downloads GFS weather data
        print "Downloading the forecast (might take a while)..."
        WeatherData.loadForecast()
        logging.info('GFS Data downloaded!')
        print "Forecast downloaded!"

        # Get relevant properties of current air mass
        self._xyWindSpeed = WeatherData.getWindSpeed(self.Lat, self.Lon, self.Alt, self.StartTime)
        logging.info('XY Wind Speed: %s', self._xyWindSpeed)

        self._WindHead = WeatherData.getWindDirection(self.Lat, self.Lon, self.Alt, self.StartTime)
        logging.info('Wind Heading: %s', self._WindHead)

        self._CurrentDensity = WeatherData.getDensity(self.Lat, self.Lon, self.Alt, self.StartTime)
        logging.info('Current Density: %s', self._CurrentDensity)

        # Appends Alt and
        self._AltDump.append(self.Alt)
        self._DensityDump.append(self._CurrentDensity)

        # Sets up dependent variables
        self._BaseLat = self.Lat
        self._BaseLon = self.Lon
        self._CurrentLat = self.Lat
        self._CurrentLon = self.Lon
        self._CurrentAlt = self.Alt
        self._CurrentTime = self.StartTime
        self._GlideAngle = numpy.arctan(float(1) / self.LD)

        # Translates the desired landing location for the UAV in x and y coordinates
        self._BaseYDiff, self._BaseXDiff = deg2m(self.ToLat - self._BaseLat, self.ToLon - self._BaseLon, self._BaseLat)
        logging.info('self._BaseXDiff: %s  self._BaseYDiff: %s', self._BaseXDiff, self._BaseYDiff)

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

            self._XZDump.append(self._XZTempSpeed)

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
                self._xyWindSpeed = WeatherData.getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = WeatherData.getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = WeatherData.getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)
            else:
                print 'UAV has travelled out of forecast region.'
                print 'New forecast has to be obtained, please wait while forecast is downloaded and processed.'
                self._BaseLat = self._CurrentLat
                self._BaseLon = self._CurrentLon

                # Intialises Weather Environment object
                WeatherData = weather.forecastEnvironment()

                # Obtains data from GFS for current flight conditions
                WeatherData.launchSiteLat = self._BaseLat
                WeatherData.launchSiteLon = self._BaseLon
                WeatherData.launchSiteElev = self._CurrentAlt
                WeatherData.dateAndTime = self._CurrentTime
                WeatherData.UTC_offset = getUTCOffset(self._BaseLat, self._BaseLon, self._CurrentTime)

                # Downloads GFS weather data
                WeatherData.loadForecast()
                logging.info('GFS Data downloaded!')
                print "Forecast downloaded!"

                # Retrieves changing variables from forecast
                self._xyWindSpeed = WeatherData.getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = WeatherData.getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = WeatherData.getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)

        while self._CurrentAlt > 0:
            # Fly for TimeStep2 seconds and recalculates dependent variables

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

            self._XZDump.append(self._XZTempSpeed)
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
                self._xyWindSpeed = WeatherData.getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = WeatherData.getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = WeatherData.getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)
                if self._CurrentAlt < 1:
                    break
            else:
                print 'UAV has travelled out of forecast region.'
                print 'New forecast has to be obtained, please wait while forecast is downloaded and processed.'
                self._BaseLat = self._CurrentLat
                self._BaseLon = self._CurrentLon

                # Intialises Weather Environment object
                WeatherData = weather.forecastEnvironment()

                # Obtains data from GFS for current flight conditions
                WeatherData.launchSiteLat = self._BaseLat
                WeatherData.launchSiteLon = self._BaseLon
                WeatherData.launchSiteElev = self._CurrentAlt
                WeatherData.dateAndTime = self._CurrentTime
                WeatherData.UTC_offset = getUTCOffset(self._BaseLat, self._BaseLon, self._CurrentTime)

                # Downloads GFS weather data
                WeatherData.loadForecast()
                logging.info('GFS Data downloaded!')
                print "Forecast downloaded!"

                # Retrieves changing variables from forecast
                self._xyWindSpeed = WeatherData.getWindSpeed(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._WindHead = WeatherData.getWindDirection(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._CurrentDensity = WeatherData.getDensity(self._CurrentLat, self._CurrentLon, self._CurrentAlt, self._CurrentTime)
                self._DensityDump.append(self._CurrentDensity)
                if self._CurrentAlt < 1:
                    break

        # Outputs
        self._FinalLatLon = m2deg(self._YDirDump[-1], self._XDirDump[-1], self.Lat)
        logging.info('Final Latitude: %s\nFinal Longitude: %s\nDescent Time: %s', self._FinalLatLon[0] + self.Lat, self._FinalLatLon[1] + self.Lon, (self._CurrentTime - self.StartTime))
        print 'Final Latitude: {}\nFinal Longitude: {}\nDescent Time: {}'.format(self._FinalLatLon[0] + self.Lat, self._FinalLatLon[1] + self.Lon, (self._CurrentTime - self.StartTime))
        WriteToFile_Simulator(self._StartTime, self._XDirDump, self._YDirDump, self._DensityDump, self._AltDump, self._TimeDump, self.TimeStep1)

    def PlotContour_LatLon(self, Fig1, Fig2, Fig3, Fig4, Fig5):
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
            plt.plot(self._BaseXDiff, self._BaseYDiff, 'bo', label='Desired End Point')
            plt.plot(0, 0, 'ro', label='Start Point')
            plt.plot(self._XDirDump, self._YDirDump, 'b-', label='Flight Path')
            plt.plot(self._XDirDump[-1], self._YDirDump[-1], 'go', label='End Point')
            plt.title('Predicted Flight Path of UAV')
            plt.xlabel('Distance in X Direction/m')
            plt.ylabel('Distance in Y Direction/m')
            plt.legend(loc=1)

        if Fig5 is True:
            fig5 = plt.figure(5)
            axes5 = fig5.gca(projection='3d')
            axes5.plot([self._BaseXDiff], [self._BaseYDiff], [0], 'bo', label='Desired End Point')
            axes5.plot([0], [0], [self.Alt], 'ro', label='Start Point')
            axes5.plot(self._XDirDump, self._YDirDump, self._AltDump, label='Flight Path')
            axes5.plot([self._XDirDump[-1]], [self._YDirDump[-1]], [self._AltDump[-1]], 'go', label='End Point')
            axes5.set_title('Predicted Flight Path of UAV in 3D')
            axes5.set_xlabel('Distance in X Direction/m')
            axes5.set_ylabel('Distance in Y Direction/m')
            axes5.set_zlabel('Distance in Z Direction/m')
            axes5.legend()

        # plt.figure(6)
        # plt.plot(self._TimeDump, self._XZDump, 'r', label='Forward Velocity')
        # plt.title('Forward Veloity Against Time')
        # plt.xlabel('Time/s')
        # plt.ylabel('Forward Velocity/m/s')
        # plt.legend()
        plt.show()
