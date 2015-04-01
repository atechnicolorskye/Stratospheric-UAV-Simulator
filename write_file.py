'''
write_file.py
WriteFile

Writes inputs and outputs into Python file for later processsing
'''

import os.path as os
import fileinput


def WriteToFile_Dist(StartTime, x_Dist, y_Dist, D_Ideal, Dist):
    '''
    Captures inputs from Fly methods and outputs them into a Python file for easy processing.
    For Distance-based algorithm.
    '''
    # Convert inputs to strings
    filename = 'plot_fly_range_' + str(StartTime).replace('(', '').replace(')', '').replace(', ', '_') + '.py'
    x_Dist, y_Dist, D_Ideal, Dist = str(x_Dist), str(y_Dist), str(D_Ideal), str(Dist)

    # Checks if file exists
    if os.isfile(filename):
        # Appends inputs to existing file
        for line in fileinput.input(filename, inplace=True):
            if 'x_Dist = ' in line:
                if 'x_Dist = []' in line:
                    print line.replace(']', x_Dist + ']'),
                else:
                    print line.replace(']', ', ' + x_Dist + ']'),
            elif 'y_Dist = ' in line:
                if 'y_Dist = []' in line:
                    print line.replace(']', y_Dist + ']'),
                else:
                    print line.replace(']', ', ' + y_Dist + ']'),
            elif 'Dist_Ideal = ' in line:
                if 'Dist_Ideal = []' in line:
                    print line.replace(']', D_Ideal + ']'),
                else:
                    print line.replace(']', ', ' + D_Ideal + ']'),
            elif 'Dist_Actual = ' in line:
                if 'Dist_Actual = []' in line:
                    print line.replace(']', Dist + ']'),
                else:
                    print line.replace(']', ', ' + Dist + ']'),
            else:
                print line,
    else:
        # Creates new file and appends inputs
        with open(filename, 'w') as x:
            x.write('import matplotlib.pyplot as plt\n\n')

            # Writes required lists for processing
            x.write('x_Dist = ' + '[' + x_Dist + ']' + '\n')
            x.write('y_Dist = ' + '[' + y_Dist + ']' + '\n')
            x.write('Dist_Ideal = ' + '[' + D_Ideal + ']' + '\n')
            x.write('Dist_Actual = ' + '[' + Dist + ']' + '\n')

            x.write('x_Head1 = ' + '[]' + '\n')
            x.write('y_Head1 = ' + '[]' + '\n')
            x.write('Head1_Ideal = ' + '[]' + '\n')
            x.write('Head1_Actual = ' + '[]' + '\n')

            x.write('x_Head2 = ' + '[]' + '\n')
            x.write('y_Head2 = ' + '[]' + '\n')
            x.write('Head2_Ideal = ' + '[]' + '\n')
            x.write('Head2_Actual = ' + '[]' + '\n')

            # Writes required code for sorting lists
            x.write('x_Dist_p = [x for (y, x) in sorted(zip(Dist_Ideal, x_Dist))]' + '\n')
            x.write('y_Dist_p = [x for (y, x) in sorted(zip(Dist_Ideal, y_Dist))]' + '\n')
            x.write('Dist_Actual_p = [x for (y, x) in sorted(zip(Dist_Ideal, Dist_Actual))]' + '\n')
            x.write('if Dist_Ideal.sort() is None:' + '\n')
            x.write('    Dist_Ideal_p = Dist_Ideal' + '\n')
            x.write('else:' + '\n')
            x.write('    Dist_Ideal_p = Dist_Ideal.sort()' + '\n')
            x.write('Dist_Diff = [y - x for (y, x) in zip(Dist_Actual_p, Dist_Ideal_p)]' + '\n\n')

            x.write('x_Head1_p = [x for (y, x) in sorted(zip(Head1_Ideal, x_Head1))]' + '\n')
            x.write('y_Head1_p = [x for (y, x) in sorted(zip(Head1_Ideal, y_Head1))]' + '\n')
            x.write('Head1_Actual_p = [x for (y, x) in sorted(zip(Head1_Ideal, Head1_Actual))]' + '\n')
            x.write('if Head1_Ideal.sort() is None:' + '\n')
            x.write('    Head1_Ideal_p = Head1_Ideal' + '\n')
            x.write('else:' + '\n')
            x.write('    Head1_Ideal_p = Head1_Ideal.sort()' + '\n')
            x.write('Head1_Diff = [y - x for (y, x) in zip(Head1_Actual_p, Head1_Ideal_p)]' + '\n\n')

            x.write('x_Head2_p = [x for (y, x) in sorted(zip(Head2_Ideal, x_Head2))]' + '\n')
            x.write('y_Head2_p = [x for (y, x) in sorted(zip(Head2_Ideal, y_Head2))]' + '\n')
            x.write('Head2_Actual_p = [x for (y, x) in sorted(zip(Head2_Ideal, Head2_Actual))]' + '\n')
            x.write('if Head2_Ideal.sort() is None:' + '\n')
            x.write('    Head2_Ideal_p = Head2_Ideal' + '\n')
            x.write('else:' + '\n')
            x.write('    Head2_Ideal_p = Head2_Ideal.sort()' + '\n')
            x.write('Head2_Diff = [y - x for (y, x) in zip(Head2_Actual_p, Head2_Ideal_p)]' + '\n\n')

            x.write('x_360 = range(0, 360)\n')
            x.write('x_0 = [0 for x in x_360]\n\n\n')

            # Writes required code to plot contour
            x.write('def PlotContours():\n    plt.figure(1)\n    plt.plot(x_Dist_p, y_Dist_p, marker=\'o\', color=\'r\', label=\'Fly_1\')\n    plt.fill(x_Dist_p, y_Dist_p, \'r\', alpha=0.5)\n    plt.plot(x_Head1_p, y_Head1_p, marker=\'o\', color=\'g\', label=\'Fly_2\')\n    plt.fill(x_Head1_p, y_Head1_p, \'g\', alpha=0.5)\n    plt.plot(x_Head2_p, y_Head2_p, marker=\'o\', color=\'b\', label=\'Fly_3\')\n    plt.fill(x_Head2_p, y_Head2_p, \'b\', alpha=0.5)\n    plt.plot(0, 0, \'wo\', label=\'Start Point\')\n    plt.title(\'Predicted Flight Path of UAV\')\n    plt.xlabel(\'Distance in X Direction/m\')\n    plt.ylabel(\'Distance in Y Direction/m\')\n    plt.legend()\n\n    plt.figure(2)\n    plt.plot(x_360, x_0, color=\'k\', linewidth=2.0, label=\'Ideal\')\n    plt.plot(Dist_Ideal_p, Dist_Diff, color=\'r\', linewidth=2.0, label=\'Fly_1\')\n    plt.plot(Head1_Ideal_p, Head1_Diff, color=\'g\', linewidth=2.0, label=\'Fly_2\')\n    plt.plot(Head2_Ideal_p, Head2_Diff, color=\'b\', linewidth=2.0, label=\'Fly_3\')\n    plt.title(\'Deviations from Desired Directions\')\n    plt.xlabel(\'Desired Direction\')\n    plt.ylabel(\'Deviation\')\n    plt.xlim(0, 359)\n    plt.legend()\n\n    plt.show()\n\n\n')

            # Writes required code to display average error
            x.write('def Error():\n    Mean_Dist_Error_List = [abs(a - b) for a, b in zip(Dist_Ideal_p, Dist_Actual_p)]\n    Mean_Dist_Error = sum(Mean_Dist_Error_List) / len(Dist_Ideal_p)\n\n    Mean_Head1_Error_List = [abs(a - b) for a, b in zip(Head1_Ideal_p, Head1_Actual_p)]\n    Mean_Head1_Error = sum(Mean_Head1_Error_List) / len(Head1_Ideal_p)\n\n    Mean_Head2_Error_List = [abs(a - b) for a, b in zip(Head2_Ideal_p, Head2_Actual_p)]\n    Mean_Head2_Error = sum(Mean_Head2_Error_List) / len(Head2_Ideal_p)\n    print Mean_Dist_Error, Mean_Head1_Error, Mean_Head2_Error\n\n')

            # Writes required code to run code
            x.write('PlotContours()\n\n')
            x.write('Error()')

    print filename
    print 'Write Completed'


def WriteToFile_Head1(StartTime, x_Dist, y_Dist, D_Ideal, Dist):
    '''
    Captures inputs from Fly methods and outputs them into a Python file for easy processing.
    For Distance-based algorithm.
    '''
    # Convert inputs to strings
    filename = 'plot_fly_range_' + str(StartTime).replace('(', '').replace(')', '').replace(', ', '_') + '.py'
    x_Dist, y_Dist, D_Ideal, Dist = str(x_Dist), str(y_Dist), str(D_Ideal), str(Dist)

    # Checks if file exists
    if os.isfile(filename):
        # Appends inputs to existing file
        for line in fileinput.input(filename, inplace=True):
            if 'x_Head1 = ' in line:
                if 'x_Head1 = []' in line:
                    print line.replace(']', x_Dist + ']'),
                else:
                    print line.replace(']', ', ' + x_Dist + ']'),
            elif 'y_Head1 = ' in line:
                if 'y_Head1 = []' in line:
                    print line.replace(']', y_Dist + ']'),
                else:
                    print line.replace(']', ', ' + y_Dist + ']'),
            elif 'Head1_Ideal = ' in line:
                if 'Head1_Ideal = []' in line:
                    print line.replace(']', D_Ideal + ']'),
                else:
                    print line.replace(']', ', ' + D_Ideal + ']'),
            elif 'Head1_Actual = ' in line:
                if 'Head1_Actual = []' in line:
                    print line.replace(']', Dist + ']'),
                else:
                    print line.replace(']', ', ' + Dist + ']'),
            else:
                print line,

        # x = open(filename, 'rb+')
        # x.seek(0, 0)
        # for line in x:
        #     if 'x_Dist = ' in line:
        #         line = line.replace(']', ', ' + x_Dist[1:])
        #     elif 'y_Dist = ' in line:
        #         line = line.replace(']', ', ' + y_Dist[1:])
        #     elif 'Dist_Ideal = ' in line:
        #         line = line.replace(']', ', ' + D_Ideal[1:])
        #     elif 'Dist_Actual = ' in line:
        #         line = line.replace(']', ', ' + Dist[1:])
        # x.readlines()
        # x.close()

    else:
        # Creates new file and appends inputs
        with open(filename, 'w') as x:
            x.write('import matplotlib.pyplot as plt\n\n')

            # Writes required lists for processing
            x.write('x_Dist = ' + '[]' + '\n')
            x.write('y_Dist = ' + '[]' + '\n')
            x.write('Dist_Ideal = ' + '[]' + '\n')
            x.write('Dist_Actual = ' + '[]' + '\n')

            x.write('x_Head1 = ' + '[' + x_Dist + ']' + '\n')
            x.write('y_Head1 = ' + '[' + y_Dist + ']' + '\n')
            x.write('Head1_Ideal = ' + '[' + D_Ideal + ']' + '\n')
            x.write('Head1_Actual = ' + '[' + Dist + ']' + '\n')

            x.write('x_Head2 = ' + '[]' + '\n')
            x.write('y_Head2 = ' + '[]' + '\n')
            x.write('Head2_Ideal = ' + '[]' + '\n')
            x.write('Head2_Actual = ' + '[]' + '\n')

            # Writes required code for sorting lists
            x.write('x_Dist_p = [x for (y, x) in sorted(zip(Dist_Ideal, x_Dist))]' + '\n')
            x.write('y_Dist_p = [x for (y, x) in sorted(zip(Dist_Ideal, y_Dist))]' + '\n')
            x.write('Dist_Actual_p = [x for (y, x) in sorted(zip(Dist_Ideal, Dist_Actual))]' + '\n')
            x.write('if Dist_Ideal.sort() is None:' + '\n')
            x.write('    Dist_Ideal_p = Dist_Ideal' + '\n')
            x.write('else:' + '\n')
            x.write('    Dist_Ideal_p = Dist_Ideal.sort()' + '\n')
            x.write('Dist_Diff = [y - x for (y, x) in zip(Dist_Actual_p, Dist_Ideal_p)]' + '\n\n')

            x.write('x_Head1_p = [x for (y, x) in sorted(zip(Head1_Ideal, x_Head1))]' + '\n')
            x.write('y_Head1_p = [x for (y, x) in sorted(zip(Head1_Ideal, y_Head1))]' + '\n')
            x.write('Head1_Actual_p = [x for (y, x) in sorted(zip(Head1_Ideal, Head1_Actual))]' + '\n')
            x.write('if Head1_Ideal.sort() is None:' + '\n')
            x.write('    Head1_Ideal_p = Head1_Ideal' + '\n')
            x.write('else:' + '\n')
            x.write('    Head1_Ideal_p = Head1_Ideal.sort()' + '\n')
            x.write('Head1_Diff = [y - x for (y, x) in zip(Head1_Actual_p, Head1_Ideal_p)]' + '\n\n')

            x.write('x_Head2_p = [x for (y, x) in sorted(zip(Head2_Ideal, x_Head2))]' + '\n')
            x.write('y_Head2_p = [x for (y, x) in sorted(zip(Head2_Ideal, y_Head2))]' + '\n')
            x.write('Head2_Actual_p = [x for (y, x) in sorted(zip(Head2_Ideal, Head2_Actual))]' + '\n')
            x.write('if Head2_Ideal.sort() is None:' + '\n')
            x.write('    Head2_Ideal_p = Head2_Ideal' + '\n')
            x.write('else:' + '\n')
            x.write('    Head2_Ideal_p = Head2_Ideal.sort()' + '\n')
            x.write('Head2_Diff = [y - x for (y, x) in zip(Head2_Actual_p, Head2_Ideal_p)]' + '\n\n')

            x.write('x_360 = range(0, 360)\n')
            x.write('x_0 = [0 for x in x_360]\n\n\n')

            # Writes required code to plot contour
            x.write('def PlotContours():\n    plt.figure(1)\n    plt.plot(x_Dist_p, y_Dist_p, marker=\'o\', color=\'r\', label=\'Fly_1\')\n    plt.fill(x_Dist_p, y_Dist_p, \'r\', alpha=0.5)\n    plt.plot(x_Head1_p, y_Head1_p, marker=\'o\', color=\'g\', label=\'Fly_2\')\n    plt.fill(x_Head1_p, y_Head1_p, \'g\', alpha=0.5)\n    plt.plot(x_Head2_p, y_Head2_p, marker=\'o\', color=\'b\', label=\'Fly_3\')\n    plt.fill(x_Head2_p, y_Head2_p, \'b\', alpha=0.5)\n    plt.plot(0, 0, \'wo\', label=\'Start Point\')\n    plt.title(\'Predicted Flight Path of UAV\')\n    plt.xlabel(\'Distance in X Direction/m\')\n    plt.ylabel(\'Distance in Y Direction/m\')\n    plt.legend()\n\n    plt.figure(2)\n    plt.plot(x_360, x_0, color=\'k\', linewidth=2.0, label=\'Ideal\')\n    plt.plot(Dist_Ideal_p, Dist_Diff, color=\'r\', linewidth=2.0, label=\'Fly_1\')\n    plt.plot(Head1_Ideal_p, Head1_Diff, color=\'g\', linewidth=2.0, label=\'Fly_2\')\n    plt.plot(Head2_Ideal_p, Head2_Diff, color=\'b\', linewidth=2.0, label=\'Fly_3\')\n    plt.title(\'Deviations from Desired Directions\')\n    plt.xlabel(\'Desired Direction\')\n    plt.ylabel(\'Deviation\')\n    plt.xlim(0, 359)\n    plt.legend()\n\n    plt.show()\n\n\n')

            # Writes required code to display average error
            x.write('def Error():\n    Mean_Dist_Error_List = [abs(a - b) for a, b in zip(Dist_Ideal_p, Dist_Actual_p)]\n    Mean_Dist_Error = sum(Mean_Dist_Error_List) / len(Dist_Ideal_p)\n\n    Mean_Head1_Error_List = [abs(a - b) for a, b in zip(Head1_Ideal_p, Head1_Actual_p)]\n    Mean_Head1_Error = sum(Mean_Head1_Error_List) / len(Head1_Ideal_p)\n\n    Mean_Head2_Error_List = [abs(a - b) for a, b in zip(Head2_Ideal_p, Head2_Actual_p)]\n    Mean_Head2_Error = sum(Mean_Head2_Error_List) / len(Head2_Ideal_p)\n    print Mean_Dist_Error, Mean_Head1_Error, Mean_Head2_Error\n\n')

            # Writes required code to run code
            x.write('PlotContours()\n\n')
            x.write('Error()')

    print filename
    print 'Write Completed'


def WriteToFile_Head2(StartTime, x_Dist, y_Dist, D_Ideal, Dist):
    '''
    Captures inputs from Fly methods and outputs them into a Python file for easy processing.
    For Distance-based algorithm.
    '''
    # Convert inputs to strings
    filename = 'plot_fly_range_' + str(StartTime).replace('(', '').replace(')', '').replace(', ', '_') + '.py'
    x_Dist, y_Dist, D_Ideal, Dist = str(x_Dist), str(y_Dist), str(D_Ideal), str(Dist)

    # Checks if file exists
    if os.isfile(filename):
        # Appends inputs to existing file
        for line in fileinput.input(filename, inplace=True):
            if 'x_Head2 = ' in line:
                if 'x_Head2 = []' in line:
                    print line.replace(']', x_Dist + ']'),
                else:
                    print line.replace(']', ', ' + x_Dist + ']'),
            elif 'y_Head2 = ' in line:
                if 'y_Head2 = []' in line:
                    print line.replace(']', y_Dist + ']'),
                else:
                    print line.replace(']', ', ' + y_Dist + ']'),
            elif 'Head2_Ideal = ' in line:
                if 'Head2_Ideal = []' in line:
                    print line.replace(']', D_Ideal + ']'),
                else:
                    print line.replace(']', ', ' + D_Ideal + ']'),
            elif 'Head2_Actual = ' in line:
                if 'Head2_Actual = []' in line:
                    print line.replace(']', Dist + ']'),
                else:
                    print line.replace(']', ', ' + Dist + ']'),
            else:
                print line,
    else:
        # Creates new file and appends inputs
        with open(filename, 'w') as x:
            x.write('import matplotlib.pyplot as plt\n\n')

            # Writes required lists for processing
            x.write('x_Dist = ' + '[]' + '\n')
            x.write('y_Dist = ' + '[]' + '\n')
            x.write('Dist_Ideal = ' + '[]' + '\n')
            x.write('Dist_Actual = ' + '[]' + '\n')

            x.write('x_Head1 = ' + '[]' + '\n')
            x.write('y_Head1 = ' + '[]' + '\n')
            x.write('Head1_Ideal = ' + '[]' + '\n')
            x.write('Head1_Actual = ' + '[]' + '\n')

            x.write('x_Head2 = ' + '[' + x_Dist + ']' + '\n')
            x.write('y_Head2 = ' + '[' + y_Dist + ']' + '\n')
            x.write('Head2_Ideal = ' + '[' + D_Ideal + ']' + '\n')
            x.write('Head2_Actual = ' + '[' + Dist + ']' + '\n')

            # Writes required code for sorting lists
            x.write('x_Dist_p = [x for (y, x) in sorted(zip(Dist_Ideal, x_Dist))]' + '\n')
            x.write('y_Dist_p = [x for (y, x) in sorted(zip(Dist_Ideal, y_Dist))]' + '\n')
            x.write('Dist_Actual_p = [x for (y, x) in sorted(zip(Dist_Ideal, Dist_Actual))]' + '\n')
            x.write('if Dist_Ideal.sort() is None:' + '\n')
            x.write('    Dist_Ideal_p = Dist_Ideal' + '\n')
            x.write('else:' + '\n')
            x.write('    Dist_Ideal_p = Dist_Ideal.sort()' + '\n')
            x.write('Dist_Diff = [y - x for (y, x) in zip(Dist_Actual_p, Dist_Ideal_p)]' + '\n\n')

            x.write('x_Head1_p = [x for (y, x) in sorted(zip(Head1_Ideal, x_Head1))]' + '\n')
            x.write('y_Head1_p = [x for (y, x) in sorted(zip(Head1_Ideal, y_Head1))]' + '\n')
            x.write('Head1_Actual_p = [x for (y, x) in sorted(zip(Head1_Ideal, Head1_Actual))]' + '\n')
            x.write('if Head1_Ideal.sort() is None:' + '\n')
            x.write('    Head1_Ideal_p = Head1_Ideal' + '\n')
            x.write('else:' + '\n')
            x.write('    Head1_Ideal_p = Head1_Ideal.sort()' + '\n')
            x.write('Head1_Diff = [y - x for (y, x) in zip(Head1_Actual_p, Head1_Ideal_p)]' + '\n\n')

            x.write('x_Head2_p = [x for (y, x) in sorted(zip(Head2_Ideal, x_Head2))]' + '\n')
            x.write('y_Head2_p = [x for (y, x) in sorted(zip(Head2_Ideal, y_Head2))]' + '\n')
            x.write('Head2_Actual_p = [x for (y, x) in sorted(zip(Head2_Ideal, Head2_Actual))]' + '\n')
            x.write('if Head2_Ideal.sort() is None:' + '\n')
            x.write('    Head2_Ideal_p = Head2_Ideal' + '\n')
            x.write('else:' + '\n')
            x.write('    Head2_Ideal_p = Head2_Ideal.sort()' + '\n')
            x.write('Head2_Diff = [y - x for (y, x) in zip(Head2_Actual_p, Head2_Ideal_p)]' + '\n\n')

            x.write('x_360 = range(0, 360)\n')
            x.write('x_0 = [0 for x in x_360]\n\n\n')

            # Writes required code to plot contour
            x.write('def PlotContours():\n    plt.figure(1)\n    plt.plot(x_Dist_p, y_Dist_p, marker=\'o\', color=\'r\', label=\'Fly_1\')\n    plt.fill(x_Dist_p, y_Dist_p, \'r\', alpha=0.5)\n    plt.plot(x_Head1_p, y_Head1_p, marker=\'o\', color=\'g\', label=\'Fly_2\')\n    plt.fill(x_Head1_p, y_Head1_p, \'g\', alpha=0.5)\n    plt.plot(x_Head2_p, y_Head2_p, marker=\'o\', color=\'b\', label=\'Fly_3\')\n    plt.fill(x_Head2_p, y_Head2_p, \'b\', alpha=0.5)\n    plt.plot(0, 0, \'wo\', label=\'Start Point\')\n    plt.title(\'Predicted Flight Path of UAV\')\n    plt.xlabel(\'Distance in X Direction/m\')\n    plt.ylabel(\'Distance in Y Direction/m\')\n    plt.legend()\n\n    plt.figure(2)\n    plt.plot(x_360, x_0, color=\'k\', linewidth=2.0, label=\'Ideal\')\n    plt.plot(Dist_Ideal_p, Dist_Diff, color=\'r\', linewidth=2.0, label=\'Fly_1\')\n    plt.plot(Head1_Ideal_p, Head1_Diff, color=\'g\', linewidth=2.0, label=\'Fly_2\')\n    plt.plot(Head2_Ideal_p, Head2_Diff, color=\'b\', linewidth=2.0, label=\'Fly_3\')\n    plt.title(\'Deviations from Desired Directions\')\n    plt.xlabel(\'Desired Direction\')\n    plt.ylabel(\'Deviation\')\n    plt.xlim(0, 359)\n    plt.legend()\n\n    plt.show()\n\n\n')

            # Writes required code to display average error
            x.write('def Error():\n    Mean_Dist_Error_List = [abs(a - b) for a, b in zip(Dist_Ideal_p, Dist_Actual_p)]\n    Mean_Dist_Error = sum(Mean_Dist_Error_List) / len(Dist_Ideal_p)\n\n    Mean_Head1_Error_List = [abs(a - b) for a, b in zip(Head1_Ideal_p, Head1_Actual_p)]\n    Mean_Head1_Error = sum(Mean_Head1_Error_List) / len(Head1_Ideal_p)\n\n    Mean_Head2_Error_List = [abs(a - b) for a, b in zip(Head2_Ideal_p, Head2_Actual_p)]\n    Mean_Head2_Error = sum(Mean_Head2_Error_List) / len(Head2_Ideal_p)\n    print Mean_Dist_Error, Mean_Head1_Error, Mean_Head2_Error\n\n')

            # Writes required code to run code
            x.write('PlotContours()\n\n')
            x.write('Error()')

    print filename
    print 'Write Completed'


def WriteToFile_Local(StartTime, x_Dist, y_Dist, x_Dist_Max, y_Dist_Max, x_Dist_Ideal, y_Dist_Ideal, D_Ideal, Dist, x_Dist_Fail, y_Dist_Fail, x_Dist_Ideal_Fail, y_Dist_Ideal_Fail, D_Ideal_Fail, Dist_Fail, Error):
    '''
    Captures inputs from Fly methods and outputs them into a Python file for easy processing.
    For Distance-based algorithm.
    '''
    # Convert inputs to strings
    filename = 'plot_fly_range_2_' + str(StartTime).replace('(', '').replace(')', '').replace(', ', '_') + '_' + str(Error).replace('.', '_') + '.py'
    x_Dist, y_Dist, x_Dist_Max, y_Dist_Max, x_Dist_Ideal, y_Dist_Ideal, D_Ideal, Dist = str(x_Dist), str(y_Dist), str(x_Dist_Max), str(y_Dist_Max), str(x_Dist_Ideal), str(y_Dist_Ideal), str(D_Ideal), str(Dist)
    x_Dist_Fail, y_Dist_Fail, x_Dist_Ideal_Fail, y_Dist_Ideal_Fail, D_Ideal_Fail, Dist_Fail = str(x_Dist_Fail), str(y_Dist_Fail), str(x_Dist_Ideal_Fail), str(y_Dist_Ideal_Fail), str(D_Ideal_Fail), str(Dist_Fail)

    # Creates new file and appends inputs
    with open(filename, 'w') as x:
        x.write('import matplotlib.pyplot as plt\n\n')

        # Writes required lists for processing
        x.write('x_Dist = ' + x_Dist + '\n')
        x.write('y_Dist = ' + y_Dist + '\n')
        x.write('x_Dist_Max = ' + x_Dist_Max + '\n')
        x.write('y_Dist_Max = ' + y_Dist_Max + '\n')
        x.write('x_Dist_Ideal = ' + x_Dist_Ideal + '\n')
        x.write('y_Dist_Ideal = ' + y_Dist_Ideal + '\n')
        x.write('Dist_Ideal = ' + D_Ideal + '\n')
        x.write('Dist_Actual = ' + Dist + '\n\n')
        x.write('x_Dist_Fail = ' + x_Dist_Fail + '\n')
        x.write('y_Dist_Fail = ' + y_Dist_Fail + '\n')
        x.write('x_Dist_Ideal_Fail = ' + x_Dist_Ideal_Fail + '\n')
        x.write('y_Dist_Ideal_Fail = ' + y_Dist_Ideal_Fail + '\n')
        x.write('Dist_Ideal_Fail = ' + D_Ideal_Fail + '\n')
        x.write('Dist_Actual_Fail = ' + Dist_Fail + '\n\n')

        x.write('plt.plot(x_Dist_Max, y_Dist_Max, \'go\', linestyle=\'-\', label=\'Flight Contour\')\n')
        x.write('plt.fill(x_Dist_Max, y_Dist_Max, \'g\', alpha=0.5)\n')
        x.write('plt.plot(x_Dist, y_Dist, \'bo\', label=\'Set Points Reached\')\n')
        x.write('plt.plot(x_Dist_Ideal_Fail, y_Dist_Ideal_Fail, \'ro\', label=\'Set Points Failed to Reach\')\n')
        x.write('plt.plot(0, 0, \'ko\', label=\'Start Point\')\n')
        x.write('plt.title(\'Flight Contour\')\n')
        x.write('plt.xlabel(\'UAV Displacement in X\')\n')
        x.write('plt.ylabel(\'UAV Displacement in Y\')\n')
        x.write('plt.xlim(-5000 - max(x_Dist_Ideal), max(x_Dist_Ideal) + 5000)\n')
        x.write('plt.ylim(-5000 - max(x_Dist_Ideal), max(x_Dist_Ideal) + 5000)\n')
        x.write('plt.legend()\n')
        x.write('plt.show()')

    print filename
    print 'Write Completed'


def WriteToFile_Simulator(StartTime, x_Dist, y_Dist, Density, Altitude, Time, TimeStep1):
    '''
    Captures inputs from Fly methods and outputs them into a Python file for easy processing.
    For Distance-based algorithm.
    '''
    # Convert inputs to string
    filename = 'plot_flight_' + str(StartTime).replace('(', '').replace(')', '').replace(', ', '_') + '_' + str(TimeStep1).replace('.', '_') + '.py'
    x_Dist, y_Dist, Density, Altitude, Time = str(x_Dist), str(y_Dist), str(Density), str(Altitude), str(Time)

    with open(filename, 'w') as x:
        x.write('x_Dist = ' + x_Dist + '\n')
        x.write('y_Dist = ' + y_Dist + '\n')
        x.write('Density = ' + Density + '\n')
        x.write('Altitude = ' + Altitude + '\n')
        x.write('Time = ' + Time + '\n')

    print filename
    print 'Write Completed'
