from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import xlrd
import pandas as pd
import numpy as np


################################################# Curve fitting ########################################################
def fitLine(data, xUnit, yUnit):
    c1 = c2 = c3 = c4 = c5 = 0

    c5 = len(data)
    for i in range(len(data)):
        c1 = c1 + (data.at[i, xUnit] * data.at[i, xUnit])
        c2 = c2 + data.at[i, xUnit]
        c3 = c3 + (data.at[i, xUnit] * data.at[i, yUnit])
        c4 = c4 + data.at[i, yUnit]

    c = ((c3 * c2) - (c1 * c4)) / ((c2 * c2) - (c1 * c5))
    m = (c3 - (c * c2)) / c1

    return m, c


def fitPlane(data, xUnit, yUnit, zUnit):
    a1 = a2 = a3 = a4 = a5 = a6 = 0
    b1 = b2 = b3 = 0

    a6 = len(data)
    for i in range(len(data)):
        a1 = a1 + (data.at[i, xUnit] * data.at[i, xUnit])
        a2 = a2 + (data.at[i, xUnit] * data.at[i, yUnit])
        a3 = a3 + data.at[i, xUnit]
        a4 = a4 + (data.at[i, yUnit] * data.at[i, yUnit])
        a5 = a5 + data.at[i, yUnit]

        b1 = b1 + (data.at[i, zUnit] * data.at[i, xUnit])
        b2 = b2 + (data.at[i, zUnit] * data.at[i, yUnit])
        b3 = b3 + data.at[i, zUnit]

    coefficientMatrix = np.array([[a1, a2, a3],
                                  [a2, a4, a5],
                                  [a3, a5, a6]])
    constantMatrix = np.array([b1, b2, b3])

    xCoeff, yCoeff, const = solveLinearSystem(coefficientMatrix, constantMatrix)

    writeCoeffInText(xCoeff, yCoeff, const)


########################################################################################################################


############################################ Equations solver ##########################################################
def solveLinearSystem(coefficientMatrix, constantMatrix):
    coefficientMatrixInverse = np.linalg.inv(coefficientMatrix)

    return np.dot(coefficientMatrixInverse, constantMatrix)


########################################################################################################################


################################################## Prediction ##########################################################
def predict(NCData, toPredict):
    NCWtMeans, NCBGMeans = getNCMeans(NCData)

    xCoeff, yCoeff, const = readCoeffFromText()
    diseasePercentage = []

    for i in range(len(toPredict)):

        if toPredict.at[i, getDayValue()] == 0:
            deltaBodyWeight = float(NCWtMeans[0]) - toPredict.at[i, getBodyWeightValue()]
            deltaBloodGlucose = float(NCBGMeans[0]) - toPredict.at[i, getBloodGlucoseValue()]
        else:
            j = int(toPredict.at[i, getDayValue()] / 7)

            deltaBodyWeight = (float(NCWtMeans[j]) - float(NCWtMeans[j - 1])) - \
                              (toPredict.at[j, getBodyWeightValue()] - toPredict.at[j - 1, getBodyWeightValue()])
            deltaBloodGlucose = toPredict.at[j, getBloodGlucoseValue()] - float(NCBGMeans[j])

        percentage = (xCoeff * deltaBodyWeight) + (yCoeff * deltaBloodGlucose) + const
        if percentage < 0:
            percentage = 0
        elif percentage > 100:
            percentage = 100

        diseasePercentage.append([toPredict.at[i, 'Type'], toPredict.at[i, 'Mouse'], toPredict.at[i, getDayValue()],
                                      percentage])

    print(diseasePercentage)
    diseasePercentDataFrame = pd.DataFrame(diseasePercentage, columns = ['Type', 'Mouse', 'Day', 'Percentage'])

    return diseasePercentDataFrame


########################################################################################################################


############################################ Curve Surface Plotter #####################################################
def drawLineAndScatterPlot(data, m, c, xUnit, yUnit):
    plt.scatter(data[xUnit], data[yUnit], color="m", marker="o", s=30)

    yPred = m * data[xUnit] + c

    plt.plot(data[xUnit], yPred, color="g")

    plt.xlabel('x')
    plt.ylabel('y')

    plt.show()


def drawPlaneAndScatterPlot(data, xUnit, yUnit, zUnit):
    xCoeff, yCoeff, const = readCoeffFromText()

    def z_function(x, y):
        return (xCoeff * x) + (yCoeff * y) + const

    x = np.linspace(-110, 180, 70)
    y = np.linspace(100, 400, 70)

    X, Y = np.meshgrid(x, y)
    Z = z_function(X, Y)

    ax = plt.axes(projection = '3d')

    ax.scatter(data[xUnit], data[yUnit], data[zUnit], color = 'm', marker = 'o', s = 30)

    ax.plot_wireframe(X, Y, Z, color = 'green')
    ax.set_xlabel('Delta BodyWeight')
    ax.set_ylabel('Delta BloodGlucose')
    ax.set_zlabel('Disease Percentage')

    plt.show()


def plotDiseaseProgression(data, mouseNumber):
    tempData = getByMouse(data, mouseNumber)

    plt.scatter(tempData[getDayValue()], tempData[getPercentageValue()], color='m', marker='o', s=30)

    plt.plot(tempData[getDayValue()], tempData[getPercentageValue()])

    plt.title(tempData.iat[0, 0] + ' mouse number ' + str(mouseNumber))  # something to be done about this
    plt.xlabel(getDayValue())
    plt.ylabel('Disease percentage')
    plt.grid(True)

    plt.show()


def drawSubPlotToCompare(STD, test1, test2, test3):
    STDTempData = getByMouse(STD, 13)
    test1TempData = getByMouse(test1, 19)
    test2TempData = getByMouse(test2, 25)
    test3TempData = getByMouse(test3, 31)

    plt.subplot(2, 2, 1)
    plt.scatter(STDTempData[getDayValue()], STDTempData[getPercentageValue()], color='m', marker='o', s=30)
    plt.plot(STDTempData[getDayValue()], STDTempData[getPercentageValue()])
    plt.title('Disease progression under standard drug for rat number = ' + str(13))  # something to be done about this
    plt.ylabel('Disease percentage')
    plt.grid(True)

    plt.subplot(2, 2, 2)
    plt.scatter(test1TempData[getDayValue()], test1TempData[getPercentageValue()], color='m', marker='o', s=30)
    plt.plot(test1TempData[getDayValue()], test1TempData[getPercentageValue()])
    plt.title('Disease progression under Test-1 drug for rat number = ' + str(19))  # something to be done about this
    plt.ylabel('Disease percentage')
    plt.grid(True)

    plt.subplot(2, 2, 3)
    plt.scatter(test2TempData[getDayValue()], test2TempData[getPercentageValue()], color='m', marker='o', s=30)
    plt.plot(test2TempData[getDayValue()], test2TempData[getPercentageValue()])
    plt.title('Disease progression under Test-2 drug for rat number = ' + str(25))  # something to be done about this
    plt.xlabel(getDayValue())
    plt.ylabel('Disease percentage')
    plt.grid(True)

    plt.subplot(2, 2, 4)
    plt.scatter(test3TempData[getDayValue()], test3TempData[getPercentageValue()], color='m', marker='o', s=30)
    plt.plot(test3TempData[getDayValue()], test3TempData[getPercentageValue()])
    plt.title('Disease progression under Test-3 drug for rat number = ' + str(31))  # something to be done about this
    plt.xlabel(getDayValue())
    plt.ylabel('Disease percentage')
    plt.grid(True)

    plt.show()


########################################################################################################################


############################################## Drawing scatterplots ####################################################
def draw2DScatterPlotTest(data, xUnit, yUnit):
    colors = ['green', 'orange', 'brown', 'blue', 'red', 'black']
    colorMap = []

    for i in range(6):
        for j in range(5):
            colorMap.append(colors[i])

    plt.scatter(data[xUnit], data[yUnit], color=colorMap[:], marker='o', s=30)

    plt.title(data.iat[0, 0])
    plt.xlabel(xUnit)
    plt.ylabel(yUnit)
    plt.grid(True)

    plt.show()


def draw2DScatterPlot(data, xUnit, yUnit):
    plt.scatter(data[xUnit], data[yUnit], color='m', marker='o', s=30)

    plt.title(data.iat[0, 0])  # something to be done about this
    plt.xlabel(xUnit)
    plt.ylabel(yUnit)
    plt.grid(True)

    plt.show()


def draw3DScatterPlot(data, xUnit, yUnit, zUnit):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(data[xUnit], data[yUnit], data[zUnit], color='m', marker='o', s=30)

    plt.title(data.iat[0, 0])
    ax.set_xlabel(xUnit)
    ax.set_ylabel(yUnit)
    ax.set_zlabel(zUnit)

    plt.show()


########################################################################################################################


#################################################### NCMeans ###########################################################
def getNCMeans(NCData):
    NCBloodGlucoseMean = []
    NCBodyWeightMean = []

    for i in range(0, 29, 7):
        NCTempData = getByDay(NCData, i)

        bloodGlucose = NCTempData[['BloodGlucose']].mean()
        NCBloodGlucoseMean.append(bloodGlucose)

        bodyWeight = NCTempData[['BodyWeight']].mean()
        NCBodyWeightMean.append(bodyWeight)

    return NCBodyWeightMean, NCBloodGlucoseMean


########################################################################################################################

################################################# Data massaging #######################################################
def massageData(NCData, DCData):
    tempData = DCData.copy()

    NCBodyWeightMean, NCBloodGlucoseMean = getNCMeans(NCData)

    for i in range(len(DCData)):
        if tempData.at[i, 'Day'] == 0:
            tempData.at[i, 'BloodGlucose'] = DCData.at[i, 'BloodGlucose'] - NCBloodGlucoseMean[0]
            tempData.at[i, 'BodyWeight'] = NCBodyWeightMean[0] - DCData.at[i, 'BodyWeight']

        else:
            j = int(tempData.at[i, 'Day'] / 7)
            tempData.at[i, 'BloodGlucose'] = DCData.at[i, 'BloodGlucose'] - NCBloodGlucoseMean[j]
            tempData.at[i, 'BodyWeight'] = (NCBodyWeightMean[j] - NCBodyWeightMean[j - 1]) - \
                                           (DCData.at[i, 'BodyWeight'] - DCData.at[i - 1, 'BodyWeight'])

    return tempData


########################################################################################################################


################################################# File reading #########################################################
def readData():
    loc = ("E:\College\Projects\DiseaseProgression\Data\DiabetesData.xlsx")

    wb = xlrd.open_workbook(loc)
    bodyWeightSheet = wb.sheet_by_index(0)
    bloodGlucoseSheet = wb.sheet_by_index(1)
    percentageSheet = wb.sheet_by_index(2)

    NCData = []
    DCData = []
    STDData = []
    T1 = []
    T2 = []
    T3 = []

    for i in range(1, 37):
        day = 0

        if bodyWeightSheet.cell_value(i, 0) == 'NC':
            for j in range(2, 7):
                NCData.append(['NC', i, day, bodyWeightSheet.cell(i, j).value,
                               bloodGlucoseSheet.cell(i, j).value, percentageSheet.cell(i, j).value])
                day = day + 7

        elif bodyWeightSheet.cell_value(i, 0) == 'DC':
            for j in range(2, 7):
                DCData.append(['DC', i, day, bodyWeightSheet.cell(i, j).value,
                               bloodGlucoseSheet.cell(i, j).value, percentageSheet.cell(i, j).value])
                day = day + 7

        elif bodyWeightSheet.cell_value(i, 0) == 'STD':
            for j in range(2, 7):
                STDData.append(['STD', i, day, bodyWeightSheet.cell(i, j).value,
                               bloodGlucoseSheet.cell(i, j).value])
                day = day + 7

        elif bodyWeightSheet.cell_value(i, 0) == 'TEST-1':
            for j in range(2, 7):
                T1.append(['TEST-1', i, day, bodyWeightSheet.cell(i, j).value,
                               bloodGlucoseSheet.cell(i, j).value])
                day = day + 7

        elif bodyWeightSheet.cell_value(i, 0) == 'TEST-2':
            for j in range(2, 7):
                T2.append(['TEST-2', i, day, bodyWeightSheet.cell(i, j).value,
                               bloodGlucoseSheet.cell(i, j).value])
                day = day + 7

        elif bodyWeightSheet.cell_value(i, 0) == 'TEST-3':
            for j in range(2, 7):
                T3.append(['TEST-3', i, day, bodyWeightSheet.cell(i, j).value,
                               bloodGlucoseSheet.cell(i, j).value])
                day = day + 7

    NCDataFrame = pd.DataFrame(NCData, columns=['Type', 'Mouse', 'Day', 'BodyWeight',
                                                'BloodGlucose', 'Percentage'])
    DCDataFrame = pd.DataFrame(DCData, columns=['Type', 'Mouse', 'Day', 'BodyWeight',
                                                'BloodGlucose', 'Percentage'])
    STDDataFrame = pd.DataFrame(STDData, columns=['Type', 'Mouse', 'Day', 'BodyWeight',
                                                  'BloodGlucose'])
    T1DataFrame = pd.DataFrame(T1, columns=['Type', 'Mouse', 'Day', 'BodyWeight',
                                            'BloodGlucose'])
    T2DataFrame = pd.DataFrame(T2, columns=['Type', 'Mouse', 'Day', 'BodyWeight',
                                            'BloodGlucose'])
    T3DataFrame = pd.DataFrame(T3, columns=['Type', 'Mouse', 'Day', 'BodyWeight',
                                            'BloodGlucose'])

    return NCDataFrame, DCDataFrame, STDDataFrame, T1DataFrame, T2DataFrame, T3DataFrame


########################################################################################################################


############################################## Read write text file ####################################################
def writeCoeffInText(xCoeff, yCoeff, const):
    file = open('E:\College\Projects\DiseaseProgression\Generated Data\TrainedFunction\TrainedFunction.txt', 'w')

    file.write(str(xCoeff) + ' \n')
    file.write(str(yCoeff) + ' \n')
    file.write(str(const) + ' \n')

    file.close()


def readCoeffFromText():
    coeff = []

    with open('E:\College\Projects\DiseaseProgression\Generated Data\TrainedFunction\TrainedFunction.txt') as file:
        for line in file:
            coeff.append(float(line))
    file.close()

    xCoeff = coeff[0]
    yCoeff = coeff[1]
    const = coeff[2]

    return xCoeff, yCoeff, const


def writeNCMeansInText(NCWtMeans, NCBGMeans):
    file = open('E:\College\Projects\DiseaseProgression\Generated Data\Means\Means.txt', 'w')

    for i in range(len(NCWtMeans)):
        file.write(str(NCWtMeans[i]) + str(NCBGMeans[i]))

    file.close()


########################################################################################################################


############################################### Dataframe operations ###################################################
def joinDataFrames(data1, data2):
    tempData = pd.concat([data1, data2], ignore_index = True)

    return tempData


def getByMouse(data, mouse):
    tempData = data.query('Mouse == ' + str(mouse), inplace=False)

    return tempData


def getByDay(data, day):
    tempData = data.query('Day == ' + str(day), inplace=False)

    return tempData


########################################################################################################################


##################################################### Utility ##########################################################
def getDayValue():
    return 'Day'


def getBodyWeightValue():
    return 'BodyWeight'


def getBloodGlucoseValue():
    return 'BloodGlucose'


def getPercentageValue():
    return 'Percentage'

########################################################################################################################
