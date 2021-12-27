import DiseaseProgressionLib as DLib
import UserInterface as UI


def main():
    NCData, DCData, STD, test1, test2, test3 = DLib.readData()  # get dataframes

    massagedData = DLib.massageData(NCData, DCData)  # massage data

    # training data
    DLib.fitPlane(massagedData, DLib.getBodyWeightValue(), DLib.getBloodGlucoseValue(),
                                          DLib.getPercentageValue())

    while(1):
        action = int(input('Press 1 for displaying trained function \n'
                           'Press 2 for comparing disease progression after inducing different drugs \n'
                           'Press 3 to exit \n'))
        if action == 1:
            DLib.drawPlaneAndScatterPlot(massagedData, DLib.getBodyWeightValue(), DLib.getBloodGlucoseValue(),
                                         DLib.getPercentageValue())
        elif action == 2:
            STDData = DLib.predict(NCData, STD)
            test1Data = DLib.predict(NCData, test1)
            test2Data = DLib.predict(NCData, test2)
            test3Data = DLib.predict(NCData, test3)
            DLib.drawSubPlotToCompare(STDData, test1Data, test2Data, test3Data)
        elif action == 3:
            break


main()