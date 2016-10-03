import numpy as np
import matplotlib.pyplot as plt
import sys
from PyQt4 import QtGui, QtCore
import csv
from safeEval import safeEval

SPINNAMES = ['0', '1/2', '1', '3/2', '2', '5/2', '3', '7/2', '4', '9/2', '5', '11/2', '6', '13/2', '7', '15/2', '8']
SPINCOLORS = ['white', 'blue', 'orange', 'green', 'yellow', 'red', 'lime', 'olive', 'lightBlue', 'pink']
GAMMASCALE = 100/42.576
with open("IsotopeProperties") as isoFile:
    isoList = [line.strip().split('\t') for line in isoFile]
isoList = isoList[1:]
N = len(isoList)
nameList = []
atomNumList = np.zeros(N)
atomMassList = np.zeros(N)
spinList = np.zeros(N)
abundanceList = np.zeros(N)
magMomentList = np.zeros(N)
gammaList = np.zeros(N)
qList = np.zeros(N)
freqRatioList = np.zeros(N)
refSampleList = []
sampleConditionList = []
linewidthFactorList = np.zeros(N)
dpList = np.zeros(N)
dcList = np.zeros(N)

for i in range(N):
    isoN = isoList[i]
    atomNumList[i] = int(isoN[0])
    nameList = np.append(nameList, isoN[1])
    atomMassList[i] = int(isoN[2])
    if isoN[3] == '-':
        spinList[i] = np.nan
    else:
        spinList[i] = isoN[3]
    if isoN[4] == '-':
        abundanceList[i] = np.nan
    else:
        abundanceList[i] = isoN[4]
    if isoN[5] == '-':
        magMomentList[i] = np.nan
    else:
        magMomentList[i] = isoN[5]
    if isoN[6] == '-':
        gammaList[i] = np.nan
    else:
        gammaList[i] = isoN[6]
    if isoN[7] == '-':
        qList[i] = np.nan
    else:
        qList[i] = isoN[7]
    if isoN[8] == '-':
        freqRatioList[i] = np.nan
    else:
        freqRatioList[i] = isoN[8]
    refSampleList = np.append(refSampleList, isoN[9])
    sampleConditionList = np.append(sampleConditionList, isoN[10])
    if isoN[11] == '-':
        linewidthFactorList[i] = np.nan
    else:
        linewidthFactorList[i] = isoN[11]
    if isoN[12] == '-':
        dpList[i] = np.nan
    else:
        dpList[i] = isoN[12]
    if isoN[13] == '-':
        dcList[i] = np.nan
    else:
        dcList[i] = isoN[13]

# Create a list of structures containing the isotope information
ATOMNUM = 116
MASTERISOTOPELIST = []
LONGEST = 0
for i in range(ATOMNUM):
    select = atomNumList == (i+1)
    LONGEST = np.max((LONGEST, np.sum(select)))
    isotopeEntries = {'name': nameList[select],
                      'mass': atomMassList[select],
                      'spin': spinList[select],
                      'abundance': abundanceList[select],
                      'magMoment': magMomentList[select],
                      'gamma': gammaList[select],
                      'q': qList[select],
                      'freqRatio': freqRatioList[select],
                      'refSample': refSampleList[select],
                      'sampleCondition': sampleConditionList[select],
                      'linewidthFactor': linewidthFactorList[select],
                      'dp': dpList[select],
                      'dc': dcList[select]}
    if len(nameList[select]) > 0:
        MASTERISOTOPELIST.append(isotopeEntries)
    else:
        MASTERISOTOPELIST.append(None)
nameList = sorted(set(nameList))

class PeriodicTable(QtGui.QWidget):
    
    def __init__(self):
        super(PeriodicTable, self).__init__()
        self.freqConst = 6
        self.windowList = []
        self.resetIso()
        self.initUI()
        self.upd()
        
    def resetIso(self):
        self.isoSelect = np.zeros(ATOMNUM)
        for i in range(ATOMNUM):
            if MASTERISOTOPELIST[i] is not None:
                self.isoSelect[i] = np.nanargmax(MASTERISOTOPELIST[i]['abundance'])
            else:
                self.isoSelect[i] = None
        
    def initUI(self):        
        grid = QtGui.QGridLayout()
        self.setLayout(grid)
        count1 = 0
        count2 = 0
        groupList = []
        self.labelList = []
        self.freqEditList = []
        grid.addWidget(PtQLabel('B<sub>0</sub>[T]:'), 0, 2)
        self.b0Entry = PtQLineEdit()
        self.b0Entry.returnPressed.connect(self.setB0)
        grid.addWidget(self.b0Entry, 0, 3)
        grid.addWidget(PtQLabel('Spin:'), 0, 4)
        for i in range(1, len(SPINCOLORS)):
            legendEntry = PtQLineEdit(SPINNAMES[i])
            legendEntry.setStyleSheet('border-style: outset; border-width: 2px; border-color: ' + SPINCOLORS[i] + ';')
            legendEntry.setReadOnly(True)
            grid.addWidget(legendEntry, 0, i+4)
        for i in range(ATOMNUM):
            # groupList.append(QtGui.QGroupBox(str(i+1)))
            groupList.append(QtGui.QWidget())
            groupList[-1].mouseDoubleClickEvent=lambda arg, i=i: self.openWindow(arg, i)
            grid.addWidget(groupList[-1], count1+1, count2)
            count2 += 1
            if count1 is 0 and count2 is 1:
                count2 = 17
            elif (count1 is 1 or count1 is 2) and count2 is 2:
                count2 = 12
            elif (count1 is 5 or count1 is 6) and count2 is 2:
                count1 += 2
            elif (count1 is 7 or count1 is 8) and count2 is 16:
                count2 = 2
                count1 -= 2
            elif count2 is 18:
                count2 = 0
                count1 += 1
            grid2 = QtGui.QGridLayout()
            grid2.setSpacing(0)
            grid2.setMargin(0)
            groupList[i].setLayout(grid2)
            groupList[i].setStyleSheet('background-color: white;')
            self.labelList.append(PtQLabel())
            grid2.addWidget(self.labelList[-1], 0, 0)
            self.freqEditList.append(PtQLineEdit())
            self.freqEditList[-1].returnPressed.connect(lambda i=i: self.setFreq(i))
            grid2.addWidget(self.freqEditList[-1], 1, 0)
        self.setWindowTitle('NMR table')
        self.show()

    def upd(self):
        self.updWindows()
        self.b0Entry.setText('%0.2f' %(self.freqConst*GAMMASCALE))
        for i in range(ATOMNUM):
            if MASTERISOTOPELIST[i] is not None:
                self.labelList[i].setText(str(i+1) + ': <sup>' + str(int(MASTERISOTOPELIST[i]['mass'][int(self.isoSelect[i])])) + '</sup>' + MASTERISOTOPELIST[i]['name'][int(self.isoSelect[i])])
                self.freqEditList[i].setText('%0.2f' %(self.freqConst*MASTERISOTOPELIST[i]['freqRatio'][int(self.isoSelect[i])]))
                self.freqEditList[i].setStyleSheet('border-style: outset; border-width: 2px; border-color: ' + SPINCOLORS[int(2*MASTERISOTOPELIST[i]['spin'][int(self.isoSelect[i])])] + ';')
            else:
                self.labelList[i].setText(str(i+1) + ': ')
                self.freqEditList[i].setText('')

    def updWindows(self):
        for win in self.windowList:
            win.upd()
                
    def openWindow(self, event, n):
        self.windowList.append(DetailWindow(self, n))

    def removeWindow(self, win):
        self.windowList.remove(win)

    def setFreq(self, n):
        if MASTERISOTOPELIST[n] is not None:
            val = safeEval(self.freqEditList[n].text())
            if val is not None:
                self.freqConst = val / MASTERISOTOPELIST[n]['freqRatio'][int(self.isoSelect[n])] 
                self.upd()

    def setB0(self):
        val = safeEval(self.b0Entry.text())
        if val is not None:
            self.freqConst = val / GAMMASCALE
            self.upd()


class DetailWindow(QtGui.QWidget):
    
    def __init__(self, parent, n=0):
        super(DetailWindow, self).__init__()
        self.father = parent
        self.n = n
        self.initUI()
        self.atomSelect(n+1)
        self.show()

    def initUI(self):
        self.setWindowTitle('Details')
        self.grid = QtGui.QGridLayout()
        self.setLayout(self.grid)
        self.grid.addWidget(PtQLabel('N:'), 0, 0)
        self.nSpinBox = QtGui.QSpinBox()
        self.nSpinBox.setMinimum(1)
        self.nSpinBox.setMaximum(ATOMNUM)
        self.nSpinBox.setValue(self.n)
        self.nSpinBox.valueChanged.connect(self.atomSelect)
        self.grid.addWidget(self.nSpinBox, 0, 1)
        self.grid.addWidget(PtQLabel('Name:'), 0, 2)
        self.nameLabel = QtGui.QComboBox()
        self.nameLabel.addItems(list(nameList))
        self.nameLabel.currentIndexChanged[str].connect(self.atomSelectName)
        self.grid.addWidget(self.nameLabel, 0, 3)
        self.grid.addWidget(PtQLabel('Mass:'), 1, 1)
        self.grid.addWidget(PtQLabel('Spin:'), 1, 2)
        self.grid.addWidget(PtQLabel('Abundance:'), 1, 3)
        self.grid.addWidget(PtQLabel('Gamma:'), 1, 4)
        self.grid.addWidget(PtQLabel('Q:'), 1, 5)
        self.grid.addWidget(PtQLabel('Freq:'), 1, 6)
        self.grid.addWidget(PtQLabel('Sample:'), 1, 7)
        self.grid.addWidget(PtQLabel('Condition:'), 1, 8)
        self.grid.addWidget(PtQLabel('Linewidth:'), 1, 9)
        self.grid.addWidget(PtQLabel('dp:'), 1, 10)
        self.grid.addWidget(PtQLabel('dc:'), 1, 11)
        self.buttongroup = QtGui.QButtonGroup()
        self.buttongroup.buttonClicked.connect(self.changeSelect)
        self.radiobuttons = []
        self.massLabels = []
        self.spinLabels = []
        self.abundanceLabels = []
        self.gammaLabels = []
        self.qLabels = []
        self.freqEntries = []
        self.sampleLabels = []
        self.conditionLabels = []
        self.linewidthLabels = []
        self.dpLabels = []
        self.dcLabels = []
        for i in range(LONGEST):
            self.radiobuttons.append(QtGui.QRadioButton())
            self.buttongroup.addButton(self.radiobuttons[-1], i)
            self.grid.addWidget(self.radiobuttons[-1], i+2, 0)
            self.massLabels.append(PtQLabel())
            self.massLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.massLabels[-1], i+2, 1)
            self.spinLabels.append(PtQLabel())
            self.spinLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.spinLabels[-1], i+2, 2)
            self.abundanceLabels.append(PtQLabel())
            self.abundanceLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.abundanceLabels[-1], i+2, 3)
            self.gammaLabels.append(PtQLabel())
            self.gammaLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.gammaLabels[-1], i+2, 4)
            self.qLabels.append(PtQLabel())
            self.qLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.qLabels[-1], i+2, 5)
            self.freqEntries.append(PtQLineEdit())
            self.freqEntries[-1].setFixedWidth(self.freqEntries[-1].sizeHint().width())
            self.freqEntries[-1].returnPressed.connect(lambda i=i: self.setFreq(i))
            self.grid.addWidget(self.freqEntries[-1], i+2, 6)
            self.sampleLabels.append(PtQLabel())
            self.sampleLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.sampleLabels[-1], i+2, 7)
            self.conditionLabels.append(PtQLabel())
            self.conditionLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.conditionLabels[-1], i+2, 8)
            self.linewidthLabels.append(PtQLabel())
            self.linewidthLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.linewidthLabels[-1], i+2, 9)
            self.dpLabels.append(PtQLabel())
            self.dpLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.dpLabels[-1], i+2, 10)
            self.dcLabels.append(PtQLabel())
            self.dcLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.dcLabels[-1], i+2, 11)
        self.grid.setRowStretch(LONGEST+2, 1)
        self.grid.setColumnStretch(12, 1)

    def upd(self):
        self.atomSelect(self.n + 1)
        
    def atomSelectName(self, name):
        for i in range(len(MASTERISOTOPELIST)):
            var = MASTERISOTOPELIST[i]
            if var is not None:
                if name == var['name'][0]:
                    self.nSpinBox.setValue(i+1)
                    return
        
    def atomSelect(self, n):
        self.n = n-1
        atomProp = MASTERISOTOPELIST[self.n]
        if atomProp is None:
            self.display(0)
            return
        num = len(atomProp['mass'])
        index = self.nameLabel.findText(atomProp['name'][0])
        self.nameLabel.setCurrentIndex(index)
        self.radiobuttons[int(self.father.isoSelect[self.n])].setChecked(True)
        for i in range(num):
            self.massLabels[i].setText(str(int(atomProp['mass'][i])))
            self.spinLabels[i].setText(SPINNAMES[int(2*atomProp['spin'][i])])
            self.abundanceLabels[i].setText(str(atomProp['abundance'][i]))
            self.gammaLabels[i].setText(str(atomProp['gamma'][i]))
            self.qLabels[i].setText(str(atomProp['q'][i]))
            self.freqEntries[i].setText(str(self.father.freqConst*atomProp['freqRatio'][i]))
            self.sampleLabels[i].setText(str(atomProp['refSample'][i]))
            self.conditionLabels[i].setText(str(atomProp['sampleCondition'][i]))
            self.linewidthLabels[i].setText(str(atomProp['linewidthFactor'][i]))
            self.dpLabels[i].setText(str(atomProp['dp'][i]))
            self.dcLabels[i].setText(str(atomProp['dc'][i]))
        self.display(num)

    def changeSelect(self, n):
        self.father.isoSelect[self.n] = self.buttongroup.checkedId()
        self.father.upd()

    def setFreq(self, i):
        val = safeEval(self.freqEntries[i].text())
        if val is not None:
            self.father.freqConst = val / MASTERISOTOPELIST[self.n]['freqRatio'][i]
            self.father.upd()
            self.atomSelect(self.n+1)
        
    def display(self, num):
        for i in range(num):
            self.radiobuttons[i].show()
            self.massLabels[i].show()
            self.spinLabels[i].show()
            self.abundanceLabels[i].show()
            self.gammaLabels[i].show()
            self.qLabels[i].show()
            self.freqEntries[i].show()
            self.sampleLabels[i].show()
            self.conditionLabels[i].show()
            self.linewidthLabels[i].show()
            self.dpLabels[i].show()
            self.dcLabels[i].show()
        for i in range(num, LONGEST):
            self.radiobuttons[i].hide()
            self.massLabels[i].hide()
            self.spinLabels[i].hide()
            self.abundanceLabels[i].hide()
            self.gammaLabels[i].hide()
            self.qLabels[i].hide()
            self.freqEntries[i].hide()
            self.sampleLabels[i].hide()
            self.conditionLabels[i].hide()
            self.linewidthLabels[i].hide()
            self.dpLabels[i].hide()
            self.dcLabels[i].hide()

    def closeEvent(self, event):
        super(DetailWindow, self).closeEvent(event)
        self.father.removeWindow(self)


class PtQLabel(QtGui.QLabel):
    def __init__(self, parent=None):
        QtGui.QLabel.__init__(self, parent)
        self.setAlignment(QtCore.Qt.AlignCenter)


class PtQLineEdit(QtGui.QLineEdit):
    def __init__(self, parent=None):
        QtGui.QLineEdit.__init__(self, parent)
        self.setAlignment(QtCore.Qt.AlignCenter)


def main():
    app = QtGui.QApplication(sys.argv)
    ex = PeriodicTable()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
