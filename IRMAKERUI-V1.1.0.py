import numpy as np

from scipy.io.wavfile import read, write
from scipy.signal import convolve, spectrogram

from PyQt5.QtWidgets import *
from PyQt5 import QtCore, uic
import pyqtgraph as pg
from pathlib import Path

from Sweep_Window import Sweep_Window
from funcs import next_power_of_2,smooth

import sys

# Main Window
class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setWindowTitle('IR Maker')
        # load ui made with qt designer
        uic.loadUi('design_main_window.ui',self)

        # get elements for the engine
        # sweep section
        self.browsesweep_button = self.findChild(QPushButton, 'browseSweep_button')
        self.beg_freq = self.findChild(QLineEdit, 'begfreq_txtedit')
        self.beg_freq.setPlaceholderText('in Hz')
        self.end_freq = self.findChild(QLineEdit, 'endfreq_txtedit')
        self.end_freq.setPlaceholderText('in Hz')
        self.sweepgen_button = self.findChild(QPushButton, 'sweepgen_button')
        self.sweep_path = ''
        # response section
        self.loadfolder_button = self.findChild(QPushButton,'selectfolder_button')
        self.files_list = self.findChild(QListView,'files_list')
        self.fileModel = QFileSystemModel()
        self.files_list.setModel(self.fileModel)

        # output section
        self.autosave_radio = self.findChild(QRadioButton,'autosave_radio')
        self.customsave_radio = self.findChild(QRadioButton,'customsave_radio')
        self.browseout_button = self.findChild(QPushButton,'browseout_button')
        self.bitdepth_combo = self.findChild(QComboBox,'bitdepth_combo')
        self.mpt_checkbox = self.findChild(QCheckBox,'mpt_checkbox')
        # execute section
        self.createir_button = self.findChild(QPushButton,'compute_button')
        self.playir_button = self.findChild(QPushButton,'play_ir')
        # plots
        self.ir_plot = self.findChild(pg.PlotWidget, "temporal_plot")
        self.ir_plot.setBackground('w')
        self.spectral_plot = self.findChild(pg.PlotWidget,"spectral_plot")
        self.spectral_plot.setBackground('w')

        # engine
        # sweep section
        self.browsesweep_button.clicked.connect(lambda: self.openFileSweep())
        self.sweepgen_button.clicked.connect(lambda: self.sweep())
        # response section
        self.loadfolder_button.clicked.connect(lambda: self.openFileResponse())

        self.show()
    
    
    # functions section
    def openFileSweep(self):
        fileName, _ = QFileDialog.getOpenFileName(self,"Select sweep file", "","*wav")
        test = read(fileName)[1] 
        if test.ndim == 2 :
            _,c = np.shape(test)
            while c > 2 :
                warn=QMessageBox()
                warn.setText('File contains too much channels. Please choose a mono or stereo file only.')
                warn.setTitle('Error')
                fileName, _ = QFileDialog.getOpenFileName(self,"Please choose a mono or stereo file only", "",'*.wav')
                test = read(fileName)[1] 
                if test.ndim == 1:
                    c=1
                else :
                    _,c = np.shape(test)
        self.sweep_path=fileName
        
        self.browsesweep_button.setText(Path(fileName).name)
        return

    def openFileResponse(self):
        self.fileModel.setFilter(QtCore.QDir.NoDotAndDotDot |  QtCore.QDir.Files)
        self.fileModel.setNameFilters(['*.wav'])
        self.fileModel.setNameFilterDisables(False)
        folderpath = QFileDialog.getExistingDirectory(self, 'Select recordings folder')
        self.files_list.setRootIndex(self.fileModel.index(folderpath))
        self.loadfolder_button.setText(Path(folderpath).name)
        return
    
    def sweep(self):
        dialog=Sweep_Window(self)
        dialog.show()
        return
    



# run app
app = QApplication(sys.argv)
main_window = MainWindow()
app.exec_()