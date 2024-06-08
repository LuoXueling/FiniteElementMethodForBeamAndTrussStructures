import sys
from PyQt5.QtWidgets import QWidget, QToolTip, QPushButton, QApplication, QDesktopWidget, QMainWindow, \
    QAction, QPushButton, QGridLayout, QLabel, QLineEdit, QFileDialog, QTextEdit
from PyQt5.QtGui import QFont
from PyQt5 import QtGui
import os
import sys
import numpy as np
import time
import newStructure
import plotStructure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar


# 命令行输出重定向
class EmittingStream():
    def __init__(self, handle):
        self.textWritten = handle

    def write(self, text):
        self.textWritten(str(text))

    def flush(self):
        pass


# 命令行输入重定向
class receiveStream():
    def __init__(self, handle):
        self.textReceive = handle

    def readline(self):
        return self.textReceive()


class MyWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.cmdFlag = False
        self.initUI()

    def initUI(self):
        QToolTip.setFont(QFont('SansSerif', 10))

        # 显示菜单
        self.initMenu()

        widget = QWidget()

        self.grid = QGridLayout()
        self.grid.setSpacing(10)

        # 输入文件显示框
        inputfile = QLabel('Input file')
        self.inputfileEdit = QLineEdit()
        self.grid.addWidget(inputfile, 1, 0)
        self.grid.addWidget(self.inputfileEdit, 1, 1, 1, 5)
        # 新建文件按钮
        newfilebtn = QPushButton('New')
        newfilebtn.clicked.connect(self.newStruc)
        self.grid.addWidget(newfilebtn, 1, 6)
        # 选择文件按钮
        inputfilebtn = QPushButton('Select')
        inputfilebtn.clicked.connect(self.selectfile)
        self.grid.addWidget(inputfilebtn, 1, 7)
        # 作图按钮
        drawbtn = QPushButton('Plot')
        drawbtn.clicked.connect(self.draw)
        self.grid.addWidget(drawbtn, 1, 8)
        # 计算按钮
        runbtn = QPushButton('Compute')
        runbtn.clicked.connect(self.start)
        self.grid.addWidget(runbtn, 1, 9)
        # 保存图片按钮
        runbtn = QPushButton('Save fig')
        runbtn.clicked.connect(self.saveFig)
        self.grid.addWidget(runbtn, 1, 10)

        # 控制台信息输出行
        self.sysout = QTextEdit()
        self.sysout.setReadOnly(True)
        self.sysout.setObjectName('sysout')
        self.grid.addWidget(self.sysout, 2, 6, 10, 7)
        sys.stdout = EmittingStream(self.outputWritten)  # 重定向输出
        # sys.stderr = EmittingStream(self.outputWritten)
        # 命令输入行
        command = QLabel('  Command  >>')
        self.grid.addWidget(command, 12, 6)
        self.cmdEdit = QLineEdit()
        self.cmdEdit.returnPressed.connect(self.cmdProcess)
        self.grid.addWidget(self.cmdEdit, 12, 7, 1, 6)
        sys.stdin = receiveStream(self.inputReceive)  # 重定向输入

        # 将网格加入mainwindow
        widget.setLayout(self.grid)
        self.setCentralWidget(widget)

        # 控制窗口大小与位置
        self.resize(1600, 800)
        self.center()

        # self.setGeometry(300, 300, 300, 200)
        self.setWindowTitle('CSM')
        self.show()

    # 控制窗口显示在屏幕中心的方法
    def center(self):
        # 获得窗口
        qr = self.frameGeometry()
        # 获得屏幕中心点
        cp = QDesktopWidget().availableGeometry().center()
        # 显示到屏幕中心
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    # 初始化菜单栏
    def initMenu(self):
        # 显示菜单栏
        menubar = self.menuBar()
        # 添加菜单
        fileMenu = menubar.addMenu('File')
        aboutMenu = menubar.addMenu('About')
        # 菜单项：文件
        actionfileopen = QAction('&SelectFile', self)
        actionfileopen.setObjectName("actionfileopen")
        actionfileopen.triggered.connect(self.selectfile)
        fileMenu.addAction(actionfileopen)
        # 菜单项：新建文件
        actionnewstruc = QAction('&NewStructure', self)
        actionnewstruc.setObjectName("actionnewstruc")
        actionnewstruc.triggered.connect(self.newStruc)
        fileMenu.addAction(actionnewstruc)


    # 通过菜单选择文件
    def selectfile(self):
        print('\n--------------------------------SelectFile--------------------------------\n')
        fileName, fileType = QFileDialog.getOpenFileName(self, "SelectFile", os.getcwd(),
                                                         "Text Files(*.txt)")
        self.inputfileEdit.setText(fileName)
        print('\n----------------SelectFile ended: ' + self.inputfileEdit.text() + ' ----------------\n')

    def draw(self):
        print('\n----------------Plotting: '+self.inputfileEdit.text()+' ----------------\n')
        self.F = plotStructure.run(self.inputfileEdit.text())
        self.grid.addWidget(self.F, 2, 0, 7, 6)
        mpl_ntb = NavigationToolbar(self.F, self)
        self.grid.addWidget(mpl_ntb, 8, 0, 1, 6)
        print('\n----------------Plotting ended: ' + self.inputfileEdit.text() + ' ----------------\n')

    # 点击start按钮
    def start(self):
        print('\n----------------Computation: ' + self.inputfileEdit.text() + ' ----------------\n')
        if self.inputfileEdit.text() == '':
            print('File not selected')
            return 0
        try:
            root, filename, inpname = plotStructure.path_preprocess(self.inputfileEdit.text())
            # 开始计算
            selfpath = os.path.realpath(sys.argv[0]).replace("ui.py", "")
            exe_path = selfpath + "csm.exe"
            params = "\"" + self.inputfileEdit.text() + " " + root + inpname + " " + inpname + "\\\""
            params = params.replace('/', '\\')
            with open('./runCSM.bat', 'w') as file:
                # 不生成新
                file.write("@echo off\n")
                file.write("powershell.exe -command ^\n")
                file.write("  \"start " + exe_path + " \\" + params + "\" -NoNewWindow -Wait")
            os.system('runCSM.bat')
        except Exception as e:
            print('Unable to compute:')
            print(e)
        # 绘图
        self.F = plotStructure.run(self.inputfileEdit.text())
        self.grid.addWidget(self.F, 2, 0, 7, 6)
        mpl_ntb = NavigationToolbar(self.F, self)
        self.grid.addWidget(mpl_ntb, 8, 0, 1, 6)
        print('\n----------------Computation ended: ' + self.inputfileEdit.text() + ' ----------------\n')

    # 点击新建按钮
    def newStruc(self):
        print('\n--------------------------------New structure--------------------------------\n')
        newStructure.run()
        fileName, fileType = QFileDialog.getSaveFileName(self, "SaveFile", os.getcwd(),
                                                         "Text Files(*.txt)")
        newStructure.save(fileName)
        self.inputfileEdit.setText(fileName)
        print('\n----------------SaveFile ended: ' + self.inputfileEdit.text() + ' ----------------\n')

    def saveFig(self):
        print('\n--------------------------------SaveFig--------------------------------\n')
        fileName, fileType = QFileDialog.getSaveFileName(self, "SaveFig", os.getcwd(),
                                                         "Image Files(*.png)")
        self.F.fig.savefig(fileName, dpi=300)
        print('\n----------------SaveFig ended: ' + self.inputfileEdit.text() + ' ----------------\n')

    # 将text在textedit中输出
    def outputWritten(self, text):
        if str(text) == 'Previewing structure':
            global info, node, element, restrict, force
            F = plotStructure.plot_structure(newStructure.info, newStructure.node, element=newStructure.element,
                                             restrict=newStructure.restrict, force=newStructure.force)
            self.grid.addWidget(F, 2, 0, 7, 6)
            mpl_ntb = NavigationToolbar(F, self)
            self.grid.addWidget(mpl_ntb, 8, 0, 1, 6)
            return 0
        cursor = self.sysout.textCursor()
        cursor.movePosition(QtGui.QTextCursor.End)
        cursor.insertText(text)
        self.sysout.setTextCursor(cursor)
        self.sysout.ensureCursorVisible()

    # 接收lineedit中的输入
    def inputReceive(self):
        # 接收到回车是开始输入
        while self.cmdFlag == False:
            QApplication.processEvents()
            time.sleep(0.01)
        s = self.cmdEdit.text()
        self.cmdEdit.clear()
        print('>> ' + s)
        self.cmdFlag = False
        return s

    # 接收到回车时触发的函数，让input停止等待
    def cmdProcess(self):
        self.cmdFlag = True


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MyWindow()
    sys.exit(app.exec_())
