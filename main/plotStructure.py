# -*- coding:utf-8 -*-
# 读取输入文件和结果，绘制节点、单元、力、约束、变形、应力图

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from matplotlib import patches as mpatches
import tkinter as tk
from tkinter import filedialog
import os
import numpy as np
import warnings
import copy
from matplotlib.widgets import Cursor

warnings.filterwarnings("ignore")

plt.rcParams['figure.dpi'] = 150

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

tmppressed = -1


class MyFigure(FigureCanvas):
    def __init__(self, figsize=(5, 6)):
        # 第一步：创建一个创建Figure
        self.fig = Figure(figsize=figsize)
        # 第二步：在父类中激活Figure窗口
        super(MyFigure, self).__init__(self.fig)  # 此句必不可少，否则不能显示图形


def run(filepath, oncanvas=True):
    if oncanvas:
        mpl.use('Qt5Agg')
    info, node, element, restrict, force = read_inputfile(filepath)
    root, filename, inpname = path_preprocess(filepath)
    print('\n--------------------------------Results--------------------------------\n')
    try:
        displacement = read_displacement(info, root, inpname)
        print("Nodal displacement: ")
        for i in displacement:
            print('\t', i)
    except:
        displacement = None
    try:
        restrictForce = read_restrictForce(info, root, inpname)
        print("Constraint force: ")
        for i in restrictForce:
            print('\t', i)
    except:
        restrictForce = None
    try:
        if info['mode'] == 'bar':
            stress = read_stress(info, root, inpname)
            print("Stress: ")
            for i in stress:
                print('\t', i)
        else:
            stress = None
    except:
        stress = None

    F = plot_structure(info, node, element, restrict, force, displacement, stress, restrictForce, oncanvas=oncanvas)
    return F


def read_inputfile(inp_filepath):
    info = {}
    node = []
    element = []
    restrict = []
    force = []

    with open(inp_filepath, 'r')as file:
        line = file.readline().replace("\n", "")
        info['mode'] = line
        line = rdline(file.readline())
        info['nEleNode'] = int(line[0])
        info['nNodeFdm'] = int(line[1])
        info['nNode'] = int(line[2])
        info['nEle'] = int(line[3])
        info['nType'] = int(line[4])
        info['nRestrict'] = int(line[5])
        line = rdline(file.readline())
        info['eleType'] = list(map(int, line))
        info['eleMaterial'] = []
        for n in range(info['nType']):
            line = rdline(file.readline())
            info['eleMaterial'].append(list(map(float, line)))
        for n in range(info['nNode']):
            line = rdline(file.readline())
            node.append(list(map(float, line)))
        if info['mode'] == 'beam':  # 跳过辅助节点
            for n in range(info['nEle']):
                line = rdline(file.readline())
                node.append(list(map(float, line)))
        for n in range(info['nEle']):
            line = rdline(file.readline())
            element.append(list(map(int, line)))
        line = rdline(file.readline())
        restrict = list(map(int, line))
        line = rdline(file.readline())
        tmp = []
        for i in range(info['nNode'] * info['nNodeFdm']):
            tmp.append(float(line[i]))
            if (i + 1) % info['nNodeFdm'] == 0:
                force.append(tmp)
                tmp = []
        print('\n--------------------------------Structure information--------------------------------\n')
        print("Structure information: ", info)
        print("Coordinates of nodes: ", node)
        print("Nodes of elements: ", element)
        print("Constrained coordinates: ", restrict)
        print("Loads of nodes", force)
        print('\n-----------------------------------------------------------------------\n')
        file.close()
    return info, node, element, restrict, force


def read_displacement(info, root, inpname):
    disp_filepath = root + inpname + '//' + inpname + '_nodeDisplacement.txt'
    file = open(disp_filepath, 'r')
    displacement = []
    tmp = []
    for i in range(info['nNode'] * info['nNodeFdm']):
        line = rdline(file.readline())
        tmp.append(float(line[1]))
        if (i + 1) % info['nNodeFdm'] == 0:
            displacement.append(tmp)
            tmp = []
    return displacement


def read_restrictForce(info, root, inpname):
    resf_filepath = root + inpname + '//' + inpname + '_restrictForce.txt'
    file = open(resf_filepath, 'r')
    restrictForce = []
    tmp = []
    for i in range(info['nEle'] * (info['nNodeFdm'] * info['nEleNode'] + 1)):
        file.readline()
    for i in range(info['nNode'] * info['nNodeFdm']):
        line = rdline(file.readline())
        tmp.append(float(line[-1]))
        if (i + 1) % info['nNodeFdm'] == 0:
            restrictForce.append(tmp)
            tmp = []
    return restrictForce


def read_stress(info, root, inpname):
    disp_filepath = root + inpname + '//' + inpname + '_stress.txt'
    file = open(disp_filepath, 'r')
    stress = []
    for i in range(info['nEle']):
        line = rdline(file.readline())
        stress.append(float(line[0]))
    return stress


def plot_structure(info, node=None, element=None, restrict=None, force=None, displacement=None, stress=None,
                   restrictForce=None, oncanvas=True):
    plt.ion()
    if oncanvas == True:
        F = MyFigure()
        if stress != None:
            ax = F.fig.add_axes([0.15, 0.05, 0.8, 0.9], projection='3d')
        else:
            # fig = plt.figure(figsize=(7, 6))
            ax = F.fig.gca(projection='3d')
    else:
        F = plt.figure(figsize=(7, 6))
        if stress != None:
            ax = F.add_axes([0.15, 0.05, 0.8, 0.9], projection='3d')
        else:
            # fig = plt.figure(figsize=(7, 6))
            ax = F.add_subplot(projection='3d')
    #######################################################
    # 记录最大值，用于控制缩放
    X1 = Y1 = Z1 = X2 = Y2 = Z2 = []
    if node != None:
        X1 = np.array([node[n][0] for n in range(info['nNode'])])
        Y1 = np.array([node[n][1] for n in range(info['nNode'])])
        Z1 = np.array([node[n][2] for n in range(info['nNode'])])

    if displacement != None:
        X2 = np.array([node[n][0] + displacement[n][0] for n in range(info['nNode'])])
        Y2 = np.array([node[n][1] + displacement[n][1] for n in range(info['nNode'])])
        Z2 = np.array([node[n][2] + displacement[n][2] for n in range(info['nNode'])])
        X = np.hstack((X1, X2))
        Y = np.hstack((Y1, Y2))
        Z = np.hstack((Z1, Z2))
    else:
        X = X1
        Y = Y1
        Z = Z1
    max_range = np.max(np.array([np.max(X) - np.min(X), np.max(Y) - np.min(Y), np.max(Z) - np.min(Z)]))
    max_origin = np.max(np.array([np.max(X1) - np.min(X1), np.max(Y1) - np.min(Y1), np.max(Z1) - np.min(Z1)]))

    # 计算缩放系数
    if displacement != None:
        scale_factor = max_range / np.max(np.max(np.abs(displacement))) * 0.1
        print("Scale factor: ", scale_factor)
    else:
        scale_factor = 1

    #######################################################
    # 绘制变形后节点
    if displacement != None:
        X2 = np.array([node[n][0] + displacement[n][0] * scale_factor for n in range(info['nNode'])])
        Y2 = np.array([node[n][1] + displacement[n][1] * scale_factor for n in range(info['nNode'])])
        Z2 = np.array([node[n][2] + displacement[n][2] * scale_factor for n in range(info['nNode'])])
        if info['mode'] == 'bar':
            ax.scatter(X2, Y2, Z2, color='grey')

    global nodeplot
    # 绘制变形前节点
    if node != None:
        X1 = np.array([node[n][0] for n in range(info['nNode'])])
        Y1 = np.array([node[n][1] for n in range(info['nNode'])])
        Z1 = np.array([node[n][2] for n in range(info['nNode'])])
        if info['mode'] == 'bar':
            nodeplot = ax.scatter(X1, Y1, Z1, color='blue', zorder=0)

    # 绘制变形后单元
    if displacement != None and info['mode'] == 'bar':
        for i in range(info['nEle']):
            x = [X2[n - 1] for n in element[i][:info['nEleNode']]]
            y = [Y2[n - 1] for n in element[i][:info['nEleNode']]]
            z = [Z2[n - 1] for n in element[i][:info['nEleNode']]]
            ax.plot(x, y, z, color='grey')
    # 对于梁单元，要使用形函数作曲线
    if displacement != None and info['mode'] == 'beam':
        for i in range(info['nEle']):
            x1, x2 = [node[n - 1][0] for n in element[i][:info['nEleNode']]]
            y1, y2 = [node[n - 1][1] for n in element[i][:info['nEleNode']]]
            z1, z2 = [node[n - 1][2] for n in element[i][:info['nEleNode']]]
            u1, u2 = [displacement[n - 1][0] for n in element[i][:info['nEleNode']]]
            v1, v2 = [displacement[n - 1][1] for n in element[i][:info['nEleNode']]]
            w1, w2 = [displacement[n - 1][2] for n in element[i][:info['nEleNode']]]
            tx1, tx2 = [displacement[n - 1][3] for n in element[i][:info['nEleNode']]]
            ty1, ty2 = [displacement[n - 1][4] for n in element[i][:info['nEleNode']]]
            tz1, tz2 = [displacement[n - 1][5] for n in element[i][:info['nEleNode']]]

            linear = np.linspace(0, 1, 50)
            u = np.zeros((50))
            v = np.zeros((50))
            w = np.zeros((50))

            # A在不同方向的弯曲上第2、4列符号不同，但这两类的符号变化是一致的
            # 为了方便处理，这里让d矩阵中转角的符号改变
            # 初步的规律是：xy、yz、zx方向的弯曲为正，yx、zy、xz的弯曲为负
            # 但这一规律在x2<x1、y2<y1、z2<z1时相反，应该是因为p68最下方边界条件x不变，但v、t交换了

            l = abs(x2 - x1)
            xs = np.linspace(0, l, 50)
            flag = 1
            if x2 < x1:
                flag = -1
            if abs(l) > 1e-5:
                A = np.array([[1, 0, 0, 0,
                               0, 1, 0, 0,
                               -3 / l ** 2, -2 / l, 3 / l ** 2, -1 / l,
                               2 / l ** 3, 1 / l ** 2, -2 / l ** 3, 1 / l ** 2]]).reshape(4, 4)
                for i, x in enumerate(xs):
                    # xy
                    tmp = np.dot(np.array([1, x, x ** 2, x ** 3]).reshape(1, 4), A)
                    v[i] += np.dot(tmp, np.array([v1, flag * tz1, v2, flag * tz2]).reshape(4, 1))
                    # xz
                    tmp = np.dot(np.array([1, x, x ** 2, x ** 3]).reshape(1, 4), A)
                    w[i] += np.dot(tmp, np.array([w1, -flag * ty1, w2, -flag * ty2]).reshape(4, 1))
            # else:
            #     u += u1 * (1 - linear) + u2 * linear

            l = abs(y2 - y1)
            ys = np.linspace(0, l, 50)
            flag = 1
            if y2 < y1:
                flag = -1
            if abs(l) > 1e-5:
                A = np.array([[1, 0, 0, 0,
                               0, 1, 0, 0,
                               -3 / l ** 2, -2 / l, 3 / l ** 2, -1 / l,
                               2 / l ** 3, 1 / l ** 2, -2 / l ** 3, 1 / l ** 2]]).reshape(4, 4)
                for i, y in enumerate(ys):
                    # yz
                    tmp = np.dot(np.array([1, y, y ** 2, y ** 3]).reshape(1, 4), A)
                    w[i] += np.dot(tmp, np.array([w1, flag * tx1, w2, flag * tx2]).reshape(4, 1))
                    # yx
                    tmp = np.dot(np.array([1, y, y ** 2, y ** 3]).reshape(1, 4), A)
                    u[i] += np.dot(tmp, np.array([u1, -flag * tz1, u2, -flag * tz2]).reshape(4, 1))  ##tz2正确
            # else:
            #     v += v1 * (1 - linear) + v2 * linear

            l = abs(z2 - z1)
            zs = np.linspace(0, l, 50)
            flag = 1
            if z2 < z1:
                flag = -1
            if abs(l) > 1e-5:
                A = np.array([[1, 0, 0, 0,
                               0, 1, 0, 0,
                               -3 / l ** 2, -2 / l, 3 / l ** 2, -1 / l,
                               2 / l ** 3, 1 / l ** 2, -2 / l ** 3, 1 / l ** 2]]).reshape(4, 4)
                for i, z in enumerate(zs):
                    # zx
                    tmp = np.dot(np.array([1, z, z ** 2, z ** 3]).reshape(1, 4), A)
                    u[i] += np.dot(tmp, np.array([u1, flag * ty1, u2, flag * ty2]).reshape(4, 1))
                    # zy
                    tmp = np.dot(np.array([1, z, z ** 2, z ** 3]).reshape(1, 4), A)
                    v[i] += np.dot(tmp, np.array([v1, -flag * tx1, v2, -flag * tx2]).reshape(4, 1))
            # else:
            #     w += w1 * (1 - linear) + w2 * linear

            # 如果在某个方向上没有弯曲，则要加上刚体位移（线性项）
            # 因为在计算弯曲时会计算刚体位移
            if np.linalg.norm(u) < 1e-40:
                u += u1 * (1 - linear) + u2 * linear
            if np.linalg.norm(v) < 1e-40:
                v += v1 * (1 - linear) + v2 * linear
            if np.linalg.norm(w) < 1e-40:
                w += w1 * (1 - linear) + w2 * linear

            xs = np.linspace(x1, x2, 50)
            ys = np.linspace(y1, y2, 50)
            zs = np.linspace(z1, z2, 50)
            ax.plot(xs + u * scale_factor, ys + v * scale_factor, zs + w * scale_factor, color='grey')

    # 绘制变形前单元
    grad = np.linspace(0.0, 1.0, 100)
    rgb = plt.get_cmap('cool')(grad)[np.newaxis, :, :3][0]
    if element != None:
        for i in range(info['nEle']):
            if stress != None:
                n = int(np.floor((stress[i] - np.min(stress)) / (np.max(stress) - np.min(stress)) * (len(rgb) - 1)))
                clr = rgb[n]
            else:
                clr = 'black'
            x = [X1[n - 1] for n in element[i][:info['nEleNode']]]
            y = [Y1[n - 1] for n in element[i][:info['nEleNode']]]
            z = [Z1[n - 1] for n in element[i][:info['nEleNode']]]
            ax.plot(x, y, z, c=clr)

    #######################################################
    # 绘制力
    # https://matplotlib.org/2.0.2/mpl_toolkits/mplot3d/tutorial.html
    if force != None and node != None:
        for i, fc in enumerate(force):
            l = np.sqrt(fc[0] ** 2 + fc[1] ** 2 + fc[2] ** 2)
            ax.quiver(node[i][0], node[i][1], node[i][2], 0.3 * fc[0] / l, 0.3 * fc[1] / l, 0.3 * fc[2] / l,
                      color='red')
            if info['mode'] == 'beam':
                if abs(fc[3]) > 1e-10:
                    xs = node[i][0] * np.zeros(50) + node[i][0]
                    ys = (max_origin * 0.05) * np.cos(np.linspace(0, 2 * np.pi, 50)) + node[i][1]
                    zs = (max_origin * 0.05) * np.sin(np.linspace(0, 2 * np.pi, 50)) + node[i][2]
                    ax.plot(xs, ys, zs, color='red')
                if abs(fc[4]) > 1e-10:
                    xs = (max_origin * 0.05) * np.cos(np.linspace(0, 2 * np.pi, 50)) + node[i][0]
                    ys = node[i][1] * np.zeros(50) + node[i][1]
                    zs = (max_origin * 0.05) * np.sin(np.linspace(0, 2 * np.pi, 50)) + node[i][2]
                    ax.plot(xs, ys, zs, color='red')
                if abs(fc[5]) > 1e-10:
                    xs = (max_origin * 0.05) * np.cos(np.linspace(0, 2 * np.pi, 50)) + node[i][0]
                    ys = (max_origin * 0.05) * np.sin(np.linspace(0, 2 * np.pi, 50)) + node[i][1]
                    zs = node[i][2] * np.zeros(50) + node[i][2]
                    ax.plot(xs, ys, zs, color='red')

    # 绘制约束，用线段代替
    if restrict != None and node != None:
        for i, rs in enumerate(restrict):
            xs = ys = zs = []
            rs_node = (rs - 1) // info['nNodeFdm']
            rs_coord = rs - info['nNodeFdm'] * rs_node - 1
            if rs_coord == 0:
                # x方向约束
                xs = [node[rs_node][0], node[rs_node][0] - max_origin * 0.05]
                ys = [node[rs_node][1], node[rs_node][1]]
                zs = [node[rs_node][2], node[rs_node][2]]
            elif rs_coord == 1:
                # y方向约束
                xs = [node[rs_node][0], node[rs_node][0]]
                ys = [node[rs_node][1], node[rs_node][1] - max_origin * 0.05]
                zs = [node[rs_node][2], node[rs_node][2]]
            elif rs_coord == 2:
                # z方向约束
                xs = [node[rs_node][0], node[rs_node][0]]
                ys = [node[rs_node][1], node[rs_node][1]]
                zs = [node[rs_node][2], node[rs_node][2] - max_origin * 0.05]
            elif rs_coord == 3:  # 约束x旋转，x=node[x]
                xs = node[rs_node][0] * np.zeros(50) + node[rs_node][0]
                ys = (max_origin * 0.05) * np.cos(np.linspace(0, 2 * np.pi, 50)) + node[rs_node][1]
                zs = (max_origin * 0.05) * np.sin(np.linspace(0, 2 * np.pi, 50)) + node[rs_node][2]
            elif rs_coord == 4:  # 约束y旋转
                xs = (max_origin * 0.05) * np.cos(np.linspace(0, 2 * np.pi, 50)) + node[rs_node][0]
                ys = node[rs_node][1] * np.zeros(50) + node[rs_node][1]
                zs = (max_origin * 0.05) * np.sin(np.linspace(0, 2 * np.pi, 50)) + node[rs_node][2]
            elif rs_coord == 5:
                xs = (max_origin * 0.05) * np.cos(np.linspace(0, 2 * np.pi, 50)) + node[rs_node][0]
                ys = (max_origin * 0.05) * np.sin(np.linspace(0, 2 * np.pi, 50)) + node[rs_node][1]
                zs = node[rs_node][2] * np.zeros(50) + node[rs_node][2]
            ax.plot(xs, ys, zs, color='blue')

    #######################################################
    # 控制坐标轴比例
    # https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
    # Create cubic bounding box to simulate equal aspect ratio
    Xb = 0.5 * max_origin * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() * 4 / 3 + 0.5 * (np.max(X1) + np.min(X1))
    Yb = 0.5 * max_origin * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() * 4 / 3 + 0.5 * (np.max(Y1) + np.min(Y1))
    Zb = 0.5 * max_origin * np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() + 0.5 * (np.max(Z1) + np.min(Z1))
    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    st = 'Auto-scale factor: ' + str(round(scale_factor, 2))
    ax.set_title(st, pad=20)

    #######################################################
    # 绘制colorbar
    if stress != None:
        if oncanvas == True:
            ax2 = F.fig.add_axes([0.15, 0.15, 0.03, 0.7])
        else:
            ax2 = F.add_axes([0.15, 0.15, 0.03, 0.7])
        x = [0] * 100
        c = np.linspace(0.0, 1.0, 100)
        y = np.linspace(np.min(stress), np.max(stress), 100)
        ax2.scatter(x, y, c=c, s=300, linewidths=0.0, cmap='cool')
        ax2.set_xticks([])
        ax2.set_ylim([np.min(stress), np.max(stress)])
        ax2.set_ylabel('stress(MPa)')
    #######################################################
    # 重载mpl_toolkits\mplot3d\axis3d.py中的format_coord方法
    # 从而随时更新currentpos用于处理点击事件
    global currentpos
    currentpos = []

    def format_coord(xd, yd):
        global currentpos
        import mpl_toolkits.mplot3d.proj3d as proj3d
        """
        Given the 2D view coordinates attempt to guess a 3D coordinate.
        Looks for the nearest edge to the point and then assumes that
        the point is at the same z location as the nearest point on the edge.
        """

        if ax.M is None:
            return ''

        if ax.button_pressed in ax._rotate_btn:
            return 'azimuth={:.0f} deg, elevation={:.0f} deg '.format(
                ax.azim, ax.elev)
            # ignore xd and yd and display angles instead

        # nearest edge
        p0, p1 = min(ax.tunit_edges(),
                     key=lambda edge: proj3d._line2d_seg_dist(
                         edge[0], edge[1], (xd, yd)))

        # scale the z value to match
        x0, y0, z0 = p0
        x1, y1, z1 = p1
        d0 = np.hypot(x0 - xd, y0 - yd)
        d1 = np.hypot(x1 - xd, y1 - yd)
        dt = d0 + d1
        z = d1 / dt * z0 + d0 / dt * z1

        x, y, z = proj3d.inv_transform(xd, yd, z, ax.M)

        xs = ax.format_xdata(x)
        ys = ax.format_ydata(y)
        zs = ax.format_zdata(z)
        currentpos = [xs, ys, zs]
        currentpos = [float(n.replace(' ', '')) for n in currentpos]
        return 'x=%s, y=%s, z=%s' % (xs, ys, zs)

    ax.format_coord = format_coord

    #######################################################
    # 监测鼠标按下的活动，作高亮点并输出数据
    # print('Activating cursor. Press on figure to show nearest node. Result may not be accurate because of 3D axes.')

    def on_click(event):
        global nodeplot, pressedplot, tmppressed
        # print(currentpos)
        norms = [np.linalg.norm(np.array(k) - np.array(currentpos)) for k in node[:info['nNode']]]
        if max(norms) / min(norms) > 1.3:
            pressed = norms.index(min(norms))
            if pressed == tmppressed:
                return 0
            else:
                tmppressed = pressed
            try:
                if info['mode'] == 'bar':
                    nodeplot.remove()
                pressedplot.remove()
            except Exception as e:
                pass
            if info['mode'] == 'bar':
                X_, Y_, Z_ = map(copy.deepcopy, [X1, Y1, Z1])
                X_, Y_, Z_ = map(list, [X_, Y_, Z_])
                X_.pop(pressed)
                Y_.pop(pressed)
                Z_.pop(pressed)
                nodeplot = ax.scatter(X_, Y_, Z_, color='blue', zorder=0)
            pressedplot = ax.scatter(node[pressed][0], node[pressed][1], node[pressed][2], s=50, color='yellow',
                                     zorder=5)

            print('Node ' + str(pressed + 1) + ' :')
            print('\tCoord: ', node[pressed])
            if displacement != None:
                print('\tNodal displacement: ', displacement[pressed])
            if restrictForce != None:
                print('\tConstraint force: ', restrictForce[pressed])
            rltEle = []
            for i, ele in enumerate(element):
                if pressed + 1 in ele:
                    rltEle.append(i)
            print('\tConnected elements: ')
            for i in rltEle:
                if info['mode'] == 'bar' and stress != None:
                    print('\t\tElement ' + str(i + 1) + '(Node' + str(element[i]) + '), Stress: ' + str(stress[i]))
                elif info['mode'] == 'beam':
                    print('\t\tElement ' + str(i + 1) + '(Node' + str(element[i]) + ')')
            plt.pause(0.01)

    if oncanvas == True:
        F.fig.canvas.mpl_connect('button_press_event', on_click)
    else:
        F.canvas.mpl_connect('button_press_event', on_click)
    #######################################################

    plt.ioff()
    # plt.savefig(root + inpname + '_plot.png', dpi=300)
    plt.show()
    return F


def rdline(line):
    return line.split()


# 路径预处理
def path_preprocess(inp_filepath):
    # 字符串预处理
    st = inp_filepath.replace("\\", "/")
    st = st.split("/")
    root = ""
    for part in range(len(st) - 1):
        root += (st[part])
        root += (r"/")
    # 含后缀名
    filename = st[-1]
    # 不含后缀名
    inpname = filename.split(".")[0]
    return root, filename, inpname


if __name__ == "__main__":
    tkroot = tk.Tk()
    inp_filepath = tk.filedialog.askopenfilename(title="Select input file", filetypes=[('TextFile', '.txt')],
                                                 initialdir=os.getcwd())
    tkroot.destroy()
    run(inp_filepath, oncanvas=False)
