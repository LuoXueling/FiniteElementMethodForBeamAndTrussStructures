# 生成一个输入文件

from plotStructure import *
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog

info = {'mode': None, 'nEleNode': 0, 'nNodeFdm': 0, 'nNode': 0, 'nEle': 0, 'nType': 0, 'nRestrict': 0, 'eleType': [],
        'eleMaterial': []}
node = []
element = []
restrict = []
force = []
displacement = []
stress = []


def run(oncanvas=True):
    global info, node, element, restrict, force
    print('\n------------Computation type------------\n')
    info['mode'] = input('Type of elements (bar or beam):\n')
    if info['mode'] == 'bar':
        info['nNodeFdm'] = 3
        info['nEleNode'] = 2
    elif info['mode'] == 'beam':
        info['nNodeFdm'] = 6
        info['nEleNode'] = 2
    print('\n------------Node definition------------\n')
    # 输入节点
    while True:
        flag = input(
            'Input node ' + str(info['nNode'] + 1) + ' (3D). \'end\' to stop input,  \'c\' to cancel last input.\n')
        if flag == 'end':
            break
        elif flag == 'c':
            node.pop()
            info['nNode'] -= 1
            print(node)
            continue
        try:
            line = list(map(float, flag.split(',')))
        except Exception as e:
            print('Input error')
            continue
        if len(line) != 3:
            print('Input error')
            continue
        node.append(line)
        info['nNode'] += 1
        try:
            if not oncanvas:
                plot_structure(info, node, oncanvas=False)
            print('Previewing structure')
        except Exception as e:
            print('Input error')
            print(e)
            plt.close()
            continue
    print('\n------------Element definition------------\n')
    # 输入单元
    while True:
        flag = input('Input ID of nodes of element ' + str(
            info['nEle'] + 1) + '. \'end\' to stop input, \'c\' to cancel last input\n')
        if flag == 'end':
            break
        elif flag == 'c':
            element.pop()
            info['nEle'] -= 1
            print(element)
            continue
        try:
            line = list(map(int, flag.split(',')))
        except:
            print('Input error')
            continue
        if len(line) != info['nEleNode']:
            print('Input error')
            continue
        element.append(line)
        info['nEle'] += 1
        try:
            if not oncanvas:
                Fig = plot_structure(info, node, element=element, oncanvas=False)
            print('Previewing structure')
        except:
            print('Input error')
            plt.close()
            continue

    # 输入辅助节点
    if info['mode'] == 'beam':
        i = 0
        while i < info['nEle']:
            flag = input('Input support node of element ' + str(i + 1) + '. \'c\' to stop input.\n')
            if flag == 'c':
                element[i - 1].pop()
                node.pop()
                print(element)
                i -= 1
                continue
            try:
                line = list(map(float, flag.split(',')))
            except:
                print('Input error')
                continue
            if len(line) != 3:
                print('Input error')
                continue
            node.append(line)
            element[i].append(info['nNode'] + i + 1)
            i += 1
    print('\n------------Constraint definition------------\n')
    # 输入约束
    while True:
        flag = input('Input ID of constrained node. \'end\' to stop input, \'c\' to cancel last input.\n')
        if flag == 'end':
            break
        elif flag == 'c':
            restrict.pop()
            info['nRestrict'] -= 1
            print(restrict)
            continue
        try:
            n = int(flag)
            line = ''
            if info['mode'] == 'beam':
                line = input('Input constrained direction (x y z tx ty tz)\n')
            elif info['mode'] == 'bar':
                line = input('Input constrained direction (x y z)。\n')
            line = line.split(',')
            if 'x' in line:
                restrict.append((n - 1) * info['nNodeFdm'] + 1)
                info['nRestrict'] += 1
            if 'y' in line:
                restrict.append((n - 1) * info['nNodeFdm'] + 2)
                info['nRestrict'] += 1
            if 'z' in line:
                restrict.append((n - 1) * info['nNodeFdm'] + 3)
                info['nRestrict'] += 1
            if 'tx' in line:
                restrict.append((n - 1) * info['nNodeFdm'] + 4)
                info['nRestrict'] += 1
            if 'ty' in line:
                restrict.append((n - 1) * info['nNodeFdm'] + 5)
                info['nRestrict'] += 1
            if 'tz' in line:
                restrict.append((n - 1) * info['nNodeFdm'] + 6)
                info['nRestrict'] += 1
        except Exception as e:
            print('Input error')
            continue
        # plot_structure(info, node, element=element, restrict=restrict)
        try:
            if not oncanvas:
                plot_structure(info, node, element=element, restrict=restrict, oncanvas=False)
            print('Previewing structure')
        except Exception as e:
            print('Input error')
            plt.close()
            continue
    print('\n------------Load definition------------\n')
    # 输入力
    force = [[0.0] * info['nNodeFdm'] for i in range(info['nNode'])]
    while True:
        flag = input(
            'Input ID of node where the load is applied. \'end\' to stop input, \'c\' to reset input.\n')
        if flag == 'end':
            break
        elif flag == 'c':
            force = [[0.0] * info['nNodeFdm'] for i in range(info['nNode'])]
            print(force)
            continue
        try:
            n = 0
            n = int(flag)
            line = input('Input load vector\n')
            line = line.split(',')
            if len(line) != info['nNodeFdm']:
                print('Input error')
                continue
            for i, fc in enumerate(line):
                if fc == '0':
                    continue
                else:
                    force[(n - 1)][i] += float(fc)
        except Exception as e:
            print('Input error')
            print(e)
            continue
        try:
            if not oncanvas:
                plot_structure(info, node, element=element, restrict=restrict, force=force, oncanvas=False)
            print('Previewing structure')
        except:
            print('Input error')
            plt.close()
            continue
    print('\n------------Section definition------------\n')
    # 输入截面
    while True:
        if info['mode'] == 'bar':
            flag = input(
                'Input E and A. \'end\' to stop input, \'c\' to reset input.\n')
        elif info['mode'] == 'beam':
            print('Input E,A,G,Ay,Az,Ix,Iy,Iz,K,a1,b1,c1,a2,b2,c2.')
            print('For example: 1,1,1,1,1,1,1,1,0,0,0,0,0,0,0')
            flag = input(
                '\'end\' to stop input, \'c\' to reset input\n')
        if flag == 'end':
            break
        elif flag == 'c':
            info['eleMaterial'].pop()
            info['nType'] -= 1
            print(info['eleMaterial'])
            continue
        try:
            line = list(map(float, flag.split(',')))
        except:
            print('Input error')
            continue
        if info['mode'] == "bar" and len(line) != 2:
            print('Input error')
            continue
        if info['mode'] == "beam" and len(line) != 15:
            print('Input error')
            continue
        info['eleMaterial'].append(line)
        info['nType'] += 1
    print('\n------------Assign sections------------\n')
    # 分配截面
    print('ID of element: ' + str(info['nEle']))
    print('The number of sections: ' + str(len(info['eleMaterial'])))
    print('Information of sections: ', info['eleMaterial'])
    for i in range(info['nEle']):
        n = input(
            'Input ID of section for element ' + str(i + 1) + '. \'all k\' to assign No.k section to all elements.\n')
        if 'all' in n:
            n = n.replace('all', '')
            n = n.replace(' ', '')
            info['eleType'] = [int(n)] * info['nEle']
            break
        info['eleType'].append(int(n))


def save(filename):
    global info, node, element, restrict, force
    with open(filename, 'w')as file:
        file.write(str(info['mode']) + '\n')
        file.write(str(info['nEleNode']) + ' ' + str(info['nNodeFdm']) + ' ' + str(info['nNode']) + ' ' + str(
            info['nEle']) + ' ' + str(info['nType']) + ' ' + str(info['nRestrict']) + '\n')
        for i in info['eleType']:
            file.write(str(i) + ' ')
        file.write('\n')
        for i in info['eleMaterial']:
            for j in i:
                file.write(str(j) + ' ')
            file.write('\n')
        for n in node:
            for j in n:
                file.write(str(j) + ' ')
            file.write('\n')
        for n in element:
            for j in n:
                file.write(str(j) + ' ')
            file.write('\n')
        for i in restrict:
            file.write(str(i) + ' ')
        file.write('\n')
        for i in force:
            for j in i:
                file.write(str(j) + ' ')
        file.write('\n')

        print('\n--------------------------------Structure information--------------------------------\n')
        print("Structure information: ", info)
        print("Coordinates of nodes: ", node)
        print("Nodes of elements: ", element)
        print("Constrained coordinates: ", restrict)
        print("Loads of nodes", force)
        print('\n-----------------------------------------------------------------------\n')
        file.close()

        info = {'mode': None, 'nEleNode': 0, 'nNodeFdm': 0, 'nNode': 0, 'nEle': 0, 'nType': 0, 'nRestrict': 0,
                'eleType': [],
                'eleMaterial': []}
        node = []
        element = []
        restrict = []
        force = []
        displacement = []
        stress = []


if __name__ == "__main__":
    run(oncanvas=False)
    root = tk.Tk()
    filename = tk.filedialog.asksaveasfilename(title="Select file path", filetypes=[('TextFile', '.txt')], initialdir=os.getcwd())
    root.destroy()
    save(filename)
