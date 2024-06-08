#pragma once
#ifndef CSM_H_INCLUDED
#define CSM_H_INCLUDED
#include <iostream>
#include <malloc.h>
#include <string>
using namespace std;

string openFilePath;//完整路径
string filepath;//包含绝对路径+文件名的路径
string filename;//不包含后缀的文件名

string mode;//bar 或者 beam
int KeSize;//单元刚度矩阵大小，并不是总的自由度，因为杆节点有三个自由度
int nEleNode;//单元的节点总数
int nNodeFdm;//单个节点的自由度数
int nNode;//节点总数
int nEle;//单元总数
int nRestrict;//受约束的自由度总数
int nType, eleTypeN;//单元类别总数，单元的类别数
int N;//N=nNode*nNodeFdm
int lenK;//一维存储K的总容量
double* X;
double* Y;
double* Z;//各个节点的三维坐标
double X1, X2, Y1, Y2, Z1, Z2, b;
int** eleNodeNo;//每个单元的节点号,eleNodeNo[1][i]存放第i个单元第一个节点坐标号,eleNodeNo[2][i]存放第i个单元第二个节点坐标号
int* resFdmNo; //约束的位移号
int* LD;
int* eleType;//每个单元的类别
double** eleMaterial;// 存放杆或梁的材料、截面，对杆为E A，对梁为E A G AZ AY J JY JZ K a1 b1 c1 a2 b2 c2
double* P;
double* P1;//节点载荷
double* allNodalForce; // 整体结构节点力
double** Ke; //单元刚度矩阵
double** T; //坐标转换矩阵
double** TT;//坐标转换矩阵的转置
double* K; //整体刚度矩阵
double** KeInRef;//整体坐标系下的单元刚度阵
int* eleDispInAll;//
double L; //杆单元的长度
double sigStress;//单元应力
double* d;//结构位移矩阵
double* ue;//单元位移矩阵
double** eleDispInEle; //局部坐标下的单元位移矩阵
double** eleForceInEle;//局部坐标系下的单元节点力(即单元内力)
double** F;//单元整体坐标系中的节点力
double** I;
double* y;//解方程时用到的L.Y矩阵

//梁单元额外需要的变量
double** beamTT;//梁元位移偏心修正矩阵
double** beamTTF;//梁元内力偏心修正矩阵


void getFilename();
void scan();
void eleLength(int i);
void stiffnessMatrix_unit(int i);
void transformMatrix(int i);
void multiplyMatrix(int i);
void formLD();
void formIS(int i);
void structureMatrix();
void cholesky(int n, double* a, double* x);
void internalForce();
void Ke_T_TT_TTF_Matrix_unitbeam(int i);

void dispMat(double**& obj) {
	int r = _msize(obj) / sizeof(*obj);
	int col = _msize(obj[0]) / sizeof(*obj[0]);
	for (int i = 1; i < r; i++) {
		for (int j = 1; j < col; j++) {
			cout << obj[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void dispMat(double*& obj) {
	int r = _msize(obj) / sizeof(*obj);
	for (int i = 1; i < r; i++) {
		cout << obj[i] << " ";
	}
	cout << endl << endl;
}

void dispMat(int**& obj) {
	int r = _msize(obj) / sizeof(*obj);
	int col = _msize(obj[0]) / sizeof(*obj[0]);
	for (int i = 1; i < r; i++) {
		for (int j = 1; j < col; j++) {
			cout << obj[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void dispMat(int*& obj) {
	int r = _msize(obj) / sizeof(*obj);
	for (int i = 1; i < r; i++) {
		cout << obj[i] << " ";
	}
	cout << endl << endl;
}

void dispAK() {
	int r = _msize(K) / sizeof(*K);
	int flag = 1;
	bool zerosPrinted = false;
	ofstream file(filepath + "_K.csv");
	for (int i = 1; i <= lenK; i++) {
		int nofZeros = flag - (LD[flag] - LD[flag - 1]);
		if (!zerosPrinted) {
			for (int j = 0; flag <= N && j < nofZeros; j++) {
				file << "0" << ",";
			}
			zerosPrinted = true;
		}
		file << K[i] << ",";
		if (LD[flag] == i) {
			file << endl;
			flag++;
			zerosPrinted = false;
		}
	}
	file.close();
}

void allocSpace() {
	if (mode == "bar") {
		KeSize = 2;
	}
	else if (mode == "beam") {
		KeSize = 12;
	}
	N = nNode * nNodeFdm;
	eleType = new int[nEle + 1]();

	if (mode == "bar") {
		eleMaterial = new double* [2 + 1];
		for (int i = 0; i < 2 + 1; i++) {
			eleMaterial[i] = new double[nType + 1]();
		}
	}
	else if (mode == "beam") {
		eleMaterial = new double* [15 + 1];
		for (int i = 0; i < 15 + 1; i++) {
			eleMaterial[i] = new double[nType + 1]();
		}
	}

	if (mode == "bar") {
		X = new double[nNode + 1]();
		Y = new double[nNode + 1]();
		Z = new double[nNode + 1]();
	}
	else if (mode == "beam") {
		X = new double[nNode + nEle + 1]();
		Y = new double[nNode + nEle + 1]();
		Z = new double[nNode + nEle + 1]();
	}
	resFdmNo = new int[nRestrict + 1]();
	P = new double[N + 1]();
	P1 = new double[N + 1]();
	LD = new int[N + 1]();
	//K在Fld之后分配空间
	allNodalForce = new double[N + 1]();
	d = new double[N + 1]();
	eleDispInAll = new int[nNodeFdm * nEleNode+1]();
	ue = new double[nEleNode * nNodeFdm + 1];

	if (mode == "bar") {
		eleNodeNo = new int* [nEleNode + 1];
		for (int i = 0; i < nEleNode + 1; i++) {
			eleNodeNo[i] = new int[nEle + 1]();
		}
	}
	else if (mode == "beam") {
		eleNodeNo = new int* [nEleNode + 2];
		for (int i = 0; i < nEleNode + 2; i++) {
			eleNodeNo[i] = new int[nEle + 1]();
		}
	}

	Ke = new double* [KeSize+1];//杆单元为2+1，其他单元需要修改
	for (int i = 0; i < KeSize + 1; i++) {
		Ke[i] = new double[KeSize + 1]();
	}
	T = new double* [KeSize + 1];
	for (int i = 0; i < KeSize + 1; i++) {
		T[i] = new double[nNodeFdm * nEleNode + 1]();
	}
	TT = new double* [nNodeFdm * nEleNode + 1];
	for (int i = 0; i < nNodeFdm * nEleNode + 1; i++) {
		TT[i] = new double[KeSize + 1]();
	}
	KeInRef = new double* [nNodeFdm * nEleNode+1];
	for (int i = 0; i < nNodeFdm * nEleNode+1; i++) {
		KeInRef[i] = new double[nNodeFdm * nEleNode]();
	}
	eleDispInEle = new double* [nEle + 1];
	for (int i = 0; i < nEle + 1; i++) {
		eleDispInEle[i] = new double[KeSize+1]();
	}
	eleForceInEle = new double* [nEle + 1];
	for (int i = 0; i < nEle + 1; i++) {
		eleForceInEle[i] = new double[nNodeFdm * nEleNode + 1]();
	}
	F = new double* [nEle + 1];
	for (int i = 0; i < nEle + 1; i++) {
		F[i] = new double[N + 1]();
	}
	I = new double* [N + 1];
	for (int i = 0; i < N + 1; i++) {
		I[i] = new double[N + 1]();
	}
	y = new double[N + 1]();
	

	if (mode == "beam") {
		//梁元需要的变量
		beamTT = new double* [12 + 1];
		for (int i = 0; i < 12 + 1; i++) {
			beamTT[i] = new double[12 + 1]();
		}
		beamTTF = new double* [12 + 1];
		for (int i = 0; i < 12 + 1; i++) {
			beamTTF[i] = new double[12 + 1]();
		}
	}
}

void freeSpace() {
	delete eleType;
	delete[] eleMaterial;
	delete X;
	delete Y;
	delete Z;
	delete[] eleNodeNo;
	delete resFdmNo;
	delete P;
	delete P1;
	delete LD;
	delete K;
	delete allNodalForce;
	delete d;
	delete ue;
	delete[] Ke;
	delete[] T;
	delete[] TT;
	delete[] KeInRef;
	delete[] eleDispInEle;
	delete[] eleForceInEle;
	delete[] F;

	delete[] I;
	delete y;

	if (mode == "beam") {
		//梁元需要的变量
		delete[] beamTT;
		delete[] beamTTF;
	}
}

#endif
