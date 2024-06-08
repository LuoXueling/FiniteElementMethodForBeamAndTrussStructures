#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <windows.h>
#include <Commdlg.h>
#include <stdio.h>
#include "csm.h"

using namespace std;

//梁单元结构矩阵乘法没搞定，以及其他坐标变换
//eleNodeNo数组在梁单元应该有nEleNode+1=3+1行，第三行是辅助节点，书上说辅助节点最后编号。辅助节点在第一个节点截面的y主轴上

//主函数
int main(int argc, char* argv[]) {
	if (argc == 4) {
		openFilePath = argv[1];
		filepath = argv[2];
		filename = argv[3];
	}
	else {
		getFilename();
	}
	try {
		scan();//调用输入子函数输入各种参数
	}
	catch (std::exception & e) {
		cout << "An error occured when scanning file" << endl;
		cout << e.what() << endl;
		return -1;
	}
	try {
		structureMatrix();//计算结构刚度矩阵
	}
	catch (std::exception & e) {
		cout << "An error occured in structure matrix" << endl;
		cout << e.what() << endl;
		return -2;
	}
	try {
		cholesky(N, K, P);//计算节点位移
	}
	catch (std::exception & e) {
		cout << "An error occured when solving equation" << endl;
		cout << e.what() << endl;
		return -3;
	}
	try {
		internalForce();//计算单元内力
	}
	catch (std::exception & e) {
		cout << "An error occured when solving internal force" << endl;
		cout << e.what() << endl;
		return -4;
	}
	try {
		freeSpace();
	}
	catch (std::exception & e) {
		cout << "An error occured when freeing memory" << endl;
		cout << e.what() << endl;
		return -5;
	}
	return 0;
}

void getFilename() {
	OPENFILENAME ofn;
	char path[100];
	char name[100];
	ZeroMemory(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.hwndOwner = NULL;
	ofn.lpstrFile = (LPSTR)path;
	ofn.lpstrFile[0] = '\0';
	ofn.nMaxFile = sizeof(path);
	ofn.lpstrFilter = (LPCSTR)"Text\0*.txt\0";
	ofn.nFilterIndex = 1;
	ofn.lpstrFileTitle = name;
	ofn.nMaxFileTitle = sizeof(name);
	ofn.lpstrInitialDir = NULL;
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;
	if (GetOpenFileName(&ofn))
	{
		openFilePath = path;
		cout << openFilePath << endl;
		char* buf;
		filename = strtok_s(name, ".", &buf);
		filepath = strtok_s(path, ".", &buf);
	}
	else
	{
		cout << "cancelled" << endl;
	}
}

//数据输人子函数
void scan(){
	CreateDirectory(filepath.c_str(), NULL);

	ifstream fp6(openFilePath);
	ofstream fp5(filepath + "\\" + filename + "_displayInput.txt");

	fp6 >> mode;
	fp5 << "Element type = " << mode << endl;

	fp6 >> nEleNode;
	fp5 << "nEleNode(Number of nodes of each element) = " << nEleNode << endl;

	fp6 >> nNodeFdm;
	fp5 << "nNodeFdm(Number of DOF of each node) = " << nNodeFdm << endl;

	fp6 >> nNode;
	fp5 << "nNode(Number of nodes) = " << nNode << endl;

	fp6 >> nEle;
	fp5 << "nEle(Number of elements) = " << nEle << endl;

	fp6 >> nType;
	fp5 << "nType(Number of types of elements) =" << nType << endl;

	fp6 >> nRestrict;
	fp5 << "nRestrict(Constrained number of DOF) = " << nRestrict << endl;

	allocSpace();

	for (int i = 1; i <= nEle; i++) {
		fp6 >> eleType[i];
		fp5 << "Type of the " << i << " th element = " << eleType[i] << endl;
	}

	for (int i = 1; i <= nType; i++) {
		if (mode == "bar") {
			fp6 >> eleMaterial[1][i] >> eleMaterial[2][i];
			fp5 << "E and A of " << i << " th element" << eleMaterial[1][i] << " " << eleMaterial[2][i] << endl;
		}
		else if (mode == "beam") {
			for (int j = 1; j <= 15; j++) {
				fp6 >> eleMaterial[j][i];
			}
			fp5 << "E, A, K of " << i << " th type of element:" << eleMaterial[1][i] << " " << eleMaterial[2][i] << " " << eleMaterial[9][i] << endl;
			fp5 << "G、Ay、Az of " << i << " th type of element:" << eleMaterial[3][i] << " " << eleMaterial[4][i] << " " << eleMaterial[5][i] << endl;
			fp5 << "Ix、Iy、Iz of " << i << " th type of element:" << eleMaterial[6][i] << " " << eleMaterial[7][i] << " " << eleMaterial[8][i] << endl;
			fp5 << "a1、b1、c1 of " << i << " th type of element:" << eleMaterial[10][i] << " " << eleMaterial[11][i] << " " << eleMaterial[12][i] << endl;
			fp5 << "a2、b2、c2 of " << i << " th type of element:" << eleMaterial[13][i] << " " << eleMaterial[14][i] << " " << eleMaterial[15][i] << endl;
		}
	}

	if (mode == "bar") {
		for (int i = 1; i <= nNode; i++) {
			fp6 >> X[i] >> Y[i] >> Z[i];
			fp5 << "Coords of " << i << " th node: " << X[i] << " " << Y[i] << " " << Z[i] << endl;
		}
	}
	else if (mode == "beam") {
		for (int i = 1; i <= nNode+nEle; i++) {
			fp6 >> X[i] >> Y[i] >> Z[i];
			fp5 << "Coords of " << i << " th node (The support node is at the end): " << X[i] << " " << Y[i] << " " << Z[i] << endl;
		}
	}
	
	for (int i = 1; i <= nEle; i++) {
		fp6 >> eleNodeNo[1][i] >> eleNodeNo[2][i];
		if (mode == "beam") {
			fp6 >> eleNodeNo[3][i];
		}
		fp5 << "ID of node of " << i << " th element: " << eleNodeNo[1][i] << " " << eleNodeNo[2][i] << endl;
		if (mode == "beam") {
			fp5 << "ID of the support node: " << eleNodeNo[3][i] << endl;
		}
	}

	for (int i = 1; i <= nRestrict; i++) {
		fp6 >> resFdmNo[i];
		fp5 << "ID of " << i << " th constrained DOF: " << resFdmNo[i] << endl;
	}

	for (int i = 1; i <= nNode * nNodeFdm; i++) {
		fp6 >> P[i];
		P1[i] = P[i];
		fp5 << "Load of " << i << " th DOF: " << P[i] << endl;
	}
	fp5.close();
	fp6.close();
}

//求每个杆的长度(i表示单元号)
void eleLength(int i) {
	X2 = X[eleNodeNo[2][i]];    X1 = X[eleNodeNo[1][i]];
	Y2 = Y[eleNodeNo[2][i]];    Y1 = Y[eleNodeNo[1][i]];
	Z2 = Z[eleNodeNo[2][i]];    Z1 = Z[eleNodeNo[1][i]];
	L = sqrt((X2 - X1) * (X2 - X1) + (Y2 - Y1) * (Y2 - Y1) + (Z2 - Z1) * (Z2 - Z1));
}

//杆单元的单元刚度阵(单元坐标系下)(i表示单元号)
void stiffnessMatrix_unit(int i) {
	double A, E;
	eleLength(i);
	eleTypeN = eleType[i];
	E = eleMaterial[1][eleTypeN];  A = eleMaterial[2][eleTypeN];
	Ke[1][1] = E * A / L;    Ke[1][2] = -E * A / L;
	Ke[2][1] = -E * A / L;   Ke[2][2] = E * A / L;
	//dispMat(Ke);
}

//坐标转换矩阵(i表示单元号)
void transformMatrix(int i) {
	eleLength(i);
	T[1][1] = (X2 - X1) / L;    T[2][4] = (X2 - X1) / L;
	T[1][2] = (Y2 - Y1) / L;    T[2][5] = (Y2 - Y1) / L;
	T[1][3] = (Z2 - Z1) / L;    T[2][6] = (Z2 - Z1) / L;

	//其余元素已在声明时赋零
	for (int m = 1; m <= KeSize; m++)
		for (int n = 1; n <= (nNodeFdm * nEleNode); n++)
			TT[n][m] = T[m][n];//求出坐标转换阵的转置矩阵
	//dispMat(TT);
}

//总体坐标下的单元刚度阵(矩阵乘法)(i表示单元号)
void multiplyMatrix(int i) {
	double b;
	if (mode == "bar") {
		stiffnessMatrix_unit(i);
		transformMatrix(i);
	}
	else if (mode == "beam") {
		Ke_T_TT_TTF_Matrix_unitbeam(i);
	}
	double s[13][13]; //作矩阵乘法时的中间矩阵

	//TT*Ke
	for (int n = 1; n <= (nNodeFdm * nEleNode); n++) {
		for (int m = 1; m <= KeSize; m++) {//其他单元需要修改区间
			b = 0.0;
			for (int j = 1; j <= KeSize; j++)
				b += TT[n][j] * Ke[j][m];
			s[n][m] = b;
		}
	}

	//TT*Ke*T
	for (int m = 1; m <= (nNodeFdm * nEleNode); m++) {
		for (int n = 1; n <= (nNodeFdm * nEleNode); n++) {
			b = 0.0;
			for (int j = 1; j <= KeSize; j++)//其他单元需要修改区间
				b += s[m][j] * T[j][n];
			KeInRef[m][n] = b;
		}
	}
}

//形成LD数组，长度为nNode * nNodeFdm
void formLD() {
	int IG, J;
	LD[1] = 1;

	//按节点循环，确定与其有关的最小节点号和该点所占行的LD数组
	for (int node = 1; node <= nNode; node++) {
		IG = 100000;

		//按单元循环，确定其最小节点号
		for (int ele = 1; ele <= nEle; ele++)
			for (int j = 1; j <= nEleNode; j++) {
				//判别单元中是否含有node点
				if (eleNodeNo[j][ele] != node)continue;
				//寻找与node点有关的最小节点号放入IG
				for (int L = 1; L <= nEleNode; L++) {
					if (eleNodeNo[L][ele] >= IG)continue;
					IG = eleNodeNo[L][ele];
				}
			}

		//确定node点所含的LD数组
		for (int i = 1; i <= nNodeFdm; i++) {
			//node点所对应的刚阵的每一行
			J = nNodeFdm * (node - 1) + i;
			if (J == 1)continue;//第一行一定从0开始
			LD[J] = LD[J - 1] + nNodeFdm * (node - IG) + i;
		}

		//确定一维数组的总容量
		N = nNode * nNodeFdm;
		lenK = LD[N];
	}
	K = new double[lenK + 1]();
}

//形成eleDispInAll数组
//对于杆单元，单元位移为I(I=1,2,...,6)，它在结构位移d中的编号为eleDispInAll[I]
//对于梁单元，单元位移为I(I=1,2,...,12)，它在结构位移d中的编号为eleDispInAll[I]
void formIS(int i) {
	for (int k=1; k <= nEleNode * nNodeFdm; k++) {
		eleDispInAll[k] = (eleNodeNo[(k - 1) / nNodeFdm + 1][i] - 1) * nNodeFdm + (k - ((k - 1) / nNodeFdm) * nNodeFdm);
	}
}

//组集结构刚度阵
void structureMatrix() {
	int isL, NI, NJ, IJ;
	formLD();

	for (int m = 1; m <= nEle; m++) {
		// 形成该单元在整体坐标下的刚度矩阵
		multiplyMatrix(m);
		// 确定该单元各位移在整体的编号
		formIS(m);
		for (int i = 1; i <= (nNodeFdm * nEleNode); i++) {
			for (int j = 1; j <= (nNodeFdm * nEleNode); j++) {
				isL = eleDispInAll[i] - eleDispInAll[j];
				if (isL >= 0) {//只需要处理下三角
					// NI是该单元第i个节点对应的行 
					NI = eleDispInAll[i];
					// 如i=1，j=2，eleDispInAll[i]=1,eleDispInAll[j]=2
					// 则这一块放在K[LD[NI]+1]处
					IJ = LD[NI] - (NI - eleDispInAll[j]);
					K[IJ] += KeInRef[i][j];
				}
			}
		}
	}

	//约束处理（置大数法）
	for (int i = 1; i <= nRestrict; i++) {
		NI = resFdmNo[i]; NJ = LD[NI];
		K[NJ] = 1e25;
	}
}

void cholesky(int n, double* a, double* x) {
	int ij, ii, jj, ik, jk, kk, iig, igp, jgp, mi, mj, mij;

	for (int i = 1; i <= n; i++) {
		if (i != 1) {//i！=1时进行以下分解（i=1时,d11=k11，无须分解）
			mi = i - (LD[i] - LD[i - 1]) + 1;//第i行最左非零元素列号mi
			if (mi != i) {//第i行不只有对角元的情况
				iig = LD[i] - i;
				for (int j = mi; j <= i - 1; j++) {
					if (j != mi) {
						mj = j - (LD[j] - LD[j - 1]) + 1;//第j行最左非零元素列号mj
						igp = LD[i] - (i - j);//当前被分解元素Kij的一维地址
						if (mj < mi) mij = mi;//求mij=max{mi,mj}
						else mij = mj;
						jgp = LD[j] - j;
						if (mij <= j - 1) {//mij<=j-1时
							for (int k = mij; k <= j - 1; k++) {//式5-43中求和部分的循环
								ik = iig + k; jk = jgp + k; kk = LD[k];//分别计算Lik,Ljk,dkk的一维地址
								a[igp] -= a[ik] * a[kk] * a[jk];//累积计算式(5-43)中第二式括号内部分
							}
						}
					}
					if (j == mi) igp = LD[i - 1] + 1;//最左非零元素的一维地址
					ii = LD[j];//对角元djj的一维地址
					a[igp] = a[igp] / a[ii];//按式(5-43)最后算出Lij
					x[i] -= a[igp] * a[ii] * x[j];//分解载荷项，累积计算式(5-53)括号内部分
				}
				ij = LD[i];
				for (int k = mi; k <= i - 1; k++) {//式(5-44)中求和部分的循环
					ii = iig + k; jj = LD[k];//分别计算Lik和dkk的一维地址
					a[ij] -= a[ii] * a[ii] * a[jj];//按式(5-44)累积计算dii
				}
			}
		}
		ij = LD[i];//对角元dii的一维地址（对i行只有对角元的情况未得出过）
		x[i] = x[i] / a[ij];//按式(5-53)，完成Pi的分解
	}

	for (int i = n; i >= 2; i--) {//回代循环
		mi = i - (LD[i] - LD[i - 1]) + 1;//第i行最左非零元素列号
		if (mi == i) continue;//mi=i在式(5-56)应用范围之外，不作任何处理
		iig = LD[i] - i;
		for (int k = mi; k <= i - 1; k++) {//式5-56中k的循环
			ij = iig + k;//Lik的一维地址
			x[k] -= a[ij] * x[i];//按式5-56计算
		}
	}

	ofstream fp2(filepath + "\\" + filename + "_nodeDisplacement.txt");

	for (int i = 1; i <= n; i++) {
		fp2 << i << " " << x[i] << " " << a[i] << endl;
	}
	fp2.close();
	//dispMat(x);
}

void internalForce(){
	ofstream fp1(filepath + "\\" + filename + "_stress.txt");
	ofstream fp3(filepath + "\\" + filename +"_result.txt");
	ofstream fp4(filepath + "\\" + filename +"_restrictForce.txt");
	for (int i = 1; i <= nEle; i++) {
		//该单元的节点位移（整体坐标系）
		if (mode == "bar") {
			transformMatrix(i);
			stiffnessMatrix_unit(i);
		}
		else if (mode == "beam") {
			Ke_T_TT_TTF_Matrix_unitbeam(i);
		}
		formIS(i);
		for (int k = 1; k <= (nNodeFdm * nEleNode); k++)
			ue[k] = P[eleDispInAll[k]];//P当前存的是节点位移
		fp3 << "Nodal displacement of " << i << " th element in the global coordinate" << endl;
		for (int ii = 1; ii <= nEleNode * nNodeFdm; ii++)
			fp3 << ii << " " << ue[ii] << endl;

		//该单元局部坐标系下节点位移
		//T*ue将整体转为局部
		for (int j = 1; j <= KeSize; j++) {
			for (int k = 1; k <= nEleNode * nNodeFdm; k++)
				eleDispInEle[i][j] += T[j][k] * ue[k]; 
		}

		fp3 << "Nodal displacement of " << i << " the element in the local coordinate" << endl;
		for (int k = 1; k <= nEleNode * nNodeFdm; k++) {
			fp3 << eleDispInEle[i][k] << " ";
		}
		fp3 << endl;
		
		//该单元内力、应力
		//Ke*ue=Fe
		for (int m = 1; m <= KeSize; m++) {
			for (int j = 1; j <= KeSize; j++)
				eleForceInEle[i][m] += Ke[m][j] * eleDispInEle[i][j];
		}
		fp3 << "Nodal force of " << i << " th element in the local coordinate" << endl;
		for (int k = 1; k <= nEleNode * nNodeFdm; k++) {
			fp3 << eleForceInEle[i][k] << " ";
		}
		fp3 << endl;

		if (mode == "bar") {
			sigStress = eleForceInEle[i][2] / eleMaterial[2][eleType[i]];
			fp3 << "Internal forces and stress of " << i << " th element" << endl;
			// 杆单元的第二个节点力是Fn
			fp3 << eleForceInEle[i][2] << " " << sigStress << endl;
			fp1 << sigStress << endl;
		}
		//该单元整体坐标系下的节点力
		for (int mm = 1; mm <= (nEleNode * nNodeFdm); mm++) {
			for (int j = 1; j <= KeSize; j++)
				F[i][mm] += TT[mm][j] * eleForceInEle[i][j];
		}
		fp4 << "Nodal force of " << i << " element in the global coordinate" << endl;
		for (int ii = 1; ii <= nEleNode * nNodeFdm; ii++)
			fp4 << ii << " " << F[i][ii] << endl;
		//求解整体结构节点力
		for (int jj = 1; jj <= nEleNode * nNodeFdm; jj++)
			allNodalForce[eleDispInAll[jj]] = allNodalForce[eleDispInAll[jj]] + F[i][jj];
	}
	//求解约束反力
	for (int k = 1; k <= nNodeFdm * nNode; k++)
		P1[k] = allNodalForce[k] - P1[k];
	for (int m = 1; m <= nNode; m++) {
		for (int n = 1; n <= nNodeFdm; n++)
			fp4 << m << " th node, " << n << " th constraint force " << P1[(m - 1) * nNodeFdm + n] << endl;
	}
	fp3.close();
	fp4.close();
}

void Ke_T_TT_TTF_Matrix_unitbeam(int i){
	double X1, X2, X3, Y1, Y2, Y3, Z1, Z2, Z3, A1, B1, C1, A2, B2, C2;
	double E, A, G, AZ, AY, J, JY, JZ, K, BZ, BY, C, XX, YY, ZZ;
	X1 = X[eleNodeNo[1][i]];
	X2 = X[eleNodeNo[2][i]];
	X3 = X[eleNodeNo[3][i]];
	Y1 = Y[eleNodeNo[1][i]];
	Y2 = Y[eleNodeNo[2][i]];
	Y3 = Y[eleNodeNo[3][i]];
	Z1 = Z[eleNodeNo[1][i]];
	Z2 = Z[eleNodeNo[2][i]];
	Z3 = Z[eleNodeNo[3][i]];
	L = sqrt((X2 - X1) * (X2 - X1) + (Y2 - Y1) * (Y2 - Y1) + (Z2 - Z1) * (Z2 - Z1));
	eleTypeN = eleType[i];
	E = eleMaterial[1][eleTypeN];
	A = eleMaterial[2][eleTypeN];
	G = eleMaterial[3][eleTypeN];
	AY = eleMaterial[4][eleTypeN];
	AZ = eleMaterial[5][eleTypeN];
	J = eleMaterial[6][eleTypeN];
	JY = eleMaterial[7][eleTypeN];
	JZ = eleMaterial[8][eleTypeN];
	K = eleMaterial[9][eleTypeN];
	A1 = eleMaterial[10][eleTypeN];
	B1 = eleMaterial[11][eleTypeN];
	C1 = eleMaterial[12][eleTypeN];
	A2 = eleMaterial[13][eleTypeN];
	B2 = eleMaterial[14][eleTypeN];
	C2 = eleMaterial[15][eleTypeN];
	L = pow(L, 2);
	L = L - pow((B2 - B1), 2) + pow((C2 - C1), 2);
	L = sqrt(L) + A1 - A2;

	Ke[1][1] = E * A / L;
	Ke[7][7] = Ke[1][1];
	Ke[7][1] = -Ke[1][1];
	BZ = 12 * K * E * JY / G / AZ / L / L;
	BY = 12 * K * E * JZ / G / AY / L / L;
	C = 12 * E * JZ / (1 + BY) / L / L / L;
	Ke[2][2] = C;
	Ke[8][8] = C;
	Ke[8][2] = -C;
	Ke[8][6] = -C * L / 2;
	C = 12 * E * JY / (1 + BZ) / L / L / L;
	Ke[3][3] = C;
	Ke[9][9] = C;
	Ke[9][3] = -C;
	Ke[9][5] = C * L / 2;
	Ke[4][4] = G * J / L;
	Ke[10][10] = Ke[4][4];
	Ke[10][4] = -Ke[4][4];
	C = (4 + BZ) * E * JY / (1 + BZ) / L;
	Ke[5][5] = C;
	Ke[11][11] = C;
	C = 6 * E * JY / (1 + BZ) / L / L;
	Ke[5][3] = -C;
	Ke[11][9] = C;
	Ke[11][3] = -C;
	Ke[11][5] = (2 - BZ) * E * JY / (1 + BZ) / L;
	C = (4 + BY) * E * JZ / (1 + BY) / L;
	Ke[6][6] = C;
	Ke[12][12] = C;
	C = 6 * E * JZ / (1 + BY) / L / L;
	Ke[6][2] = C;
	Ke[12][2] = C;
	Ke[12][8] = -C;
	Ke[12][6] = (2 - BZ) * E * JZ / (1 + BZ) / L;
	for (int m = 1; m <= 12; m++)
		for (int n = 1; n <= (m - 1); n++)
			Ke[n][m] = Ke[m][n];
	//形成梁元位移偏心修正矩阵,置于数组beamTT中
	for (int m = 1; m <= 12; m++)
		beamTT[m][m] = 1;
	beamTT[1][5] = -C1;
	beamTT[1][6] = B1;
	beamTT[2][4] = -C1;
	beamTT[2][6] = -A1;
	beamTT[3][4] = -B1;
	beamTT[3][5] = A1;
	beamTT[7][11] = -C2;
	beamTT[7][12] = B2;
	beamTT[8][10] = -C2;
	beamTT[8][12] = -A2;
	beamTT[9][10] = -B2;
	beamTT[9][11] = A2;
	//形成梁元坐标转换矩阵，置于数组T中
	XX = X2 - X1;
	YY = Y2 - Y1;
	ZZ = Z2 - Z1;
	L = sqrt(XX * XX + YY * YY + ZZ * ZZ);
	C = XX / L;
	T[1][1] = C;
	T[4][4] = C;
	T[7][7] = C;
	T[10][10] = C;
	C = YY / L;
	T[1][2] = C;
	T[4][5] = C;
	T[7][8] = C;
	T[10][11] = C;
	C = ZZ / L;
	T[1][3] = C;
	T[4][6] = C;
	T[7][9] = C;
	T[10][12] = C;
	XX = X3 - X1;
	YY = Y3 - Y1;
	ZZ = Z3 - Z1;
	L = sqrt(XX * XX + YY * YY + ZZ * ZZ);
	C = XX / L;
	T[2][1] = C;
	T[5][4] = C;
	T[8][7] = C;
	T[11][10] = C;
	C = YY / L;
	T[2][2] = C;
	T[5][5] = C;
	T[8][8] = C;
	T[11][11] = C;
	C = ZZ / L;
	T[2][3] = C;
	T[5][6] = C;
	T[8][9] = C;
	T[11][12] = C;
	C = T[1][2] * T[2][3] - T[1][3] * T[2][2];
	T[3][1] = C;
	T[6][4] = C;
	T[9][7] = C;
	T[12][10] = C;
	C = T[1][3] * T[2][1] - T[1][1] * T[2][3];
	T[3][2] = C;
	T[6][5] = C;
	T[9][8] = C;
	T[12][11] = C;
	C = T[1][1] * T[2][2] - T[1][2] * T[2][1];
	T[3][3] = C;
	T[6][6] = C;
	T[9][9] = C;
	T[12][12] = C;
	//dispMat(T);
	//形成梁元内力偏心修正矩阵，置于数组beamTTF中
	for (int m = 1; m <= 12; m++)
		beamTTF[m][m] = 1;
	beamTTF[4][2] = -C1;
	beamTTF[4][3] = B1;
	beamTTF[5][1] = C1;
	beamTT[5][3] = -A1;
	beamTTF[10][8] = -C2;
	beamTTF[6][1] = -B1;
	beamTTF[10][9] = -B2;
	beamTT[6][2] = -A1;
	beamTTF[11][7] = C2;
	beamTTF[11][9] = -A2;
	beamTTF[12][7] = -B2;
	beamTTF[12][8] = A2;

	//计算节点坐标系下单元刚度矩阵

	double s[13][13]; //作矩阵乘法时的中间矩阵
	//dispMat(Ke);

	//beamTTT*Ke''
	for (int n = 1; n <= KeSize; n++) {
		for (int m = 1; m <= KeSize; m++) {
			b = 0.0;
			for (int j = 1; j <= KeSize; j++)
				b += beamTT[j][n] * Ke[j][m];
			s[n][m] = b;
		}
	}

	//beamTTT*Ke''*beamTT
	for (int m = 1; m <= KeSize; m++) {
		for (int n = 1; n <= KeSize; n++) {
			b = 0.0;
			for (int j = 1; j <= KeSize; j++)
				b += s[m][j] * beamTT[j][n];
			Ke[m][n] = b;
		}
	}
	//dispMat(Ke);

	for (int m = 1; m <= KeSize; m++)
		for (int n = 1; n <= (nNodeFdm * nEleNode); n++)
			TT[n][m] = T[m][n];//求出坐标转换阵的转置矩阵
}