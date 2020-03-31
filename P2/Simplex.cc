#include<iostream>
#include<vector> 
#include<fstream>
#include<iomanip>
#include<string>
#include<algorithm>

using namespace std;

ifstream myfile;
ofstream fileout;

#define FileIn "pm19_exercici_simplex_dades.txt"

using VI = vector<int>;
using MI = vector<VI>;
using VD = vector<double>;
using MD = vector<VD>;
using PI = pair<int,int>;

const double tol = 1e-10;
int alumne;

//Determinamos el conjunto de datos que queremos leer
void leer_conjunto(){
	cout << "Introdueixi el número d'alumme desitjat:" << endl;
	cin >> alumne;
	string s;
	for(int i = 0; i <= (alumne-1)*4; ++i){
		getline(myfile,s);
		while(s[0] != 'P')	getline(myfile,s);
	}
}


void leer_vector(VI& v, string s){
	for(int i = 0; i < s.size() ; ++i){
		if(s[i] != ' '){
			bool neg = false;
			int num = 0;
			while((i < s.size() -1) and s[i] != ' '){
				if(s[i] == '-') neg = true;
				else{
					num *= 10;
					num +=  int (s[i] - '0');
				}
				++i;
			}
			if(neg) num *= -1;
			v.push_back(num);
		}
	}
}


void leer_matriz(MI& M, string s){
	for(int i = 0; s.size() > 2; ++i){
		if(i == M.size()){
			VI v;
			leer_vector(v,s);
			M.push_back(v);
		}else leer_vector(M[i],s);
		getline(myfile,s);
	}

}

//Leemos los datos de entrada
void leer(VI& c, MI& A, VI& b, int& n, int& m){
	string s;
	for(int i = 0; i < 3; ++i) getline(myfile,s);

	//Leemos el vector c
	int i = 0;
	while(s[i] == ' ') ++i;
	if(s[i] == 'C'){
		for(int i = 0; i < 2; ++i) getline(myfile,s);
		leer_vector(c,s);
		for(int i = 0; i < 4; ++i) getline(myfile,s);
	}
	leer_vector(c,s);
	for(int i = 0; i < 3; ++i) getline(myfile,s);


	//Leemos la matriz A
	if(s[i] == 'C'){
		for(int i = 0; i < 2; ++i) getline(myfile,s);
		leer_matriz(A,s);
		for(int i = 0; i < 3; ++i) getline(myfile,s);
	}
	leer_matriz(A,s);
	for(int i = 0; i < 2; ++i) getline(myfile,s);

	//Leemos el vector b 
	if(s[i] == 'C'){
		for(int i = 0; i < 2; ++i) getline(myfile,s);
		leer_vector(b,s);
		for(int i = 0; i < 4; ++i) getline(myfile,s);
	}
	leer_vector(b,s);
	for(int i = 0; i < 2; ++i) getline(myfile,s);

	while(s[0] != 'P') getline(myfile,s);
	
	n = c.size();
	m = A.size();
}


void escribir_iter(const VD& r, const PI& q, const PI& p, const double theta, const double z, const int contador) {
	string coma = ", ";
	fileout << "    Iteració " << contador << " : ";
	fileout <<  "q = " << q.first +1 << coma;
	fileout << "rq = " << r[q.second] << coma;
	fileout << "B(p) = " << p.first +1 << coma;
	fileout << "theta* = " << theta << coma;
	fileout << "z = " << z << endl;
}


void escribir_final(const VD& xB, const VD& r, const VI& vB, const VI& vNB, const double z){
	fileout.setf(ios::fixed);
	fileout.precision(4);
	fileout << "VB* = " << endl << "  ";
	for(const int& i : vB) fileout << i+1 << "  ";
		fileout << endl << endl;
	fileout << "xB* = " << endl << "  ";
	for(const double& i : xB) fileout << i << "  ";
		fileout << endl << endl;
	fileout << "VNB* = " << endl << "  ";
	for(const int& i: vNB) fileout << i+1 << "  ";
		fileout << endl << endl;
	fileout << "r* = " << endl << "  ";
	for(const double& i : r) fileout << i << "  ";
		fileout << endl << endl;
	fileout << "z* = " << endl << "  ";
	fileout << z << endl << endl;
}


int pivot_max(MD& M, int k){
	int n = M.size();
	int maxi = k;
	double maxp = abs(M[k][k]);
	for(int i = k+1; k < n; ++k){
		if(maxp < abs(M[i][k])){
			maxi = i;
			maxp = abs(M[i][k]);
		}
	}
	return maxi;
}

//Esta funcion retorna false si la matriz es singular y true si no
bool invertible(MD M){
	int n = M.size();
	for(int i = 0; i < n; ++i){
		int maxi = pivot_max(M,i);
		if(abs(M[maxi][i]) < tol) return false;
		swap(M[i],M[maxi]);
		for(int j = i+1; j < n; ++j){
			double x = M[j][i]/M[i][i];
			M[j][i] = 0;
			for(int k = i +1; k < n; ++k) M[j][k] -= x*M[i][k];
		}
	}	
	return true;
}


void actualizar_inv(MD& Binv, const VD& dB, const PI& p, const int m){

	MD Id(m, VD(m, 0));
	for(int i = 0; i < m; ++i) Id[i][i] = 1;

	//Creamos H para actualizar la inversa de la base
	MD H(m, VD(m));
	H = Id;
	for(int i = 0; i < m; ++i) H[i][p.second] = -dB[i]/dB[p.second];
	H[p.second][p.second] = -1/dB[p.second];
	
	//Matriz auxiliar para guardar el producto de H*B⁽⁻¹⁾
	MD Maux(m, VD(m, 0));
	for(int i = 0; i < m; ++i) {
		for(int j = 0; j < m; ++j) {
			for(int k = 0; k < m; ++k) {
				Maux[i][j] += H[i][k]*Binv[k][j];				
			}
		}
	}
	
	//Binv es la matriz actualizada de la base
	Binv = Maux;
}

//Calculamos el vector de costes y miramos si la solucion es optima
void calc_1(const MD& Binv, const MI& A, VD& r, const VI& c, const VI& vB, const VI& vNB, const int m, const int nb, bool& optimo){
	VD aux(m,0);
	for(int i = 0; i < m; ++i){
		for(int j = 0; j < m; ++j){
			aux[i] += c[vB[j]]*Binv[j][i];
		}
	}
	optimo = true;
	for(int i = 0; i < nb; ++i){
		r[i] = c[vNB[i]];
		for(int j = 0; j < m; ++j){
			r[i] -= aux[j]* A[j][vNB[i]];
		}
		if(r[i] < 0) optimo = false;
	}
}


void calc_2(const MD& Binv, const MI& A, const VD& r, VD& dB, VD& xB, VI& vB, VI& vNB, PI& p, PI& q, double& theta, double& z, const int m, bool& ilimitado, bool& degenerado){
	//Calculamos la direccion basica y comprobamos si el problema es ilimitado
	dB = VD(m,0);
	ilimitado = true;
	for(int i = 0; i < m; ++i){
		for(int j = 0; j < m; ++j)
			dB[i] += -Binv[i][j]*A[j][q.first];
		if(dB[i] < 0) ilimitado = false;
	}
	if(ilimitado) return;
	//Calculamos la theta y la variable p
	bool first = true;
	for(int i = 0; i < m; ++i){
		if(dB[i] < 0){
			double a = -xB[i]/dB[i];
			if(first){
				theta = a;
				p.second = i;
				first = false;
			}
			else if(a < theta){
				theta = a;
				p.second = i;
			}
		}
	}
	if(theta == 0 )degenerado = true;

	//Actualizamos los vectores vB y vNB
	p.first = vB[p.second];
	vB[p.second] = q.first;
	vNB[q.second] = p.first;
	for(int i = 0; i < m; ++i){
		if(i == p.second) xB[i] = theta;
		else xB[i] += theta*dB[i];
	}
	z += theta*r[q.second];
}


void regla_bland(VD& r, VI& vNB, PI& q, const int nb) {
	bool first = true;
	for(int i = 0; i < nb; ++i) {
		if(r[i] < 0) {
			if(first) {
				q.first = vNB[i];
				q.second = i;
				first = false;
			}
			else if(vNB[i] < q.first) {
				q.first = vNB[i];
				q.second = i;
			}
		}
	}
}


void regla_coste_red(VD& r, VI& vNB, PI& q, const int nb) {
	bool first = true;
	for(int i = 0; i < nb; ++i) {
		if(r[i] < 0) {
			if(first) {
				q.first = vNB[i];
				q.second = i;
				first = false;
			}
			else if(r[i] < r[q.second]) {
				q.first = vNB[i];
				q.second = i;
			}
		}
	}
}

//Hacemos los calculos iniciales de la Fase I 
void FASEI(MD& Binv, MI& A1, VD& xB,VD& r, VI& c1, VI& vB, VI& vNB, const VI& b, double& z, int& nb, const int n, const int m){
 	nb = n;
	c1 = VI(n+m,0);
	for(int i = 0; i < m; ++i){
		for(int j = 0; j < m; ++j){
			if(i == j) A1[i].push_back(1);
			else A1[i].push_back(0);
		}
	}
	xB = VD(b.begin(), b.end());
	vNB = VI(n);
	vB = VI(m);
	z = 0.000;
	for(int i = 0; i < n; ++i)	vNB[i] = i;
	for(int i = 0; i < m; ++i){
		vB[i] = n+i;
		c1[i+n] = 1;
		z += xB[i];
	}

	Binv = MD(m,VD(m,0));
	for(int i = 0; i < m; ++i) Binv[i][i] = 1;
	for(int i = 0; i < m; ++i){
	}
	r = VD(nb);
}

//Hacemos los calculos iniciales de la Fase
void FASEII(MD& Binv, const MI& A, const VD& xB,VD& r, const VI& c, VI& vB, VI& vNB, double& z, int& nb, const int n, const int m, bool degenerado){
	nb = n-m;
	if(degenerado){ //Si hay degeneracion miramos si hay variables artificiales en la base, y si las hay las cambiamos por variables no artificiales
					//Comprobando que la matriz B sea invertible, es decir, que no sea singular
		for(int i = 0; i < m; ++i){
			if(vB[i] >= n){
				bool done = false;
				for(int j = 0; j < n and not done; ++j){
					if(vNB[j] < n){	
						MD auxB(m, VD(m));
						for(int k = 0; k < m and not done; ++k){
							VD auxv;
							if(k == i) auxv = VD(A[vNB[j]].begin(), A[vNB[j]].end());
							else auxv = VD(A[vB[k]].begin(), A[vB[k]].end());
							auxB[k] = auxv;
						}
						if(invertible(auxB)){
							PI cambio;
							cambio.first = vNB[j];
							cambio.second = j;
							swap(vB[i], vNB[j]);
							done = true;
							VD dB(m,0);
							for(int is = 0; is < m; ++is){
								for(int js = 0; js < m; ++js)
									dB[is] += -Binv[is][js]*A[js][vNB[j]];
							}
							actualizar_inv(Binv, dB, cambio, m);
						}
					}
				}
			}
		}
	}
	z = 0;
	for(int i = 0; i < m; ++i){
		int j = vB[i];
		z += xB[i]*c[j];
	}
	sort(vNB.begin(), vNB.end());
	for(int j = 0; j < m; ++j) vNB.pop_back();
	r = VD(nb);
}


int simplex(MD& Binv, const MI& A, VD& dB, VD& xB, VD& r, VI& vB, VI& vNB, const VI& c, const VI& b, double& theta, double& z, const int m, const int n, const int nb, const bool bland, bool& ilimitado, bool& degenerado, int& contador){
	PI p;
	PI q;
	VI base_repetida;
	bool optimo = false;
	calc_1(Binv, A, r, c, vB, vNB, m, nb, optimo);
	while(not optimo){
		if(bland)regla_bland(r, vNB, q, nb);
		else regla_coste_red(r, vNB, q, nb);
		calc_2(Binv, A, r, dB, xB, vB, vNB, p, q, theta, z, m, ilimitado, degenerado);
		if(ilimitado) return 0;
		if(not bland and theta == 0){ //Comprobamos si hay ciclado
			if(base_repetida.size() == 0){
				base_repetida = vB;
				sort(base_repetida.begin(), base_repetida.end());
			}else{
				VI auxv = vB;
				sort(auxv.begin(), auxv.end());
				if(auxv == base_repetida) return -1;	
			}
		}
		actualizar_inv(Binv, dB, p, m);
		escribir_iter(r, q, p, theta, z, contador);
		++ contador;
		calc_1(Binv, A, r, c, vB, vNB, m, nb, optimo);
	}
	PI zero(-1,-1);
	VD rzero(nb,0);
	escribir_iter(rzero, zero, zero, 0, z, contador);
	return contador;
}


int main(){
	MD Binv;
	MI A, A1;
	VD dB, xB, r;
	VI c, c1, b;
	VI vB, vNB;
	int n; 			//numero de variables iniciales
	int nb; 		//numero de variables no basicas
	int m;    		//numero de restricciones == numero de variables basicas
	double z; 		
	double theta;
	bool ilimitado, degenerado;

	myfile.open (FileIn);
	leer_conjunto();
	fileout.open("Resultados.txt");
	fileout << "Conjunt de dades " << alumne << endl << endl;
	for(int i = 1; i <5; ++i){
		fileout << "Problema " << i << endl << endl;

		VI c,b;
		MI A;
		leer(c,A,b,n,m);

		for(int regla = 0; regla < 2; ++regla){
			fileout.setf(ios::fixed);
			fileout.precision(3);
			bool bland = (regla == 0);
			A1 = A;
			fileout << "Inici simplex primal amb regla de ";
			if(regla == 0) fileout << "Bland" << endl;
			else fileout << "costos mínims" << endl << endl;
			fileout << " Fase I" << endl;
			FASEI(Binv, A1, xB, r, c1, vB, vNB, b, z, nb, n, m);
			degenerado = false;
			int contador = 1;
			int iter = simplex(Binv, A1, dB, xB, r, vB, vNB, c1, b, theta, z, m, n, nb, bland, ilimitado, degenerado, contador);
			if(iter == -1) fileout << "    Problema en ciclat" << endl << endl; 
			else if(z > tol) fileout << "    Problema infactible" << endl << endl;
			else{
				fileout << "    Solució bàsica factible trobada, iteració "<< iter << endl << endl;
				fileout << " Fase II" << endl;
				FASEII(Binv, A, xB, r, c, vB, vNB, z, nb, n, m, degenerado);
				ilimitado = false;
				degenerado = false;
				++contador;
				iter = simplex(Binv, A, dB, xB, r, vB, vNB, c, b, theta, z, m, n, nb, bland, ilimitado, degenerado, contador);
				if(ilimitado) fileout << "    Problema il.limitat" << endl << endl;
				else if(iter == -1) fileout << "    Problema en ciclat" << endl << endl;
				else{
					fileout << "    Solució òptima trobada, iteració " << iter << ", z = ";
					fileout.setf(ios::fixed);
					fileout.precision(6);
					fileout << z << endl;
					fileout << "Fi simplex primal" << endl << endl;
					escribir_final(xB, r, vB, vNB, z);
				}
			}
			
		}
	}
}