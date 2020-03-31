#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<vector>


using namespace std;

ifstream myfile;
ofstream fileout;

string flujo;
string caminos;

vector<pair<int,int>> VP;
 
vector<int> c;
vector<int> b;
int arcos, nodos;

//Lee el archivo de datos inciales
void modif_data() {
	VP = vector<pair<int,int>> (0);
	string s;
	for(int i = 0; i < 5; ++i ){
		getline(myfile,s);
	}
	int h;
	for(h = 1; s.size() > 1; ++h){
		int n = s.size();
		int k = 0;
		pair<int,int> a;
		for(int i = 0; i < n; ++i){
			if(s[i] == '1'){
				if (h == 1)VP.push_back(a);
				VP[k].first = h;
				++k;
			}
			else if(s[i] == '-'){
				if (h == 1)VP.push_back(a);
				VP[k].second = h;
				++k;
				++i;
			}else if(s[i]== '0') {
				if(h == 1)VP.push_back(a);
				++k;
			}
		}
	getline(myfile,s);
	}
	nodos = h-1;

	getline(myfile,s);
	arcos = VP.size();
	b = vector<int> (arcos);
	int i = 6;
	for(int j = 0; j < arcos; ++j){
		b[j] = s[i] - '0';
		i += 3;
	}

	getline(myfile,s);
	c = vector<int> (arcos);
	i = 6;
	for(int j = 0; j < arcos; ++j){
		c[j] = s[i] - '0';
		i += 3;
	}
}

//Genera archivos de datos compatibles con nuestro modelo de AMPL
void generar_datos(){
	fileout.open(flujo);
	fileout << "param n:=" << nodos << ";" << endl;
	fileout << "param: ARCOS: capacidad :=" << endl;
	for(int i = 0; i < arcos ; ++i){
		fileout << VP[i].first << " " << VP[i].second << " " << b[i] << endl;
	}
	fileout.close();
	fileout.open(caminos);
	fileout << "param n:=" << nodos << ";" << endl;
	fileout << "param: ARCOS: coste :=" << endl;
	for(int i = 0; i < arcos ; ++i){
		fileout << VP[i].first << " " << VP[i].second << " " << c[i] << endl;
	}
	fileout.close();
}

//Genera un archivo run que usa los ficheros de datos generados
void generar_run(){
	fileout.open("Practica1.run");
	fileout << "reset;" << endl;
	fileout << "#Problema flujo" << endl;
	fileout << "model fluxmax.mod;" << endl;
	fileout << "data " << flujo << ";" << endl;
	fileout << "option solver cplex;" << endl << "solve;" << endl;
	fileout << "display x;" << endl;
	fileout << "reset;" << endl;
	fileout << "#Problema caminos" << endl;
	fileout << "model caminsminims.mod;" << endl;
	fileout << "data " << caminos << ";" << endl;
	fileout << "option solver cplex;" << endl << "solve;" << endl;
	fileout << "display x;" << endl;
	fileout.close();
}


int main(){
	string nombre;
	cout << "Introduzca nombre del fichero:" << endl;
	//Se refiere al nombre del archivo que contiene los datos del problema
	while (cin >> nombre){

		myfile.open (nombre);
		modif_data();
		cout << endl;

		myfile.close();
		cout << "Introduzca nombre del fichero P.flujo máximo" << endl;
		cin >> flujo;
		cout << "Introduzca nombre del fichero P.caminos mínimos" << endl;
		cin >> caminos;
		generar_datos();
		generar_run();

		//Ejecuta el archivo run des de la terminal de AMPL
		system("./ampl Practica1.run");
		//NOTA: para que esto funcione es necesario que los ficheros empleados y 
		//el ejecutable de AMPL esten en el mismo directorio

		cout << "Introduzca nombre del fichero:" << endl;
	}
}