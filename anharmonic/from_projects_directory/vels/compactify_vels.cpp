#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
using namespace std;

void skip_lines(istream&,size_t);

int main(int argc,char** argv){

	if(argc==1){
		cout<<"Usage: Give the filenames on the command line."<<endl;
		return(0);
	}		
	
	string filename(argv[1]);
	string filename2;
	
	if(argc==2){
		filename2=filename+".compact";
	}
	else if(argc==3)
		filename2.assign(argv[2]);
	else{
		cout<<"Usage: Give only the filenames on the command line."<<endl;
		return(0);
	}
	
	cout<<"Writing to file "<<filename2<<"."<<endl;
	
	
	ifstream file1;
	
	double a,b;
	
	file1.open(argv[1]);
	
	string str;
	int t0,t1;
	int N;
	int *ids;
	if(file1.is_open()){ // Initializations, set number of particles etc.
		getline(file1,str);
		file1>>t0;
		cout<<"The initial time step is "<<t0<<"."<<endl;
		getline(file1,str); // The line break
		getline(file1,str);
		file1 >> N;
		cout<<"The number of atoms is "<<N<<"."<<endl;
		ids=new int[N];
		getline(file1,str); // Line break
		skip_lines(file1,5);
		
		for (int i=0;i<N;i++){
			file1>>ids[i];
			getline(file1,str); // Skip the rest of the line
			cout<<ids[i]<<endl;
		}	
		skip_lines(file1,1);
		file1>>t1;
		cout<<"The difference of time steps is "<<t1-t0<<"."<<endl;
	}
	else{
		cout<<"Could not open file "<<filename<<", exiting."<<endl;
		return 0;
	}
	
	// Rewind back to the beginning
	file1.seekg(0,ios::beg);
	
	int id,type;
	double vx,vy,vz;
	
	size_t iter=0;
	
	ofstream file2;
	file2.open(filename2.c_str());
	
	file2<<"Atoms "<<N<<endl;
	file2<<"d_timestep "<<t1-t0<<endl;
	file2<<"Atom ids:"<<endl;
	for(int i=0;i<N;i++){
		file2<<ids[i]<<endl;
	}
	file2<<"----------"<<endl;
	
	while(file1.peek() != char_traits<wchar_t>::eof()){ // The next character is not EOF
		if ((iter)%100000==0)
			cout<<"iter="<<iter<<endl;
		iter++;
		skip_lines(file1,9);

		// Start reading the velocities
		int i;
		for(i=0;i<N;i++){
			file1>>id>>type>>vx>>vy>>vz;
			file2<<vx<<'\n';
			file2<<vy<<'\n';
			file2<<vz<<'\n';
		}

		getline(file1,str); // The line break
		
	}
	
	
	file1.close();
	file2.close();
	return 0;
}

void skip_lines(istream& pStream, size_t pLines)
{
    string s;
    for (; pLines; --pLines)
        getline(pStream, s);
}
