/*This isn't the best way to do this, but it works well enough and we need it only one time */
#include <iostream>
#include <fstream>
#include <filesystem>
int main() {
	
	double estimate;
	double average;
	double error;
	double T;
	
	std::ofstream outputMag("outputMagTemp.csv");
	outputMag<<"temp,finalMag,err"<<std::endl;
	
	std::ofstream outputEnergy("outputEnergyTemp.csv");
	outputEnergy<<"temp,finalEnergy,err"<<std::endl;
	
	std::ofstream outputHeat("outputHeatTemp.csv");
	outputHeat<<"temp,finalHeat,err"<<std::endl;
	
	std::ofstream outputChi("outputChiTemp.csv");
	outputChi<<"temp,finalChi,err"<<std::endl;
	
	
	
	
	
	for(int i =0; i<16;i++){
		T=(double)i/10.0+0.5;
		std::ofstream output("input.dat");
		output<<T<<std::endl;
		output<<"50"<<std::endl;
		output<<"1.0"<<std::endl;
		output<<"0.00"<<std::endl;
		output<<"0"<<std::endl;
		output<<"20"<<std::endl;
		output<<"1000"<<std::endl;
		output<<"0"<<std::endl;


		output<<"  ReadInput << temp;"<<std::endl;
		output<<"  ReadInput << nspin;"<<std::endl;
		output<<"  ReadInput << J;"<<std::endl;
		output<<"  ReadInput << h;"<<std::endl;
		output<<"  ReadInput << metro;"<<std::endl;
		output<<"  ReadInput << nblk;"<<std::endl;
		output<<"  ReadInput << nstep;"<<std::endl;
		output<<"  ReadInput << reset"<<std::endl;
		
		
		
		system("del output.mag.0");
		system("del output.ene.0");
		system("del output.chi.0");
		system("del output.heat.0");
		system("Monte_Carlo_ISING_1D");
		
		std::ifstream inputMag("output.mag.0");
		double data;
		
		
		while(inputMag>>data){
			inputMag>>estimate;
			inputMag>>average;
			inputMag>>error;
		}
		std::cout<<"Temperatura: "<<T<<std::endl;
		outputMag<<T<<","<<average<<","<<error<<std::endl;
		
		
		std::ifstream inputEnergy("output.ene.0");
		
		while(inputEnergy>>data){
			inputEnergy>>estimate;
			inputEnergy>>average;
			inputEnergy>>error;
		}
		std::cout<<"Temperatura: "<<T<<std::endl;
		outputEnergy<<T<<","<<average<<","<<error<<std::endl;
		
		
		std::ifstream inputChi("output.chi.0");
		
		while(inputChi>>data){
			inputChi>>estimate;
			inputChi>>average;
			inputChi>>error;
		}
		std::cout<<"Temperatura: "<<T<<std::endl;
		outputChi<<T<<","<<average<<","<<error<<std::endl;
		
		
		
		std::ifstream inputHeat("output.heat.0");
		
		while(inputHeat>>data){
			inputHeat>>estimate;
			inputHeat>>average;
			inputHeat>>error;
		}
		std::cout<<"Temperatura: "<<T<<std::endl;
		outputHeat<<T<<","<<average<<","<<error<<std::endl;
		
		
		
		
		inputMag.close();
		inputEnergy.close();
		inputHeat.close();
		inputChi.close();
	}
    return 0;
}
