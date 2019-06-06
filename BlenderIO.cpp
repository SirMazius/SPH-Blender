#include "BlenderIO.h"

unordered_map<string, string> BlenderIO::parametersMap;
vector<string> BlenderIO::l_parameters { "name", "nparts", "particle_volume", "kernel_parts", "tini", "tfin", "tstep", "density", "visco", "stiff", "surften",
		"surften_threshold", "collision_restitution" };
bool BlenderIO::ReadParams(string fileName, Vec3 l_bounds[2]) {
	ifstream myFile;
	string word;

	myFile.open(fileName + ".scn", ios::in);
	if (myFile) {
		cout << "ARCHIVO ENCONTRADO .SCN" << endl;
		while (myFile >> word) {
			//find(l_parameters.begin(), l_parameters.end(), word);
			for (string s : l_parameters)
				if (word.find(s + "=") != std::string::npos) {
					//cout << word.substr(string(s+"=").size(), string(s +"=").size() - word.size()) << endl;
					std::pair<string, string> parameter(s, word.substr(string(s + "=").size(), string(s + "=").size() - word.size()));
					parametersMap.insert(parameter);
				}
		}
		myFile.close();
	} else {
		cout << "ARCHIVO NO ENCONTRADO" << endl;
	}

	myFile.open(fileName + ".prm", ios::in);
	if (myFile) {
		cout << "ARCHIVO ENCONTRADO .PRM" << endl;
		while (myFile >> word) {
			//find(l_parameters.begin(), l_parameters.end(), word);
			for (string s : l_parameters)
				if (word.find(s + "=") != std::string::npos) {
					//cout << word.substr(string(s+"=").size(), string(s +"=").size() - word.size()) << endl;
					std::pair<string, string> parameter(s, word.substr(string(s + "=").size(), string(s + "=").size() - word.size()));
					parametersMap.insert(parameter);
				}
		}
		myFile.close();
	} else {
		cout << "ARCHIVO NO ENCONTRADO" << endl;
		return false;
	}

	int count = 0;
	Vec3 vAux;
	myFile.open(fileName + ".dom", ios::in);
	if (myFile) {
		cout << "ARCHIVO ENCONTRADO .DOM" << endl;
		while (myFile >> word) {
			//find(l_parameters.begin(), l_parameters.end(), word);
			switch (count) {
			case 0:
				l_bounds[0].x = stof(word);
				count++;
				break;
			case 1:
				l_bounds[1].x = stof(word);
				count++;
				break;
			case 2:
				l_bounds[0].y = stof(word);
				count++;
				break;
			case 3:
				l_bounds[1].y = stof(word);
				count++;
				break;
			case 4:
				l_bounds[0].z = stof(word);
				count++;
				break;
			case 5:
				l_bounds[1].z = stof(word);
				break;
			}
		}
		myFile.close();
		return true;
	} else {
		cout << "ARCHIVO NO ENCONTRADO" << endl;
		return false;
	}
}

bool BlenderIO::ReadPOSVEL(string fileName, vector<Vec3> & l_pos, vector<Vec3> & l_velocity) {
	ifstream myFile;
	string word;
	int count = 0;
	int auxCounter = 0;
	int pMax = (FluidParams::nParticles + 1) * 3; //REVISAR
	//int pMax = (FluidParams::nParticles + 1) * 3;
	myFile.open(fileName + ".dat", ios::in);

	if (myFile) {
		cout << "ARCHIVO ENCONTRADO .DAT" << endl;
		Vec3 vAux;
		while (myFile >> word) {

			//cout << stof(word) << endl;
			if (auxCounter != 0 && auxCounter < pMax) {

				switch (count) {
				case 0:
					vAux.x = stof(word);
					count++;
					auxCounter++;
					break;
				case 1:
					vAux.y = stof(word);
					count++;
					auxCounter++;
					break;
				case 2:
					vAux.z = stof(word);
					l_pos.push_back(vAux);
					count = 0;
					auxCounter++;
					break;
				}

			} else if (auxCounter != 0) {

				switch (count) {
				case 0:
					vAux.x = stof(word);
					count++;
					auxCounter++;
					break;
				case 1:
					vAux.y = stof(word);
					count++;
					auxCounter++;
					break;
				case 2:
					vAux.z = stof(word);
					l_velocity.push_back(vAux);
					count = 0;
					auxCounter++;
					break;
				}

			} else {
				auxCounter++;
			}

		}
		myFile.close();
		return true;
	} else {
		cout << "ARCHIVO NO ENCONTRADO" << endl;
		return false;
	}
}

int countDigits(int i) {
	int count = 0;
	while (i > 0) {
		i /= 10;
		count++;
	}
	return count;
}
void BlenderIO::WritePOSVEL(string name, float t, int iteration, vector<Vec3>  l_positions, vector<Vec3>  l_velocity)
//void BlenderIO::WritePOSVEL(string name, float t, int iteration, vector<Vec3> & l_positions, vector<Vec3> & l_velocity)
{

	ofstream myfile;

	int nZeros = 4 - countDigits(iteration);
	string sIteration = string(nZeros, '0').append(to_string(iteration));

	myfile.open(name + "_" + sIteration + ".dat");
	myfile << to_string(t) << "\n";
	for (Vec3 & p : l_positions) {
		myfile << to_string(p.x) << " ";
		myfile << to_string(p.y) << " ";
		myfile << to_string(p.z) << " ";
	}
	myfile << "\n";
	for (Vec3 & v : l_velocity) {
		myfile << to_string(v.x) << " ";
		myfile << to_string(v.y) << " ";
		myfile << to_string(v.z) << " ";
	}
	myfile.close();
	cout << "WritePOSVEL" << endl;
}

void BlenderIO::WriteExcelData(string fileName, float deepDensity, float surfaceDensity, float executionTime) {

	ofstream myFileOutput;
	ifstream myFileInput;
	myFileInput.open(fileName + ".csv");

	if (!myFileInput) {
		myFileInput.close();
		myFileOutput.open(fileName + ".csv");
//		cout << "CREAMOS" << endl;
		myFileOutput << deepDensity << "," << surfaceDensity << "\n";
		myFileOutput.close();
	} else {
//		cout << "Ya está creado" << endl;
		myFileOutput.open(fileName + ".csv", fstream::app);
		myFileOutput << deepDensity << "," << surfaceDensity << "\n";
		myFileOutput.close();
	}
}

void BlenderIO::WriteHeightDensityData(vector<Vec3> l_positions, vector<float> l_density, int iteration) {
	ofstream myfile;
	int size = l_positions.size();
	const string fileName = "HeightDensity";

	myfile.open(fileName + to_string(iteration) + ".csv");

	for (int i = 0; i < size; i++) {
		myfile << to_string(l_positions.at(i).x) << "," << to_string(l_positions.at(i).y) << "," << to_string(l_positions.at(i).z) << "," << l_density.at(i) << "\n";
	}

	myfile.close();
}

void BlenderIO::WriteDensityPerParticle(vector<Vec3> l_positions, vector<float> l_density, int iteration) {
	ofstream myFile;
	string fileName = "particleDensity";
	int size = l_positions.size();
	myFile.open(fileName + to_string(iteration) + ".csv");

	for (int i = 0; i < size; i++) {
		myFile << to_string(l_density.at(i)) << ",";
	}

	myFile.close();

}

void BlenderIO::WriteEigenVectorsNEigenValues(const vector<vector<int>> & l_neighbors, const vector<Vec3> & l_positions, int iteration) {
	ofstream myFile;
	string fileName = "eVs";
	myFile.open(fileName + to_string(iteration) + ".csv");
	int count = l_positions.size();
	float h = FluidParams::kernelRadius;
	for (int i = 0; i < count; i++) {
		vector<Vector3d> l_points;

		for (int j : l_neighbors.at(i)) {
			Vec3 vAux;
			Vec3::vDirector(l_positions.at(j), l_positions.at(i), vAux);
			if (vAux.mag() <= h) {
				Vector3d point(l_positions.at(j).x, l_positions.at(j).y, l_positions.at(j).z);
				l_points.push_back(point);
			}
		}

		int auxCounter = l_points.size();
		MatrixXd points(auxCounter, 3);

		for (int j = 0; j < auxCounter; j++)
			points.row(j) = l_points.at(j);

		MatrixXd a = points.transpose();
		VectorXd mean = points.colwise().mean();
		MatrixXd var = (a).colwise() - mean;
		MatrixXd t = var * var.transpose();
		MatrixXd dst = t * (1 / (a.cols() - 1)); //S
	//
		EigenSolver<MatrixXd> es(dst);
		VectorXd eigenValues = es.eigenvalues().real();
		MatrixXd eigenVectors = es.eigenvectors().real();
//		myFile << to_string(eigenValues[0]) << "," << to_string(eigenValues[1]) << "," << to_string(eigenValues[2]) << to_string(eigenVectors[0]) << "," << to_string(eigenVectors[1]) << "," << to_string(eigenVectors[2]) << "\n";
	}
	myFile.close();
}

void BlenderIO::WriteForces(vector<Vec3> & l_internalForce, vector<Vec3> & l_pressureForce,
			vector<Vec3> & l_externalForce, int iteration) {
	ofstream myFile;
	string fileName = "forces";
	myFile.open(fileName + to_string(iteration) + ".csv");
	int count = l_internalForce.size();

	for (int i = 0; i < count; i++) {
		myFile << l_internalForce.at(i).x << "," << l_internalForce.at(i).y << "," << l_internalForce.at(i).z << ","
				<< l_pressureForce.at(i).x << "," << l_pressureForce.at(i).y << "," << l_pressureForce.at(i).z << ","
				<< l_externalForce.at(i).x << "," << l_externalForce.at(i).y << "," << l_externalForce.at(i).z <<  iteration << "\n";
	}
	myFile.close();
}

