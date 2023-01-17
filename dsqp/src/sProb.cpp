/*MIT License

Copyright (c) 2023 Goesta Stomberg, Henrik Ebel, Timm Faulwasser, Peter Eberhard

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#include "sProb.hpp"

sProb::sProb(){
    Nagents = 0;
    H = nullptr;
    g = nullptr;
    Aeq = nullptr;
    beq = nullptr;
    Aineq = nullptr;
    bineq = nullptr;
    A = nullptr;
    lb = nullptr;
    ub = nullptr;
}


//Constructor
sProb::sProb(unsigned int nagents)
{
    Nagents = nagents;
    H       = new MatrixXd[Nagents];
    g       = new VectorXd[Nagents];
    Aeq     = new MatrixXd[Nagents];
    beq     = new VectorXd[Nagents];
    Aineq   = new MatrixXd[Nagents];
    bineq   = new VectorXd[Nagents];
    A       = new MatrixXd[Nagents];
    ub      = new VectorXd[Nagents];
    lb      = new VectorXd[Nagents];
}

//Destructor
sProb::~sProb()
{
    delete [] H;
    delete [] g;
    delete [] Aeq;
    delete [] beq;
    delete [] Aineq;
    delete [] bineq;
    delete [] A;
    delete [] ub;
    delete [] lb;
}

//Copy constructor
sProb::sProb(const sProb& other)
{
    Nagents = other.Nagents;
    H       = new MatrixXd[Nagents];
    g       = new VectorXd[Nagents];
    Aeq     = new MatrixXd[Nagents];
    beq     = new VectorXd[Nagents];
    Aineq   = new MatrixXd[Nagents];
    bineq   = new VectorXd[Nagents];
    A       = new MatrixXd[Nagents];
    ub      = new VectorXd[Nagents];
    lb      = new VectorXd[Nagents];

    for (unsigned int i = 0; i < Nagents; i++)
    {
        H[i] = other.H[i];
        g[i] = other.g[i];
        Aeq[i] = other.Aeq[i];
        beq[i] = other.beq[i];
        Aineq[i] = other.Aineq[i];
        bineq[i] = other.bineq[i];
        A[i] = other.A[i];
        ub[i] = other.ub[i];
        lb[i] = other.lb[i];
    }
}

//operator
sProb& sProb::operator=(const sProb& other)
{
    Nagents = other.Nagents;
    H = other.H;
    g = other.g;
    Aeq = other.Aeq;
    beq = other.beq;
    Aineq = other.Aineq;
    bineq = other.bineq;
    A = other.A;
    ub = other.ub;
    lb = other.lb;

    return *this;
}


using Eigen::MatrixXd;
using namespace std;

int sProb::csvRead(MatrixXd& outputMatrix, const std::string& fileName, const std::streamsize dPrec) {
    //This method is by Carlo Cappello, see https://www.youtube.com/watch?v=m118or4f0FE
	ifstream inputData;
	inputData.open(fileName);
	cout.precision(dPrec);
	if (!inputData)
		return -1;
	string fileline, filecell;
	unsigned int prevNoOfCols = 0, noOfRows = 0, noOfCols = 0;
	while (getline(inputData, fileline)) {
		noOfCols = 0;
		stringstream linestream(fileline);
		while (getline(linestream, filecell, ',')) {
			try {
				stod(filecell);
			}
			catch (...) {
				return -1;
			}
			noOfCols++;
		}
		if (noOfRows++ == 0)
			prevNoOfCols = noOfCols;
		if (prevNoOfCols != noOfCols)
			return -1;
	}
	inputData.close();
	outputMatrix.resize(noOfRows, noOfCols);
	inputData.open(fileName);
	noOfRows = 0;
	while (getline(inputData, fileline)) {
		noOfCols = 0;
		stringstream linestream(fileline);
		while (getline(linestream, filecell, ',')) {
			outputMatrix(noOfRows, noOfCols++) = stod(filecell);
		}
		noOfRows++;
	}
	return 0;
}

void sProb::read_sProb(const std::string& folderName, unsigned int Nagents){

    int error;
    MatrixXd tmp;

    for (unsigned int i = 0; i < Nagents; i++){
        std::string fileName = folderName + "/H" + std::to_string(i) + ".csv";
        error = csvRead(H[i],fileName,20);
      
        fileName = folderName + "/g" + std::to_string(i) + ".csv";
        error = csvRead(tmp,fileName,20);
        if (tmp.rows() == 0 || tmp.cols() == 0 || error !=0){
            g[i] = VectorXd::Zero(0);
        } else {
            g[i] = tmp.leftCols(1);
        }

        fileName = folderName + "/Aeq" + std::to_string(i) + ".csv";
        error = csvRead(Aeq[i],fileName,20);
       
        fileName = folderName + "/beq" + std::to_string(i) + ".csv";
        error = csvRead(tmp,fileName,20);
        if (tmp.rows() == 0 || tmp.cols() == 0 || error !=0){
            beq[i] = VectorXd::Zero(0);
        } else {
            beq[i] = tmp.leftCols(1);
        }

        fileName = folderName + "/Aineq" + std::to_string(i) + ".csv";
        error = csvRead(Aineq[i],fileName,20);
        
        fileName = folderName + "/bineq" + std::to_string(i) + ".csv";
        error = csvRead(tmp,fileName,20);
        if (tmp.rows() == 0 || tmp.cols() == 0 || error !=0){
            bineq[i] = VectorXd::Zero(0);
        } else {
            bineq[i] = tmp.leftCols(1);
        }

        fileName = folderName + "/A" + std::to_string(i) + ".csv";
        error = csvRead(A[i],fileName,20);
       
        fileName = folderName + "/ub" + std::to_string(i) + ".csv";
        error = csvRead(tmp,fileName,20);
        if (tmp.rows() == 0 || tmp.cols() == 0 || error !=0){
            ub[i] = VectorXd::Zero(0);
        } else {
            ub[i] = tmp.leftCols(1);
        }        

        fileName = folderName + "/lb" + std::to_string(i) + ".csv";
        error = csvRead(tmp,fileName,20);

        if (tmp.rows() == 0 || tmp.cols() == 0 || error !=0){
            lb[i] = VectorXd::Zero(0);
        } else {
            lb[i] = tmp.leftCols(1);
        }
    }
  

}


void sProb::read_AA(const std::string& folderName, unsigned int Nagents){

    int error;
    MatrixXd tmp;

    for (unsigned int i = 0; i < Nagents; i++){
        std::string fileName = folderName + "/A" + std::to_string(i+1) + ".csv";
        error = csvRead(A[i],fileName,20);
    }   

}

void sProb::read_ublb(const std::string& folderName, unsigned int Nagents){

    int error;
    MatrixXd tmp;

    std::string fileName = "";

    for (unsigned int i = 0; i < Nagents; i++){
        fileName = folderName + "/ub" + std::to_string(i+1) + ".csv";
        error = csvRead(tmp,fileName,20);
        if (tmp.rows() == 0 || tmp.cols() == 0 || error !=0){
            ub[i] = VectorXd::Zero(0);
        } else {
            ub[i] = tmp.leftCols(1);
        }   

        fileName = folderName + "/lb" + std::to_string(i+1) + ".csv";
        error = csvRead(tmp,fileName,20);
        if (tmp.rows() == 0 || tmp.cols() == 0 || error !=0){
            lb[i] = VectorXd::Zero(0);
        } else {
            lb[i] = tmp.leftCols(1);
        }
    }

}