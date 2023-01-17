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

#ifndef SPROB_H
#define SPROB_H

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <string>

using namespace Eigen;

/*
QP definition:

min sum_i 0.5*x[i].'*H[i]*x[i] + g[i].'x[i]

s.t.
    Aeq[i]*x[i] = beq[i]
    Aineq[i]*x[i] <= bineq[i]
    sum_i A[i]x[i] = 0
*/

class sProb
{
    public:
        MatrixXd *H;
        VectorXd *g;
        MatrixXd *Aeq;
        VectorXd *beq;
        MatrixXd *Aineq;
        VectorXd *bineq;
        MatrixXd *A;
        VectorXd *ub;
        VectorXd *lb;
        unsigned int Nagents;

        sProb();                            //default constructor
        sProb(unsigned int Nagents);        //constructor
        ~sProb();                           //destructor
        sProb(const sProb& other);          //copy constructor
        sProb& operator=(const sProb& other);      //operator=

        int csvRead(MatrixXd& outputMatrix, const std::string& fileName, const std::streamsize dPrec);
        void read_sProb(const std::string& folderName, unsigned int Nagents);
        void read_AA(const std::string& folderName, unsigned int Nagents);
        void read_ublb(const std::string& folderName, unsigned int Nagents);
};

#endif
