/*
  Copyright 2019 Equinor ASA

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdio>
#include <cstring>
#include <sys/time.h>
#include <cmath>
#include <iostream>

#include <opm/simulators/linalg/bda/BlockedMatrix.hpp>

using bda::BlockedMatrix;

namespace bda
{

BlockedMatrix *allocateBlockedMatrix(int Nb, int nnzbs) {
    BlockedMatrix *mat = new BlockedMatrix();

    mat->nnzValues  = new Block[nnzbs];
    mat->colIndices = new int[nnzbs];
    mat->rowPointers = new int[Nb + 1];
    mat->Nb = Nb;
    mat->nnzbs = nnzbs;

    return mat;
}

void freeBlockedMatrix(BlockedMatrix **mat) {
    if (*mat) {
        delete[] (*mat)->nnzValues;
        delete[] (*mat)->colIndices;
        delete[] (*mat)->rowPointers;
        delete (*mat);
        *mat = NULL;
    }
}

BlockedMatrix *soft_copyBlockedMatrix(BlockedMatrix *mat) {
    BlockedMatrix *res = new BlockedMatrix();
    res->nnzValues = mat->nnzValues;
    res->colIndices = mat->colIndices;
    res->rowPointers = mat->rowPointers;
    res->Nb = mat->Nb;
    res->nnzbs = mat->nnzbs;
    return res;
}


/*Sort a row of matrix elements from a blocked CSR-format.*/

void sortBlockedRow(int *colIndices, Block *data, int left, int right) {
    int l = left;
    int r = right;
    int middle = colIndices[(l + r) >> 1];
    double lDatum[BLOCK_SIZE * BLOCK_SIZE];
    do {
        while (colIndices[l] < middle)
            l++;
        while (colIndices[r] > middle)
            r--;
        if (l <= r) {
            int lColIndex = colIndices[l];
            colIndices[l] = colIndices[r];
            colIndices[r] = lColIndex;
            memcpy(lDatum, data + l, sizeof(double) * BLOCK_SIZE * BLOCK_SIZE);
            memcpy(data + l, data + r, sizeof(double) * BLOCK_SIZE * BLOCK_SIZE);
            memcpy(data + r, lDatum, sizeof(double) * BLOCK_SIZE * BLOCK_SIZE);

            l++;
            r--;
        }
    } while (l < r);

    if (left < r)
        sortBlockedRow(colIndices, data, left, r);

    if (right > l)
        sortBlockedRow(colIndices, data, l, right);
}


// LUMat->nnzValues[ik] = LUMat->nnzValues[ik] - (pivot * LUMat->nnzValues[jk]) in ilu decomposition
// a = a - (b * c)
void blockMultSub(Block a, Block b, Block c)
{
    for (int row = 0; row < BLOCK_SIZE; row++) {
        for (int col = 0; col < BLOCK_SIZE; col++) {
            double temp = 0.0;
            for (int k = 0; k < BLOCK_SIZE; k++) {
                temp += b[BLOCK_SIZE * row + k] * c[BLOCK_SIZE * k + col];
            }
            a[BLOCK_SIZE * row + col] -= temp;
        }
    }
}

/*Perform a 3x3 matrix-matrix multiplicationj on two blocks*/

void blockMult(Block mat1, Block mat2, Block resMat) {
    for (int row = 0; row < BLOCK_SIZE; row++) {
        for (int col = 0; col < BLOCK_SIZE; col++) {
            double temp = 0;
            for (int k = 0; k < BLOCK_SIZE; k++) {
                temp += mat1[BLOCK_SIZE * row + k] * mat2[BLOCK_SIZE * k + col];
            }
            resMat[BLOCK_SIZE * row + col] = temp;
        }
    }
}


/* Calculate the inverse of a block. This function is specific for only 3x3 block size.*/

void blockInvert3x3(Block mat, Block res) {
    // code generated by maple, copied from DUNE
    double t4  = mat[0] * mat[4];
    double t6  = mat[0] * mat[5];
    double t8  = mat[1] * mat[3];
    double t10 = mat[2] * mat[3];
    double t12 = mat[1] * mat[6];
    double t14 = mat[2] * mat[6];

    double det = (t4 * mat[8] - t6 * mat[7] - t8 * mat[8] +
                  t10 * mat[7] + t12 * mat[5] - t14 * mat[4]);
    double t17 = 1.0 / det;

    res[0] =  (mat[4] * mat[8] - mat[5] * mat[7]) * t17;
    res[1] = -(mat[1] * mat[8] - mat[2] * mat[7]) * t17;
    res[2] =  (mat[1] * mat[5] - mat[2] * mat[4]) * t17;
    res[3] = -(mat[3] * mat[8] - mat[5] * mat[6]) * t17;
    res[4] =  (mat[0] * mat[8] - t14) * t17;
    res[5] = -(t6 - t10) * t17;
    res[6] =  (mat[3] * mat[7] - mat[4] * mat[6]) * t17;
    res[7] = -(mat[0] * mat[7] - t12) * t17;
    res[8] =  (t4 - t8) * t17;
}


} // end namespace bda