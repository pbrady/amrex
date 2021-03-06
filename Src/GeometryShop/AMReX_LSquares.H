/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */

#ifndef _LSQUARES_H_
#define _LSQUARES_H_

#include <vector>
#include "AMReX_REAL.H"
#include "AMReX_RealVect.H"

namespace amrex
{
  ///
  /**
     This class solves least squares problems
  */
  class LSquares
  {
  public:
    void LeastSquares(Real** A,
                      std::vector<Real>&x,
                      const std::vector<Real>&rhs);

    /// gaussian elimination with partial pivoting
    int gaussElim(Real**A,
                  std::vector<Real>& rhs);

    void swapRows(Real** A,
                  const int& rowi,
                  const int& rowj,
                  const int&numberOfCols);

    void swapRows(std::vector<Real>& rhs,
                  const int& currRow,
                  const int& pivot);

    int findPivot(Real** A,
                  const int& currCol,
                  const int& currRow,
                  const int& numRows,
                  int& pivot);

    void addRows(std::vector<Real>& rhs,
                 const int& rowi,
                 const Real& alpha,
                 const int& rowj);

    void addRows(Real** A,
                 const int& rowi,
                 const Real& alpha,
                 const int& rowj,
                 const int& numberOfCols);

    void timesBeta(std::vector<Real>&rhs,
                   const int& currRow,
                   const Real& Beta);

    void timesBeta(Real** A,
                   const int& rowi,
                   const Real& Beta,
                   const int& numberOfcols);

    void transpose(Real** a_A,
                   Real ** a_Atrans,
                   const int& a_numRows,
                   const int& a_numCols);

    void matMul(Real** a_A,
                Real** a_B,
                Real** a_C,
                const int& a_numRowsA,
                const int& a_numColsA,
                const int& a_numColsB);

    void backSolve(Real** a_A,
                   const std::vector<Real>& a_rhs,
                   const int& a_numArows,
                   std::vector<Real>& a_x);

    void AtimesX(Real** A,
                 const std::vector<Real>&x,
                 const int& numRowsA,
                 std::vector<Real>& Ax);

    void allocArray(const int& rows,
                    const int& cols,
                    Real**& A);

    void freeArray(const int& rows,
                   const int& cols,
                   Real**& A);

  };
}

#endif
