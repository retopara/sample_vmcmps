// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_inv
//! @{


//! immediate inverse of a matrix, storing the result in a dense matrix
template<typename eT>
inline
void
op_inv::apply(Mat<eT>& out, const Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  // no need to check for aliasing, due to:
  // - auxlib::inv() copies A to out before inversion
  // - for 2x2 and 3x3 matrices the code is alias safe
  
  auxlib::inv(out, A);
  }



//! immediate inverse of T1, storing the result in a dense matrix
template<typename T1>
inline
void
op_inv::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const strip_diagmat<T1> strip(X.m);
  
  if(strip.do_diagmat == true)
    {
    op_inv::apply_diag(out, strip.M);
    }
  else
    {
    auxlib::inv(out, X.m);
    }
  }



template<typename T1>
inline
void
op_inv::apply_diag(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const diagmat_proxy_check<T1> A(X.get_ref(), out);
  
  const u32 N = A.n_elem;
  
  out.set_size(N,N);
  
  for(u32 col=0; col<N; ++col)
    {
    for(u32 row=0; row<col; ++row)   { out.at(row,col) = eT(0); }
    
    out.at(col,col) = eT(1) / A[col];
    
    for(u32 row=col+1; row<N; ++row) { out.at(row,col) = eT(0); }
    }
  
  }



//! inverse of T1 (triangular matrices)
template<typename T1>
inline
void
op_inv_tr::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv_tr>& X)
  {
  arma_extra_debug_sigprint();
  
  auxlib::inv_tr(out, X.m, X.aux_u32_a);
  }



//! @}
