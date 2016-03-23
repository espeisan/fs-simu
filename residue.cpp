#include "common.hpp"

#define CONTRACT1(i,       size_i)                         for (int i = 0; i < size_i; ++i)
#define CONTRACT2(i,j,     size_i, size_j)                 for (int i = 0; i < size_i; ++i) for (int j = 0; j < size_j; ++j)
#define CONTRACT3(i,j,k,   size_i, size_j, size_k)         for (int i = 0; i < size_i; ++i) for (int j = 0; j < size_j; ++j) for (int k = 0; k < size_k; ++k)
#define CONTRACT4(i,j,k,l, size_i, size_j, size_k, size_l) for (int i = 0; i < size_i; ++i) for (int j = 0; j < size_j; ++j) for (int k = 0; k < size_k; ++k) for (int l = 0; l < size_l; ++l)
 


int epsilon(int i, int j, int k)  // permutation function
{
  if(i==1 && j==2 && k==3) return  1;
  if(i==2 && j==3 && k==1) return  1;
  if(i==3 && j==1 && k==2) return  1;
  if(i==3 && j==2 && k==1) return -1;
  if(i==1 && j==3 && k==2) return -1;
  if(i==2 && j==1 && k==3) return -1;
return 0;
}

template<class D, class T>
D determinant(T const& a, int dim)  //determinant
{
  if (dim==1)
    return a(0,0);
  else
  if (dim==2)
    return a(0,0)*a(1,1)-a(0,1)*a(1,0);
  else
  if (dim==3)
    return a(0,0)*(a(1,1)*a(2,2)-a(1,2)*a(2,1))+a(0,1)*(a(1,2)*a(2,0)-a(1,0)*a(2,2))+a(0,2)*(a(1,0)*a(2,1)-a(1,1)*a(2,0));
  else
  {
    printf("double determinant(Tensor const& a, int dim): invalid dim, get %d\n", dim);
    throw;
  }
}

template<class TensorType, class Double>  //inverse matrix
void invert_a(TensorType & a, int dim)
{
  if (dim==1)
  {
    a(0,0)=1./a(0,0);
  }
  else
  if (dim==2)
  {
    Double const det = a(0,0)*a(1,1)-a(0,1)*a(1,0);

    Double const inv00 = a(1,1)/det;
    Double const inv01 = -a(0,1)/det;
    Double const inv10 = -a(1,0)/det;
    Double const inv11 = a(0,0)/det;

    a(0,0) = inv00;
    a(0,1) = inv01;
    a(1,0) = inv10;
    a(1,1) = inv11;
  }
  else if (dim==3)
  {
    Double const det = a(0,0)*(a(1,1)*a(2,2)-a(1,2)*a(2,1))+a(0,1)*(a(1,2)*a(2,0)-a(1,0)*a(2,2))+a(0,2)*(a(1,0)*a(2,1)-a(1,1)*a(2,0));

    Double const inv00 = ( a(1,1)*a(2,2)-a(1,2)*a(2,1) )/det;
    Double const inv01 = ( a(0,2)*a(2,1)-a(0,1)*a(2,2) )/det;
    Double const inv02 = ( a(0,1)*a(1,2)-a(0,2)*a(1,1) )/det;
    Double const inv10 = ( a(1,2)*a(2,0)-a(1,0)*a(2,2) )/det;
    Double const inv11 = ( a(0,0)*a(2,2)-a(0,2)*a(2,0) )/det;
    Double const inv12 = ( a(0,2)*a(1,0)-a(0,0)*a(1,2) )/det;
    Double const inv20 = ( a(1,0)*a(2,1)-a(1,1)*a(2,0) )/det;
    Double const inv21 = ( a(0,1)*a(2,0)-a(0,0)*a(2,1) )/det;
    Double const inv22 = ( a(0,0)*a(1,1)-a(0,1)*a(1,0) )/det;

    a(0,0) = inv00;
    a(0,1) = inv01;
    a(0,2) = inv02;
    a(1,0) = inv10;
    a(1,1) = inv11;
    a(1,2) = inv12;
    a(2,0) = inv20;
    a(2,1) = inv21;
    a(2,2) = inv22;

  }
  else
  {
    printf("invalid dim, try to run again dumb \n");
    throw;
  }


}

template <typename Derived>
void getProjectorMatrix(MatrixBase<Derived> & P, int n_nodes, int const* nodes, Vec const& Vec_x_, double t, AppCtx const& app)
{
  int const dim = app.dim;
  Mesh const* mesh = &*app.mesh;
  //DofHandler const* dof_handler = &*app.dof_handler;
  std::vector<int> const& dirichlet_tags  = app.dirichlet_tags;
  //std::vector<int> const& neumann_tags    = app.neumann_tags  ;
  //std::vector<int> const& interface_tags  = app.interface_tags;
  std::vector<int> const& solid_tags      = app.solid_tags    ;
  std::vector<int> const& triple_tags     = app.triple_tags   ;
  //std::vector<int> const& periodic_tags   = app.periodic_tags ;
  std::vector<int> const& feature_tags    = app.feature_tags  ;
  //Vec const& Vec_normal = app.Vec_normal;

  P.setIdentity();

  Tensor I(dim,dim);
  Tensor Z(dim,dim);
  Vector X(dim);
  Vector normal(dim);
  int    dofs[dim];
  int    tag;
  Point const* point;

  I.setIdentity();
  Z.setZero();

  // NODES
  for (int i = 0; i < n_nodes; ++i)
  {
    point = mesh->getNodePtr(nodes[i]);
    tag = point->getTag();
    //m = point->getPosition() - mesh->numVerticesPerCell();
    //cell = mesh->getCellPtr(point->getIncidCell());

    if (is_in(tag,feature_tags))
    {
      app.getNodeDofs(&*point, DH_MESH, VAR_M, dofs);
      VecGetValues(Vec_x_, dim, dofs, X.data());
      P.block(i*dim,i*dim,dim,dim)  = feature_proj(X,t,tag);
      continue;
    }
    else
    if (is_in(tag,solid_tags) || is_in(tag, triple_tags))
    {
      app.getNodeDofs(&*point, DH_MESH, VAR_M, dofs);
      VecGetValues(Vec_x_, dim, dofs, X.data());

      normal = solid_normal(X,t,tag);

      P.block(i*dim,i*dim,dim,dim)  = I - normal*normal.transpose();
      //P.block(i*dim,i*dim,dim,dim) = Z;

    }
    else
    if (is_in(tag,dirichlet_tags))
    {
      P.block(i*dim,i*dim,dim,dim) = Z;
    }


  } // end nodes
} // end getProjectorMatrix



// ******************************************************************************
//                            FORM FUNCTION
// ******************************************************************************
PetscErrorCode AppCtx::formFunction(SNES /*snes*/, Vec Vec_up_k, Vec Vec_fun)
{
  
  double utheta = AppCtx::utheta;
  
  if (is_bdf2)
  {
    if (time_step == 0)
      if(!is_bdf_euler_start)
        utheta = 0.5;
  }
  else if (is_bdf3)
  {
    if (time_step <= 1)
      utheta = 0.5;
  }
  
  bool const compact_bubble = false; // eliminate bubble from convective term
  

  //PetscErrorCode      ierr;

  int null_space_press_dof=-1;

  int iter;

  SNESGetIterationNumber(snes,&iter);

  if (!iter)
  {
    converged_times=0;
  }

  if (force_pressure && (iter<2))
  {
    Vector X(dim);
    Vector X_new(dim);
    if (behaviors & BH_Press_grad_elim)
    {
      cell_iterator cell = mesh->cellBegin();
      dof_handler[DH_UNKS].getVariable(VAR_P).getCellDofs(&null_space_press_dof, &*cell);
      // fix the initial guess
      VecSetValue(Vec_up_k, null_space_press_dof, 0.0, INSERT_VALUES);
    }
    else
    {
      point_iterator point = mesh->pointBegin();
      while (!mesh->isVertex(&*point))
        ++point;
      int x_dofs[3];
      dof_handler[DH_MESH].getVariable(VAR_M).getVertexDofs(x_dofs, &*point);
      VecGetValues(Vec_x_1, dim, x_dofs, X_new.data());
      VecGetValues(Vec_x_0, dim, x_dofs, X.data());
      X = .5*(X+X_new);
      dof_handler[DH_UNKS].getVariable(VAR_P).getVertexDofs(&null_space_press_dof, &*point);
      // fix the initial guess
      VecSetValue(Vec_up_k, null_space_press_dof, p_exact(X,current_time+.5*dt,point->getTag()), INSERT_VALUES);
    }

    Assembly(Vec_up_k);

  }

  // checking:
  if (null_space_press_dof < 0 && force_pressure==1 && (iter<2))
  {
    cout << "force_pressure: something is wrong ..." << endl;
    throw;
  }


  Mat *JJ = &Mat_Jac;



  //PetscErrorCode      ierr;

  VecZeroEntries(Vec_fun);
  MatZeroEntries(*JJ);

  // LOOP NAS CÉLULAS Parallel (uncomment it)
//#ifdef FEP_HAS_OPENMP
//  FEP_PRAGMA_OMP(parallel default(none) shared(Vec_up_k,Vec_fun,cout,null_space_press_dof,JJ,utheta,iter))
//#endif
  {
    VectorXd          FUloc(n_dofs_u_per_cell);  // U subvector part of F
    VectorXd          FPloc(n_dofs_p_per_cell);

    /* local data */
    int                 tag;
    MatrixXd            u_coefs_c_mid_trans(dim, n_dofs_u_per_cell/dim);  // n+utheta  // trans = transpost
    MatrixXd            u_coefs_c_old(n_dofs_u_per_cell/dim, dim);        // n
    MatrixXd            u_coefs_c_old_trans(dim,n_dofs_u_per_cell/dim);   // n
    MatrixXd            u_coefs_c_new(n_dofs_u_per_cell/dim, dim);        // n+1
    MatrixXd            u_coefs_c_new_trans(dim,n_dofs_u_per_cell/dim);   // n+1

    MatrixXd            du_coefs_c_old(n_dofs_u_per_cell/dim, dim);        // n
    MatrixXd            du_coefs_c_old_trans(dim,n_dofs_u_per_cell/dim);   // n
    
    MatrixXd            du_coefs_c_vold(n_dofs_u_per_cell/dim, dim);        // n-1
    MatrixXd            du_coefs_c_vold_trans(dim,n_dofs_u_per_cell/dim);   // n-1

    MatrixXd            v_coefs_c_mid(nodes_per_cell, dim);        // mesh velocity; n
    MatrixXd            v_coefs_c_mid_trans(dim,nodes_per_cell);   // mesh velocity; n

    VectorXd            p_coefs_c_new(n_dofs_p_per_cell);  // n+1
    VectorXd            p_coefs_c_old(n_dofs_p_per_cell);  // n
    VectorXd            p_coefs_c_mid(n_dofs_p_per_cell);  // n

    //MatrixXd            x_coefs_c_mid(nodes_per_cell, dim);       // n+utheta
    MatrixXd            x_coefs_c_mid_trans(dim, nodes_per_cell); // n+utheta
    MatrixXd            x_coefs_c_new(nodes_per_cell, dim);       // n+1
    MatrixXd            x_coefs_c_new_trans(dim, nodes_per_cell); // n+1
    MatrixXd            x_coefs_c_old(nodes_per_cell, dim);       // n
    MatrixXd            x_coefs_c_old_trans(dim, nodes_per_cell); // n

    Tensor              F_c_mid(dim,dim);       // n+utheta
    Tensor              invF_c_mid(dim,dim);    // n+utheta
    Tensor              invFT_c_mid(dim,dim);   // n+utheta

    Tensor              F_c_old(dim,dim);       // n
    Tensor              invF_c_old(dim,dim);    // n
    Tensor              invFT_c_old(dim,dim);   // n

    Tensor              F_c_new(dim,dim);       // n+1
    Tensor              invF_c_new(dim,dim);    // n+1
    Tensor              invFT_c_new(dim,dim);   // n+1


    //Tensor              F_c_new(dim,dim);       // n+1
    //Tensor              invF_c_new(dim,dim);    // n+1
    //Tensor              invFT_c_new(dim,dim);   // n+1
    //Tensor              F_c_old(dim,dim);       // n
    //Tensor              invF_c_old(dim,dim);    // n
    //Tensor              invFT_c_old(dim,dim);   // n

    /* All variables are in (n+utheta) by default */

    MatrixXd            dxphi_c(n_dofs_u_per_cell/dim, dim);
    MatrixXd            dxphi_c_new(dxphi_c);
    MatrixXd            dxpsi_c(n_dofs_p_per_cell, dim);
    MatrixXd            dxpsi_c_new(dxpsi_c);
    MatrixXd            dxqsi_c(nodes_per_cell, dim);
    Vector              dxbble(dim);
    Vector              dxbble_new(dim);
    Tensor              dxU(dim,dim);   // grad u
    Tensor              dxU_old(dim,dim);   // grad u
    Tensor              dxU_new(dim,dim);   // grad u
    Tensor              dxUb(dim,dim);  // grad u bble
    Vector              dxP_new(dim);   // grad p
    Vector              Xqp(dim);
    Vector              Xqp_old(dim);
    Vector              Xc(dim);  // cell center; to compute CR element
    Vector              Uqp(dim);
    Vector              Ubqp(dim); // bble
    Vector              Uqp_old(dim);  // n
    Vector              Uqp_new(dim);  // n+1
    Vector              dUqp_old(dim);  // n
    Vector              dUqp_vold(dim);  // n
    Vector              Vqp(dim);
    Vector              Uconv_qp(dim);
    Vector              dUdt(dim);
    double              Pqp_new;
    double              Pqp;
    double              bble_integ=0;
    //VectorXd          FUloc(n_dofs_u_per_cell); // subvetor da função f (parte de U)
    //VectorXd          FPloc(n_dofs_p_per_cell);     // subvetor da função f (parte de P)
    VectorXi            cell_nodes(nodes_per_cell);
    double              J_mid;
    double              J_new, J_old;
    double              JxW_mid;
    double              JxW_new, JxW_old;
    double              weight;
    double              visc=-1; // viscosity
    double              cell_volume;
    double              hk2;
    double              tauk=0;
    double              delk=0;
    double              delta_cd;
    double              rho;

    MatrixXd            Aloc(n_dofs_u_per_cell, n_dofs_u_per_cell);
    MatrixXd            Gloc(n_dofs_u_per_cell, n_dofs_p_per_cell);
    MatrixXd            Dloc(n_dofs_p_per_cell, n_dofs_u_per_cell);
    MatrixXd            Eloc(n_dofs_p_per_cell, n_dofs_p_per_cell);   // GSL, BC
    MatrixXd            Cloc(n_dofs_u_per_cell, n_dofs_p_per_cell);   // GSL
    Tensor              iBbb(dim, dim);                               // BC, i : inverse ..it is not the inverse to CR element
    MatrixXd            Bbn(dim, n_dofs_u_per_cell);                  // BC
    MatrixXd            Bnb(n_dofs_u_per_cell, dim);                  // BC
    MatrixXd            Dpb(n_dofs_p_per_cell, dim);                  // BC
    MatrixXd            Gbp(dim, n_dofs_p_per_cell);                  // BC
    MatrixXd            Gnx(n_dofs_u_per_cell, dim);                  // CR ;; suffix x means p gradient
    Vector              FUb(dim);                                     // BC
    Vector              FPx(dim); // pressure gradient

    Vector              force_at_mid(dim);
    Vector              Res(dim);                                     // residue
    Tensor              dResdu(dim,dim);                              // residue derivative
    Tensor const        I(Tensor::Identity(dim,dim));
    Vector              vec(dim);     // temp
    Tensor              Ten(dim,dim); // temp

    VectorXi            mapU_c(n_dofs_u_per_cell);
    VectorXi            mapU_r(n_dofs_u_per_corner);
    VectorXi            mapP_c(n_dofs_p_per_cell);
    VectorXi            mapP_r(n_dofs_p_per_corner);
    // mesh velocity
    VectorXi            mapM_c(dim*nodes_per_cell);
    VectorXi            mapM_f(dim*nodes_per_facet);
    VectorXi            mapM_r(dim*nodes_per_corner);

    MatrixXd            Prj(n_dofs_u_per_cell,n_dofs_u_per_cell); // projector matrix
    //VectorXi            cell_nodes(nodes_per_cell);


    const int tid = omp_get_thread_num();
    const int nthreads = omp_get_num_threads();

    cell_iterator cell = mesh->cellBegin(tid,nthreads);
    cell_iterator cell_end = mesh->cellEnd(tid,nthreads);

    //cell_iterator cell = mesh->cellBegin();
    //cell_iterator cell_end = mesh->cellEnd();
    for (; cell != cell_end; ++cell)
    {

      tag = cell->getTag();

      // mapeamento do local para o global:
      //
      dof_handler[DH_MESH].getVariable(VAR_M).getCellDofs(mapM_c.data(), &*cell);  //cout << mapM_c << endl;  //unk. global ID's
      dof_handler[DH_UNKS].getVariable(VAR_U).getCellDofs(mapU_c.data(), &*cell);  //cout << mapU_c << endl;
      dof_handler[DH_UNKS].getVariable(VAR_P).getCellDofs(mapP_c.data(), &*cell);  //cout << mapP_c << endl;

      if ((is_bdf2 && time_step > 0) || (is_bdf3 && time_step > 1))
        VecGetValues(Vec_v_1, mapM_c.size(), mapM_c.data(), v_coefs_c_mid.data());
      else
        VecGetValues(Vec_v_mid, mapM_c.size(), mapM_c.data(), v_coefs_c_mid.data());  //size of vector mapM_c
      VecGetValues(Vec_x_0,     mapM_c.size(), mapM_c.data(), x_coefs_c_old.data());
      VecGetValues(Vec_x_1,     mapM_c.size(), mapM_c.data(), x_coefs_c_new.data());
      VecGetValues(Vec_up_0,    mapU_c.size(), mapU_c.data(), u_coefs_c_old.data());
      VecGetValues(Vec_up_k,    mapU_c.size(), mapU_c.data(), u_coefs_c_new.data());
      VecGetValues(Vec_up_k,    mapP_c.size(), mapP_c.data(), p_coefs_c_new.data());
      VecGetValues(Vec_up_0,    mapP_c.size(), mapP_c.data(), p_coefs_c_old.data());

      VecGetValues(Vec_dup,    mapU_c.size(), mapU_c.data(), du_coefs_c_old.data()); // bdf2,bdf3
      if (is_bdf3)
        VecGetValues(Vec_dup_0,  mapU_c.size(), mapU_c.data(), du_coefs_c_vold.data()); // bdf3

      // get nodal coordinates of the old and new cell
      mesh->getCellNodesId(&*cell, cell_nodes.data());
      //mesh->getNodesCoords(cell_nodes.begin(), cell_nodes.end(), x_coefs_c.data());
      //x_coefs_c_trans = x_coefs_c_mid_trans;

      v_coefs_c_mid_trans = v_coefs_c_mid.transpose();
      x_coefs_c_old_trans = x_coefs_c_old.transpose();
      x_coefs_c_new_trans = x_coefs_c_new.transpose();
      u_coefs_c_old_trans = u_coefs_c_old.transpose();
      u_coefs_c_new_trans = u_coefs_c_new.transpose();
  
      du_coefs_c_old_trans = du_coefs_c_old.transpose(); // bdf2
      if (is_bdf3)
        du_coefs_c_vold_trans = du_coefs_c_vold.transpose(); // bdf3

      u_coefs_c_mid_trans = utheta*u_coefs_c_new_trans + (1.-utheta)*u_coefs_c_old_trans;
      x_coefs_c_mid_trans = utheta*x_coefs_c_new_trans + (1.-utheta)*x_coefs_c_old_trans;
      p_coefs_c_mid       = utheta*p_coefs_c_new       + (1.-utheta)*p_coefs_c_old;

      // test: erase me. Interpolated mesh velocity
      if (false)
      {
        //if (iter==0 && (&*cell)==mesh->getCellPtr(0))
        //{
        //  printf("COEFSSSSSSSSS: \n");
        //  cout << v_coefs_c_mid_trans << endl << endl;
        //}
        for (int j = 0; j < (int)v_coefs_c_mid_trans.cols(); ++j)
        {
          for (int c = 0; c < dim; ++c)
            Xqp(c) = x_coefs_c_mid_trans(c,j);
          Vqp = v_exact(Xqp, current_time+dt/2., tag);
          for (int c = 0; c < dim; ++c)
            v_coefs_c_mid_trans(c,j) = Vqp(c);
        }
        //if (iter==0 && (&*cell)==mesh->getCellPtr(0))
        //{
        //  printf("EXACTTTTTTT: \n");
        //  cout << v_coefs_c_mid_trans << endl << endl;
        //}
        //if (iter==0 && (&*cell)==mesh->getCellPtr(0))
        //{
        //  printf("COORD: \n");
        //  cout << x_coefs_c_old_trans << endl << endl;
        //}        
        //for (int i = 0; i < v_coefs_c_mid_trans.rows(); ++i)
        //  for (int j = 0; j < v_coefs_c_mid_trans.cols(); ++j)
        //    v_coefs_c_mid_trans(i,j) += dt*dt;
      }

      visc = muu(tag);
      rho  = pho(Xqp,tag);
      Aloc.setZero();
      Gloc.setZero();
      Dloc.setZero();
      FUloc.setZero();
      FPloc.setZero();
      Eloc.setZero();
      double ddt_factor;
      if (is_bdf2 && time_step > 0)
        ddt_factor = 1.5;
      else
      if (is_bdf3 && time_step > 1)
        ddt_factor = 11./6.;
      else
        ddt_factor = 1.;


      if (behaviors & BH_bble_condens_PnPn) // reset matrices
      {
        iBbb.setZero();
        Bnb.setZero();
        Gbp.setZero();
        FUb.setZero();
        Bbn.setZero();
        Dpb.setZero();
      }

      if(behaviors & BH_GLS)
      {
        cell_volume = 0;
        for (int qp = 0; qp < n_qpts_cell; ++qp) {
          F_c_mid = x_coefs_c_mid_trans * dLqsi_c[qp];
          J_mid = determinant(F_c_mid,dim);
          cell_volume += J_mid * quadr_cell->weight(qp);
        }  //cout << J_mid << " " << cell_volume << endl;

        hk2 = cell_volume / pi; // element size
        double const uconv = (u_coefs_c_old - v_coefs_c_mid).lpNorm<Infinity>();

        tauk = 4.*visc/hk2 + 2.*rho*uconv/sqrt(hk2);
        tauk = 1./tauk;
        if (dim==3)
          tauk *= 0.1;

        delk = 4.*visc + 2.*rho*uconv*sqrt(hk2);
        //delk = 0;

        Eloc.setZero();
        Cloc.setZero();
      }
      if (behaviors & BH_bble_condens_CR)
      {
        bble_integ = 0;
        Gnx.setZero();
        iBbb.setZero();
        Bnb.setZero();
        FUb.setZero();
        FPx.setZero();
        Bbn.setZero();

        cell_volume = 0;
        Xc.setZero();
        for (int qp = 0; qp < n_qpts_cell; ++qp) {
          F_c_mid = x_coefs_c_mid_trans * dLqsi_c[qp];
          J_mid = determinant(F_c_mid,dim);
          Xqp  = x_coefs_c_mid_trans * qsi_c[qp];
          cell_volume += J_mid * quadr_cell->weight(qp);
          Xc += J_mid * quadr_cell->weight(qp) * Xqp;
        }
        Xc /= cell_volume;
      }


      // Quadrature
      for (int qp = 0; qp < n_qpts_cell; ++qp)
      {
        //F_c_new = x_coefs_c_new_trans * dLqsi_c[qp];
        //inverseAndDet(F_c_new,dim,invF_c_new,J_new);
        //invFT_c_new= invF_c_new.transpose();

        //F_c_old = x_coefs_c_old_trans * dLqsi_c[qp];
        //inverseAndDet(F_c_old,dim,invF_c_old,J_old);
        //invFT_c_old= invF_c_old.transpose();

        F_c_mid = x_coefs_c_mid_trans * dLqsi_c[qp];  // (dim x nodes_per_cell) (nodes_per_cell x dim)
        F_c_old = x_coefs_c_old_trans * dLqsi_c[qp];
        F_c_new = x_coefs_c_new_trans * dLqsi_c[qp];
        inverseAndDet(F_c_mid,dim,invF_c_mid,J_mid);
        inverseAndDet(F_c_old,dim,invF_c_old,J_old);
        inverseAndDet(F_c_new,dim,invF_c_new,J_new);
        invFT_c_mid= invF_c_mid.transpose();
        invFT_c_old= invF_c_old.transpose();
        invFT_c_new= invF_c_new.transpose();

        dxphi_c_new = dLphi_c[qp] * invF_c_new;
        dxphi_c     = dLphi_c[qp] * invF_c_mid;
        dxpsi_c_new = dLpsi_c[qp] * invF_c_new;
        dxpsi_c     = dLpsi_c[qp] * invF_c_mid;
        dxqsi_c     = dLqsi_c[qp] * invF_c_mid;

        dxP_new  = dxpsi_c.transpose() * p_coefs_c_new;
        dxU      = u_coefs_c_mid_trans * dxphi_c;       // n+utheta
        dxU_old  = u_coefs_c_old_trans * dLphi_c[qp] * invF_c_old;       // n
        dxU_new  = u_coefs_c_new_trans * dLphi_c[qp] * invF_c_new;       // n+1

        Xqp      = x_coefs_c_mid_trans * qsi_c[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
        Xqp_old  = x_coefs_c_old_trans * qsi_c[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
        Uqp      = u_coefs_c_mid_trans * phi_c[qp]; //n+utheta
        Uqp_new  = u_coefs_c_new_trans * phi_c[qp]; //n+utheta
        Uqp_old  = u_coefs_c_old_trans * phi_c[qp]; //n+utheta
        Pqp_new  = p_coefs_c_new.dot(psi_c[qp]);
        Pqp      = p_coefs_c_mid.dot(psi_c[qp]);
        Vqp      = v_coefs_c_mid_trans * qsi_c[qp];
        //Vqp = v_exact(Xqp_old, current_time, tag);
        //Vqp = v_exact(Xqp, current_time+dt/2., tag);
        Uconv_qp = Uqp - Vqp;
        //Uconv_qp = Uqp_old;
        dUdt     = (Uqp_new-Uqp_old)/dt;
    
        if (is_bdf2 && time_step > 0)
        {
          dUqp_old  = du_coefs_c_old_trans * phi_c[qp]; //n+utheta
          dUdt = 1.5*dUdt - .5*dUqp_old;
        }
        else
        if (is_bdf3 && time_step > 1)
        {
          dUqp_old   = du_coefs_c_old_trans  * phi_c[qp];
          dUqp_vold  = du_coefs_c_vold_trans * phi_c[qp];
          dUdt = 11./6.*dUdt - 7./6.*dUqp_old + 1./3.*dUqp_vold;
        }


        force_at_mid = force(Xqp,current_time+utheta*dt,tag);

        weight = quadr_cell->weight(qp);
        JxW_mid = J_mid*weight;
        JxW_old = J_old*weight;
        JxW_new = J_new*weight;
//cout << "J_mid = " << J_mid << endl;
        //~ if (mesh->getCellId(&*cell) == 0)
        //~ {
          //~ printf("cHEcKKKKKKKKKKK!!\n");
          //~ cout << "x coefs mid:" << endl;
          //~ cout << x_coefs_c_mid_trans.transpose() << endl;
        //~ }
        if (J_mid < 1.e-14)
        {
          FEP_PRAGMA_OMP(critical)
          //if (tid==0)
          {
            printf("in formCellFunction:\n");
            std::cout << "erro: jacobiana da integral não invertível: ";
            std::cout << "J_mid = " << J_mid << endl;
            cout << "trans(f) matrix:\n" << F_c_mid << endl;
            cout << "x coefs mid:" << endl;
            cout << x_coefs_c_mid_trans.transpose() << endl;
            cout << "-----" << endl;
            cout << "cell id: " << mesh->getCellId(&*cell) << endl;
            cout << "cell Contig id: " << mesh->getCellContigId( mesh->getCellId(&*cell) ) << endl;
            cout << "cell nodes:\n" << cell_nodes.transpose() << endl;
            cout << "mapM :\n" << mapM_c.transpose() << endl;
            throw;
          }
        }

        for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
        {
          for (int c = 0; c < dim; ++c)
          {
            FUloc(i*dim + c) += JxW_mid*
                    ( rho*(unsteady*dUdt(c) + has_convec*Uconv_qp.dot(dxU.row(c)))*phi_c[qp][i] + // aceleração
                      visc*dxphi_c.row(i).dot(dxU.row(c) + dxU.col(c).transpose())  //rigidez
                    ) -
                    JxW_mid*force_at_mid(c)*phi_c[qp][i] -// força +
                    //JxW_new*Pqp_new*dxphi_c_new(i,c);  // pressão
                    //JxW_mid*Pqp*dxphi_c(i,c);  // pressão
                    JxW_mid*Pqp_new*dxphi_c(i,c);  // pressão

            for (int j = 0; j < n_dofs_u_per_cell/dim; ++j)
            {
              for (int d = 0; d < dim; ++d)
              {
                delta_cd = c==d;
                Aloc(i*dim + c, j*dim + d) += JxW_mid*
                                              ( has_convec*phi_c[qp][i]*utheta *rho*( delta_cd*Uconv_qp.dot(dxphi_c.row(j))  +  dxU(c,d)*phi_c[qp][j] )   // advecção
                                             + ddt_factor* unsteady* delta_cd*rho*phi_c[qp][i]*phi_c[qp][j]/dt     // time derivative
                                              + utheta*visc*( delta_cd * dxphi_c.row(i).dot(dxphi_c.row(j)) + dxphi_c(i,d)*dxphi_c(j,c))   ); // rigidez

              }
            }
            for (int j = 0; j < n_dofs_p_per_cell; ++j)
            {
              Gloc(i*dim + c,j) -= JxW_mid * psi_c[qp][j]* dxphi_c(i,c);
              //Gloc(i*dim + c,j) -= utheta*JxW_mid * psi_c[qp][j]* dxphi_c(i,c);
              //Gloc(i*dim + c,j) -= JxW_new * psi_c[qp][j]* dxphi_c_new(i,c);
              //Dloc(j, i*dim + c) -= JxW_new * psi_c[qp][j]*  dxphi_c_new(i,c);
              Dloc(j, i*dim + c) -= utheta*JxW_mid * psi_c[qp][j]*  dxphi_c(i,c);
            }

          }

          //FUloc.segment(i*dim,dim) += JxW_mid*
          //                        (   phi_c[qp][i]*rho* (  dUdt + has_convec* dxU*Uconv_qp  ) // aceleração
          //                          + visc* ( dxU + dxU.transpose() )*dxphi_c.row(i).transpose()       //rigidez
          //                          - Pqp_new*dxphi_c.row(i).transpose() // pressão
          //                          - phi_c[qp][i]* force_at_mid   ); // força

        }
        for (int i = 0; i < n_dofs_p_per_cell; ++i)
          FPloc(i) -= JxW_mid* dxU.trace()*psi_c[qp][i];
          //FPloc(i) -= JxW_new* dxU_new.trace() *psi_c[qp][i];

        // ----------------
        //
        //  STABILIZATION
        //
        //  ----------------
        if (behaviors & (BH_bble_condens_PnPn | BH_bble_condens_CR))
        {
          dxbble = invFT_c_mid * dLbble[qp];

          for (int c = 0; c < dim; c++)
          {
            for (int j = 0; j < n_dofs_u_per_cell/dim; j++)
            {
              for (int d = 0; d < dim; d++)
              {
                delta_cd = c==d;
    
                if (compact_bubble)
                {
                  Bbn(c, j*dim + d) += JxW_mid*
                                       ( utheta*visc*(delta_cd * dxphi_c.row(j).dot(dxbble) + dxphi_c(j,c)*dxbble(d)) ); // rigidez

                  Bnb(j*dim + d, c) += JxW_mid*
                                       ( utheta*visc*(delta_cd * dxphi_c.row(j).dot(dxbble) + dxphi_c(j,c)*dxbble(d)) ); // rigidez  
                }
                else
                {
                  Bbn(c, j*dim + d) += JxW_mid*
                                       ( has_convec*bble[qp]*utheta *rho*( delta_cd*Uconv_qp.dot(dxphi_c.row(j)) + dxU(c,d)*phi_c[qp][j] ) // convective
                                       + unsteady*delta_cd*rho*bble[qp]*phi_c[qp][j]/dt // time derivative
                                       + utheta*visc*(delta_cd * dxphi_c.row(j).dot(dxbble) + dxphi_c(j,c)*dxbble(d)) ); // rigidez

                  Bnb(j*dim + d, c) += JxW_mid*
                                       ( has_convec*phi_c[qp][j]*utheta *rho*( delta_cd*Uconv_qp.dot(dxbble) ) // convective
                                       + delta_cd*rho*phi_c[qp][j]*bble[qp]/dt * unsteady // time derivative
                                       + utheta*visc*(delta_cd * dxphi_c.row(j).dot(dxbble) + dxphi_c(j,c)*dxbble(d)) ); // rigidez
                }
              }
            }
            if (behaviors & BH_bble_condens_PnPn)
              for (int j = 0; j < n_dofs_p_per_cell; ++j)
                Gbp(c, j) -= JxW_mid*psi_c[qp][j]*dxbble(c);
          }

          for (int c = 0; c < dim; c++)
          {
            for (int d = 0; d < dim; d++)
            {
              delta_cd = c==d;
              
              if (compact_bubble)
              {
                iBbb(c, d) += JxW_mid*
                              ( utheta*visc*(delta_cd* dxbble.dot(dxbble) + dxbble(d)*dxbble(c)) ); // rigidez  
              }
              else
              {
                iBbb(c, d) += JxW_mid*
                              ( has_convec*bble[qp]*utheta *rho*( delta_cd*Uconv_qp.dot(dxbble) ) // convective
                              + delta_cd*rho*bble[qp]*bble[qp]/dt * unsteady // time derivative
                              + utheta*visc*(delta_cd* dxbble.dot(dxbble) + dxbble(d)*dxbble(c)) ); // rigidez  
              }
              
              
            }
            if (compact_bubble)
            {
              FUb(c) += JxW_mid*
                        ( visc*dxbble.dot(dxU.row(c) + dxU.col(c).transpose()) - //rigidez
                          Pqp_new*dxbble(c) - // pressão
                          force_at_mid(c)*bble[qp] ); // força  
            }
            else
            {
              FUb(c) += JxW_mid*
                        ( bble[qp]*rho*(dUdt(c)*unsteady + has_convec*Uconv_qp.dot(dxU.row(c))) + // time derivative + convective
                          visc*dxbble.dot(dxU.row(c) + dxU.col(c).transpose()) - //rigidez
                          Pqp_new*dxbble(c) - // pressão
                          force_at_mid(c)*bble[qp] ); // força  
            }
            
          }
        }
        else
        if(behaviors & BH_GLS)
        {
          Res = rho*( dUdt * unsteady + has_convec*dxU*Uconv_qp) + dxP_new - force_at_mid;

          for (int j = 0; j < n_dofs_u_per_cell/dim; ++j)
          {
            dResdu = unsteady*(ddt_factor*rho*phi_c[qp][j]/dt)*I + has_convec*rho*utheta*( phi_c[qp][j]*dxU + Uconv_qp.dot(dxphi_c.row(j))*I );

            for (int i = 0; i < n_dofs_p_per_cell; ++i)
            {
              vec = dxpsi_c.row(i).transpose();
              vec = dResdu.transpose()*vec;
              vec = -JxW_mid*tauk* vec;
              for (int d = 0; d < dim; d++)
                Dloc(i, j*dim + d) += vec(d);

              // atençao nos indices
              vec = JxW_mid*tauk* has_convec*rho*Uconv_qp.dot(dxphi_c.row(j))* dxpsi_c.row(i).transpose();
              for (int d = 0; d < dim; d++)
                Cloc(j*dim + d,i) += vec(d);
            }

            for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
            {
              // supg term
              Ten = JxW_mid*tauk* has_convec*( utheta*rho*phi_c[qp][j]*Res*dxphi_c.row(i) + rho*Uconv_qp.dot(dxphi_c.row(i))*dResdu );
              // divergence term
              Ten+= JxW_mid*delk*utheta*dxphi_c.row(i).transpose()*dxphi_c.row(j);

              for (int c = 0; c < dim; ++c)
                for (int d = 0; d < dim; ++d)
                  Aloc(i*dim + c, j*dim + d) += Ten(c,d);
            }
          }

          for (int i = 0; i < n_dofs_p_per_cell; ++i)
            for (int j = 0; j < n_dofs_p_per_cell; ++j)
              Eloc(i,j) -= tauk*JxW_mid * dxphi_c.row(i).dot(dxphi_c.row(j));


          for (int i = 0; i < n_dofs_u_per_cell/dim; i++)
          {
            vec = JxW_mid*( has_convec*tauk* rho* Uconv_qp.dot(dxphi_c.row(i)) * Res + delk*dxU.trace()*dxphi_c.row(i).transpose() );

            for (int c = 0; c < dim; c++)
              FUloc(i*dim + c) += vec(c);

          }
          for (int i = 0; i < n_dofs_p_per_cell; ++i)
            FPloc(i) -= JxW_mid *tauk* dxpsi_c.row(i).dot(Res);
            //FPloc(i) -= JxW_mid *tauk* dxpsi_c.row(i).dot(dxP_new - force_at_mid); // somente laplaciano da pressao
        }

        if (behaviors & BH_bble_condens_CR)
        {
          bble_integ += JxW_mid*bble[qp];

          for (int c = 0; c < dim; ++c)
          {
            for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
              for (int j = 0; j < dim; ++j) // pressure gradient
                Gnx(i*dim + c,j) -= JxW_mid* (Xqp(j) - Xc(j))*dxphi_c(i,c);

            FPx(c) -= JxW_mid* dxU.trace()*(Xqp(c) - Xc(c));
          }
        }


      } // fim quadratura

      //Dloc += utheta*Gloc.transpose();

      //
      // stabilization
      //
      if ((behaviors & BH_bble_condens_PnPn) && !compact_bubble)
      {
        //iBbb = iBbb.inverse().eval();
        invert(iBbb,dim);

        FUloc = FUloc - Bnb*iBbb*FUb;
        FPloc = FPloc - utheta*Gbp.transpose()*iBbb*FUb;

        Dpb = utheta*Gbp.transpose();

        // correções com os coeficientes da bolha

        Ubqp = -utheta*iBbb*FUb; // U bolha no tempo n+utheta

        for (int qp = 0; qp < n_qpts_cell; ++qp)
        {
          F_c_mid = x_coefs_c_mid_trans * dLqsi_c[qp];
          inverseAndDet(F_c_mid, dim, invF_c_mid,J_mid);
          invFT_c_mid= invF_c_mid.transpose();

          Uqp = u_coefs_c_mid_trans * phi_c[qp]; //n+utheta
          dxbble = invFT_c_mid * dLbble[qp];
          dxUb = Ubqp*dxbble.transpose();

          weight = quadr_cell->weight(qp);
          JxW_mid = J_mid*weight;

          for (int j = 0; j < n_dofs_u_per_cell/dim; ++j)
          {
            for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
            {
              Ten = has_convec*JxW_mid*rho*utheta*phi_c[qp][i]* phi_c[qp][j] * dxUb; // advecção

              for (int c = 0; c < dim; ++c)
                for (int d = 0; d < dim; ++d)
                  Aloc(i*dim + c, j*dim + d) += Ten(c,d);
            }
            Ten = has_convec*JxW_mid* rho*utheta* bble[qp] * phi_c[qp][j] *dxUb; // advecção
            for (int c = 0; c < dim; ++c)
              for (int d = 0; d < dim; ++d)
                Bbn(c, j*dim + d) += Ten(c,d);
          }
        } // fim quadratura 2 vez

        Aloc -= Bnb*iBbb*Bbn;
        Gloc -= Bnb*iBbb*Gbp;
        Dloc -= Dpb*iBbb*Bbn;
        Eloc = -Dpb*iBbb*Gbp;

      }
      if(behaviors & BH_GLS)
      {
        Gloc += Cloc;
      }
      if(behaviors & BH_bble_condens_CR)
      {
        Ubqp.setZero();
        //for (int i = 0; i < Gnx.cols(); ++i)
        // for (int j = 0; j < Gnx.rows(); ++j)
        // Ubqp(i) += Gnx(j,i)*u_coefs_c_new(j);
        //Ubqp /= -bble_integ;
        //Ubqp *= utheta;

        for (int c = 0; c < dim; ++c)
          for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
            for (int j = 0; j < dim; ++j) // pressure gradient
              //Gnx(i*dim + c,j) -= JxW_mid* (Xqp(j) - Xc(j))*dxphi_c(i,c);
            Ubqp(j) += Gnx(i*dim + c,j) * u_coefs_c_mid_trans(c,i);

        Ubqp /= -bble_integ;
        Ubqp *= utheta;


        //Ubqp = -Gnx.transpose()*u_coefs_c_new; // U bolha no tempo n+utheta

        for (int qp = 0; qp < n_qpts_cell; ++qp)
        {
          F_c_mid = x_coefs_c_mid_trans * dLqsi_c[qp];
          inverseAndDet(F_c_mid, dim, invF_c_mid,J_mid);
          invFT_c_mid= invF_c_mid.transpose();

          Uqp = u_coefs_c_mid_trans * phi_c[qp]; //n+utheta
          dxbble = invFT_c_mid * dLbble[qp];
          dxUb = Ubqp*dxbble.transpose();

          weight = quadr_cell->weight(qp);
          JxW_mid = J_mid*weight;

          for (int j = 0; j < n_dofs_u_per_cell/dim; ++j)
          {
            for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
            {
              Ten = has_convec*JxW_mid*rho*utheta*phi_c[qp][i]* phi_c[qp][j] * dxUb; // advecção

              for (int c = 0; c < dim; ++c)
                for (int d = 0; d < dim; ++d)
                  Aloc(i*dim + c, j*dim + d) += Ten(c,d);
            }
            Ten = has_convec*JxW_mid* rho*utheta* bble[qp] * phi_c[qp][j] *dxUb; // advecção
            for (int c = 0; c < dim; ++c)
              for (int d = 0; d < dim; ++d)
                Bbn(c, j*dim + d) += Ten(c,d);
          }
        } // fim quadratura 2 vez

        double const a = 1./(bble_integ*bble_integ);
        double const b = 1./bble_integ;
        Aloc += utheta*a*Gnx*iBbb*Gnx.transpose() - utheta*b*Bnb*Gnx.transpose() - b*Gnx*Bbn;

        FUloc += a*Gnx*iBbb*FPx - b*Bnb*FPx - b*Gnx*FUb;
      }


//cout << "\n" << FUloc << endl; cout << "\n" << Aloc << endl; cout << "\n" << Gloc << endl; cout << "\n" << Dloc << endl;
      // Projection - to force non-penetrarion bc
      mesh->getCellNodesId(&*cell, cell_nodes.data());
      getProjectorMatrix(Prj, nodes_per_cell, cell_nodes.data(), Vec_x_1, current_time+dt, *this);
//cout << Prj << endl;
      FUloc = Prj*FUloc;
      Aloc = Prj*Aloc*Prj;
      Gloc = Prj*Gloc;
      Dloc = Dloc*Prj;
//cout << "\n" << FUloc << endl; cout << "\n" << Aloc << endl; cout << "\n" << Gloc << endl; cout << "\n" << Dloc << endl;
      if (force_pressure)
      {
        for (int i = 0; i < mapP_c.size(); ++i)
        {
          if (mapP_c(i) == null_space_press_dof)
          {
            Gloc.col(i).setZero();
            Dloc.row(i).setZero();
            FPloc(i) = 0;
            Eloc.col(i).setZero();
            Eloc.row(i).setZero();
            break;
          }
        }
      }

#ifdef FEP_HAS_OPENMP
      FEP_PRAGMA_OMP(critical)
#endif
      {
        VecSetValues(Vec_fun, mapU_c.size(), mapU_c.data(), FUloc.data(), ADD_VALUES);
        VecSetValues(Vec_fun, mapP_c.size(), mapP_c.data(), FPloc.data(), ADD_VALUES);
        MatSetValues(*JJ, mapU_c.size(), mapU_c.data(), mapU_c.size(), mapU_c.data(), Aloc.data(),  ADD_VALUES);
        MatSetValues(*JJ, mapU_c.size(), mapU_c.data(), mapP_c.size(), mapP_c.data(), Gloc.data(),  ADD_VALUES);
        MatSetValues(*JJ, mapP_c.size(), mapP_c.data(), mapU_c.size(), mapU_c.data(), Dloc.data(),  ADD_VALUES);
        MatSetValues(*JJ, mapP_c.size(), mapP_c.data(), mapP_c.size(), mapP_c.data(), Eloc.data(),  ADD_VALUES);
      }
    }  //end for cell


  }  //end LOOP NAS CÉLULAS parallel

  // LOOP NAS FACES DO CONTORNO (Neumann)
  //~ FEP_PRAGMA_OMP(parallel default(none) shared(Vec_up_k,Vec_fun,cout))
  {
    int                 tag;
    bool                is_neumann;
    bool                is_surface;
    bool                is_solid;
    //MatrixXd           u_coefs_c_new(n_dofs_u_per_facet/dim, dim);
    //VectorXd           p_coefs_f(n_dofs_p_per_facet);
    MatrixXd            u_coefs_f_mid_trans(dim, n_dofs_u_per_facet/dim);  // n+utheta
    MatrixXd            u_coefs_f_old(n_dofs_u_per_facet/dim, dim);        // n
    MatrixXd            u_coefs_f_new(n_dofs_u_per_facet/dim, dim);        // n+1
    MatrixXd            u_coefs_f_old_trans(dim,n_dofs_u_per_facet/dim);   // n
    MatrixXd            u_coefs_f_new_trans(dim,n_dofs_u_per_facet/dim);   // n+1

    MatrixXd            x_coefs_f_mid_trans(dim, n_dofs_v_per_facet/dim); // n+utheta
    MatrixXd            x_coefs_f_new(n_dofs_v_per_facet/dim, dim);       // n+1
    MatrixXd            x_coefs_f_new_trans(dim, n_dofs_v_per_facet/dim); // n+1
    MatrixXd            x_coefs_f_old(n_dofs_v_per_facet/dim, dim);       // n
    MatrixXd            x_coefs_f_old_trans(dim, n_dofs_v_per_facet/dim); // n

    MatrixXd            noi_coefs_f_new(n_dofs_v_per_facet/dim, dim);  // normal interpolada em n+1
    MatrixXd            noi_coefs_f_new_trans(dim, n_dofs_v_per_facet/dim);  // normal interpolada em n+1

    Tensor              F_f_mid(dim,dim-1);       // n+utheta
    Tensor              invF_f_mid(dim-1,dim);    // n+utheta
    Tensor              fff_f_mid(dim-1,dim-1);   // n+utheta; fff = first fundamental form
    //Tensor              invFT_c_mid(dim,dim);   // n+utheta

    MatrixXd            Aloc_f(n_dofs_u_per_facet, n_dofs_u_per_facet);
    VectorXd            FUloc(n_dofs_u_per_facet);

    MatrixXd            tmp(n_dofs_u_per_facet,n_dofs_u_per_facet);

    VectorXi            mapU_f(n_dofs_u_per_facet);
    VectorXi            mapP_f(n_dofs_p_per_facet);
    VectorXi            mapM_f(dim*nodes_per_facet);

    MatrixXd            Prj(n_dofs_u_per_facet,n_dofs_u_per_facet);

    MatrixXd            dxphi_f(n_dofs_u_per_facet/dim, dim);
    Tensor              dxU_f(dim,dim);   // grad u
    Vector              Xqp(dim);
    Vector              Xqp2(dim);
    Vector              Xqp_new(dim);
    Vector              Xqp_old(dim);
    Vector              Uqp(dim);
    Vector              Uqp_new(dim);
    Vector              Uqp_old(dim);
    //VectorXd          FUloc(n_dofs_u_per_facet);
    VectorXi            facet_nodes(nodes_per_facet);
    Vector              normal(dim);
    Vector              noi(dim); // normal interpolada
    Vector              some_vec(dim);
    double              J_mid=0,JxW_mid;
    double              weight=0;
    //double              visc;
    //double              rho;
    Vector              Uqp_solid(dim);

    Vector              traction_(dim);

    //~ const int tid = omp_get_thread_num();
    //~ const int nthreads = omp_get_num_threads();
//~
    //~ facet_iterator facet = mesh->facetBegin(tid,nthreads);
    //~ facet_iterator facet_end = mesh->facetEnd(tid,nthreads);

    // LOOP NAS FACES DO CONTORNO
    facet_iterator facet = mesh->facetBegin();
    facet_iterator facet_end = mesh->facetEnd();  // the next if controls the for that follows
    if (neumann_tags.size() != 0 || interface_tags.size() != 0 || solid_tags.size() != 0)
    for (; facet != facet_end; ++facet)
    {
      tag = facet->getTag();
      is_neumann = is_in(tag, neumann_tags);
      is_surface = is_in(tag, interface_tags);
      is_solid   = is_in(tag, solid_tags);

      //if ((!is_neumann))
      if ((!is_neumann) && (!is_surface) && (!is_solid))
      //PetscFunctionReturn(0);
        continue;

      // mapeamento do local para o global:
      //
      dof_handler[DH_UNKS].getVariable(VAR_U).getFacetDofs(mapU_f.data(), &*facet);  cout << mapU_f << endl << endl;  //unk. global ID's
      dof_handler[DH_UNKS].getVariable(VAR_P).getFacetDofs(mapP_f.data(), &*facet);  //cout << mapP_f << endl;
      dof_handler[DH_MESH].getVariable(VAR_M).getFacetDofs(mapM_f.data(), &*facet);  //cout << mapM_f << endl;

      VecGetValues(Vec_normal,  mapM_f.size(), mapM_f.data(), noi_coefs_f_new.data());
      VecGetValues(Vec_x_0,     mapM_f.size(), mapM_f.data(), x_coefs_f_old.data());
      VecGetValues(Vec_x_1,     mapM_f.size(), mapM_f.data(), x_coefs_f_new.data());
      VecGetValues(Vec_up_0,    mapU_f.size(), mapU_f.data(), u_coefs_f_old.data());
      VecGetValues(Vec_up_k ,   mapU_f.size(), mapU_f.data(), u_coefs_f_new.data());

      // get nodal coordinates of the old and new cell
      mesh->getFacetNodesId(&*facet, facet_nodes.data());

      x_coefs_f_old_trans = x_coefs_f_old.transpose();
      x_coefs_f_new_trans = x_coefs_f_new.transpose();
      u_coefs_f_old_trans = u_coefs_f_old.transpose();
      u_coefs_f_new_trans = u_coefs_f_new.transpose();
      noi_coefs_f_new_trans = noi_coefs_f_new.transpose();

      u_coefs_f_mid_trans = utheta*u_coefs_f_new_trans + (1.-utheta)*u_coefs_f_old_trans;
      x_coefs_f_mid_trans = utheta*x_coefs_f_new_trans + (1.-utheta)*x_coefs_f_old_trans;

      FUloc.setZero();
      Aloc_f.setZero();

      //visc = muu(tag);
      //rho  = pho(Xqp,tag);

      //noi_coefs_f_new_trans = x_coefs_f_mid_trans;


      for (int qp = 0; qp < n_qpts_facet; ++qp)
      {

        F_f_mid   = x_coefs_f_mid_trans * dLqsi_f[qp];

        if (dim==2)
        {
          normal(0) = +F_f_mid(1,0);
          normal(1) = -F_f_mid(0,0);
          normal.normalize();
        }
        else
        {
          normal = cross(F_f_mid.col(0), F_f_mid.col(1));
          normal.normalize();
        }

        fff_f_mid.resize(dim-1,dim-1);
        fff_f_mid  = F_f_mid.transpose()*F_f_mid;
        J_mid     = sqrt(fff_f_mid.determinant());
        invF_f_mid = fff_f_mid.inverse()*F_f_mid.transpose();

        weight  = quadr_facet->weight(qp);
        JxW_mid = J_mid*weight;
        Xqp     = x_coefs_f_mid_trans * qsi_f[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
        dxphi_f = dLphi_f[qp] * invF_f_mid;
        dxU_f   = u_coefs_f_mid_trans * dxphi_f; // n+utheta
        Uqp     = u_coefs_f_mid_trans * phi_f[qp];
        noi     = noi_coefs_f_new_trans * qsi_f[qp];

        if (is_neumann)
        {
          //Vector no(Xqp);
          //no.normalize();
          //traction_ = utheta*(traction(Xqp,current_time+dt,tag)) + (1.-utheta)*traction(Xqp,current_time,tag);
          traction_ = traction(Xqp, normal, current_time + dt*utheta,tag);
          //traction_ = (traction(Xqp,current_time,tag) +4.*traction(Xqp,current_time+dt/2.,tag) + traction(Xqp,current_time+dt,tag))/6.;

          for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
          {
            for (int c = 0; c < dim; ++c)
            {
              FUloc(i*dim + c) -= JxW_mid * traction_(c) * phi_f[qp][i] ; // força
            }
          }
        }

        if (is_surface)
        {
          //Vector no(Xqp);
          //no.normalize();
          for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
          {
            for (int c = 0; c < dim; ++c)
            {
//              FUloc(i*dim + c) += JxW_mid *gama(Xqp,current_time,tag)*(dxphi_f(i,c) + (unsteady*dt) *dxU_f.row(c).dot(dxphi_f.row(i))); // correto
              FUloc(i*dim + c) += JxW_mid *gama(Xqp,current_time,tag)*dxphi_f(i,c); //inicialmente descomentado
              //FUloc(i*dim + c) += JxW_mid *gama(Xqp,current_time,tag)*normal(c)* phi_f[qp][i];
              //for (int d = 0; d < dim; ++d)
              //  FUloc(i*dim + c) += JxW_mid * gama(Xqp,current_time,tag)* ( (c==d?1:0) - noi(c)*noi(d) )* dxphi_f(i,d) ;
              //FUloc(i*dim + c) += JxW_mid * gama(Xqp,current_time,tag)* ( unsteady*dt *dxU_f.row(c).dot(dxphi_f.row(i)));
            }
          }

          if (false) // semi-implicit term //inicialmente false
          {
            for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
              for (int j = 0; j < n_dofs_u_per_facet/dim; ++j)
                for (int c = 0; c < dim; ++c)
                  Aloc_f(i*dim + c, j*dim + c) += utheta*JxW_mid* (unsteady*dt) *gama(Xqp,current_time,tag)*dxphi_f.row(i).dot(dxphi_f.row(j));
          }

        }

        if (is_solid)
        {
          Uqp_solid = solid_veloc(Xqp, current_time+utheta*dt, tag);

          for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
          {
            for (int c = 0; c < dim; ++c)
            {
              FUloc(i*dim + c) += JxW_mid *beta_diss()*(Uqp(c)-Uqp_solid(c))*phi_f[qp][i];
              //FUloc(i*dim + c) += x_coefs_f_old_trans.norm()*beta_diss()*Uqp(c)*phi_f[qp][i];

              for (int j = 0; j < n_dofs_u_per_facet/dim; ++j)
                  Aloc_f(i*dim + c, j*dim + c) += utheta*JxW_mid *beta_diss()*phi_f[qp][j]*phi_f[qp][i];

            }
          }
        }

      } // end quadratura


      // Projection - to force non-penetrarion bc
      mesh->getFacetNodesId(&*facet, facet_nodes.data());
      getProjectorMatrix(Prj, nodes_per_facet, facet_nodes.data(), Vec_x_1, current_time+dt, *this);

      FUloc = Prj*FUloc;
      Aloc_f = Prj*Aloc_f*Prj;

      //~ FEP_PRAGMA_OMP(critical)
      {
        VecSetValues(Vec_fun, mapU_f.size(), mapU_f.data(), FUloc.data(), ADD_VALUES);
        MatSetValues(*JJ, mapU_f.size(), mapU_f.data(), mapU_f.size(), mapU_f.data(), Aloc_f.data(),  ADD_VALUES);
      }

    }  //end for facet


  } // end LOOP NAS FACES DO CONTORNO (Neumann)
  

  // LINHA DE CONTATO
  //FEP_PRAGMA_OMP(parallel shared(Vec_up_k,Vec_fun,cout) default(none))
  {
    // will be useful I hope
    //Real const eps = std::numeric_limits<Real>::epsilon();
    //Real const eps_root = pow(eps,1./3.);
    //double            h;
    //volatile double   hh;

    int              tag;
    bool             is_triple;

    VectorXd         FUloc(n_dofs_u_per_corner);
    MatrixXd         Aloc_r(n_dofs_u_per_corner, n_dofs_u_per_corner);

    VectorXi         mapU_r(n_dofs_u_per_corner);
    VectorXi         mapP_r(n_dofs_p_per_corner);

    MatrixXd         Prj(n_dofs_u_per_corner,n_dofs_u_per_corner);
    VectorXi         corner_nodes(nodes_per_corner);

    bool                gen_error = false;
    //MatrixXd             u_coefs_r_mid(n_dofs_u_per_corner/dim, dim);
    MatrixXd            u_coefs_r_mid_trans(dim, n_dofs_u_per_corner/dim);  // n+utheta
    MatrixXd            u_coefs_r_old(n_dofs_u_per_corner/dim, dim);        // n
    MatrixXd            u_coefs_r_new(n_dofs_u_per_corner/dim, dim);        // n+1
    MatrixXd            x_coefs_r_mid_trans(dim, nodes_per_corner);
    MatrixXd            x_coefs_r_new(nodes_per_corner, dim);
    MatrixXd            x_coefs_r_old(nodes_per_corner, dim);
    Tensor              F_r_mid(dim,dim-2);
    Tensor              invF_r_mid(dim-2,dim);
    MatrixXd            dxphi_r(n_dofs_u_per_corner/dim, dim);
    Tensor              dxU_r(dim,dim);   // grad u
    Vector              Xqp(dim);
    Vector              Uqp(dim);
    //VectorXd          FUloc(n_dofs_u_per_corner);
    //MatrixXd         Aloc_r(n_dofs_u_per_corner, n_dofs_u_per_corner);
    Vector              normal(dim);
    Vector              line_normal(dim);
    Vector              solid_point(dim); // ponto na superfície do sólido .. ele é único
    Vector              point_a(dim); // ponto na linha triplice
    Vector              point_b(dim); // ponto na linha triplice
    Vector              ifacet_normal(dim); // ponto na linha triplice
    double              line_normal_sign = 0; // +1 or -1
    double              J_mid=0, JxW_mid;
    double              weight=0;
    double              gama_mid;
    //double              visc;
    //double              rho;
    int                 iCs[FEPIC_MAX_ICELLS];
    int                 eiCs[FEPIC_MAX_ICELLS];
    int                 *iCs_end;
    int                 *iCs_it;
    Cell                *fluid_cell;

    //VectorXi            mapU_r(n_dofs_u_per_corner);
    //VectorXi            mapP_r(n_dofs_p_per_corner);
    VectorXi            mapM_r(dim*nodes_per_corner);

    //const int tid = omp_get_thread_num();
    //const int nthreads = omp_get_num_threads();
    //const int n_corner_colors = mesh->numCornerColors();


    // LOOP NAS ARESTAS DA LINHA TRIPLICE
    CellElement * corner;

    if (triple_tags.size() != 0)
    for (int _r = 0; _r < n_corners_total; ++_r)
    {
      if (dim==2)
      {
        corner = mesh->getNodePtr(_r);
        if (!mesh->isVertex(corner))
          continue;
      }
      else
        corner = mesh->getCornerPtr(_r);
      if (corner->isDisabled())
        continue;

      tag = corner->getTag();
      is_triple = is_in(tag,triple_tags);
      if (!is_triple)
        continue;

      FUloc.setZero();
      Aloc_r.setZero();

      mesh->getCornerNodesId(&*corner, corner_nodes.data());
     // mesh->getNodesCoords(corner_nodes.begin(), corner_nodes.end(), x_coefs_r_mid.data());

      dof_handler[DH_UNKS].getVariable(VAR_U).getCornerDofs(mapU_r.data(), &*corner);
      dof_handler[DH_UNKS].getVariable(VAR_P).getCornerDofs(mapP_r.data(), &*corner);
      dof_handler[DH_MESH].getVariable(VAR_M).getCornerDofs(mapM_r.data(), &*corner);

      VecGetValues(Vec_x_0,     mapM_r.size(), mapM_r.data(), x_coefs_r_old.data());
      VecGetValues(Vec_x_1,     mapM_r.size(), mapM_r.data(), x_coefs_r_new.data());
      VecGetValues(Vec_up_0,    mapU_r.size(), mapU_r.data(), u_coefs_r_old.data());
      VecGetValues(Vec_up_k,    mapU_r.size(), mapU_r.data(), u_coefs_r_new.data());

      u_coefs_r_mid_trans = utheta*u_coefs_r_new.transpose() + (1.-utheta)*u_coefs_r_old.transpose();
      x_coefs_r_mid_trans = utheta*x_coefs_r_new.transpose() + (1.-utheta)*x_coefs_r_old.transpose();

      //visc = muu(tag);
      //rho  = pho(Xqp,tag);

      if (dim==3)
      {
        // encontrando o ponto da superfície sólido. com ele, é calculado uma normal
        // que corrigi o sinal de line_normal
        iCs_end = mesh->edgeStar(corner, iCs, eiCs);
        if (iCs_end == iCs)
        {
          printf("ERROR!: no icell found\n");
          throw;
        }
        gen_error = true;
        for (iCs_it = iCs; iCs_it != iCs_end ; ++iCs_it)
        {
          fluid_cell = mesh->getCellPtr(*iCs_it);

          for (int kk = 0; kk < mesh->numVerticesPerCell(); ++kk)
          {
            int const nodekk_id = fluid_cell->getNodeId(kk);
            Point const* pp      = mesh->getNodePtr(fluid_cell->getNodeId(kk) );
            const int    tag_aux = pp->getTag();

            if ((is_in(tag_aux, solid_tags) || is_in(tag_aux,feature_tags)) && !is_in(tag_aux, triple_tags))
            {
              if (corner_nodes(0) != nodekk_id && corner_nodes(1) != nodekk_id)
              {
                gen_error = false;
                pp->getCoord(solid_point.data(), dim);
                break;
              }
            }
          }
          if (gen_error==false)
            break;

        }
        if (gen_error)
        {
          printf("ERROR!: solid point not found\n");
          cout << "corner id: " << (mesh->getCellPtr(corner->getIncidCell())->getCornerId(corner->getPosition())) << endl;
          cout << "first icell : " << (*iCs) << endl;
          mesh->getNodePtr(corner_nodes(0))->getCoord(point_a.data(),dim);
          mesh->getNodePtr(corner_nodes(1))->getCoord(point_b.data(),dim);
          cout << "point a = " << point_a[0] << " " << point_a[1] << " " << point_a[2] << "\n";
          cout << "point b = " << point_b[0] << " " << point_b[1] << " " << point_b[2] << "\n";
          cout << "number of cells = " << (iCs_end-iCs) << endl;
          throw;
        }


        mesh->getNodePtr(corner_nodes(0))->getCoord(point_a.data(),dim);
        mesh->getNodePtr(corner_nodes(1))->getCoord(point_b.data(),dim);

        // se (a-c) cross (b-c) dot solid_normal > 0, então line_normal_sign = 1, se não, =-1
        Xqp = (point_a+point_b)/2.;
        point_a -= solid_point;
        point_b -= solid_point;
        normal  = solid_normal(Xqp, current_time, tag);
        line_normal_sign = point_a(0)*point_b(1)*normal(2) + point_a(1)*point_b(2)*normal(0) + point_a(2)*point_b(0)*normal(1)
                          -point_a(0)*point_b(2)*normal(1) - point_a(1)*point_b(0)*normal(2) - point_a(2)*point_b(1)*normal(0);
        if (line_normal_sign>0)
          line_normal_sign = 1;
        else
          line_normal_sign = -1;


      }
      else // dim==2
      {
        Point * point = mesh->getNodePtr(corner_nodes[0]);
        Point * sol_point = NULL;
        Point * sol_point_2;
        int iVs[FEPIC_MAX_ICELLS];
        int *iVs_end, *iVs_it;
        Vector aux(dim);

        iVs_end = mesh->connectedVtcs(point, iVs);

        // se esse nó está na linha, então existe um vértice vizinho que está no sólido
        for (iVs_it = iVs; iVs_it != iVs_end ; ++iVs_it)
        {
          sol_point_2 = mesh->getNodePtr(*iVs_it);
          if ( is_in(sol_point_2->getTag(), solid_tags) )
            if( !sol_point || (sol_point_2->getTag() > sol_point->getTag()) )
              sol_point = sol_point_2;
        }
        if (!sol_point)
        {
          //FEP_PRAGMA_OMP(critical)
          {
            printf("ERRO: ponto na linha tríplice não tem um vértice vizinho no sólido");
            throw;
          }
        }
        point->getCoord(Xqp.data(),dim);

        normal = solid_normal(Xqp, current_time, tag);

        sol_point->getCoord(aux.data(),dim);


        // choose a line_normal candidate
        line_normal(0) = -normal(1);
        line_normal(1) =  normal(0);

        // check orientation
        if (line_normal.dot(Xqp-aux) < 0)
          line_normal *= -1;
      }

      
      for (int qp = 0; qp < n_qpts_corner; ++qp)
      {
        if (dim==3)
        {
          F_r_mid   = x_coefs_r_mid_trans * dLqsi_r[qp];
          J_mid = F_r_mid.norm();
          weight  = quadr_corner->weight(qp);
          JxW_mid = J_mid*weight;
        }
        else
        {
          J_mid = 1;
          weight  = 1;
          JxW_mid = 1;
        }
        //invF_r_mid = F_r_mid.transpose()/(J_mid*J_mid);


        Xqp = x_coefs_r_mid_trans * qsi_r[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
        //dxphi_r = dLphi_r[qp] * invF_r_mid;
        //dxU_r   = u_coefs_r_mid_trans * dxphi_r; // n+utheta
        Uqp  = u_coefs_r_mid_trans * phi_r[qp];

        gama_mid = gama(Xqp, current_time+dt*utheta, tag);

        if (dim==3)
        {
          normal  = solid_normal(Xqp, current_time, tag);
          line_normal(0)= F_r_mid(0,0);
          line_normal(1)= F_r_mid(1,0);
          line_normal(2)= F_r_mid(2,0);
          line_normal *= line_normal_sign;
          line_normal = cross(line_normal, normal);
          line_normal.normalize();
        }

        double const Uqp_norm = abs(line_normal.dot(Uqp));
        for (int i = 0; i < n_dofs_u_per_corner/dim; ++i)
        {
          for (int c = 0; c < dim; ++c)
          {
            FUloc(i*dim + c) += JxW_mid*(-gama_mid*cos_theta0() + zeta(Uqp_norm,0)*line_normal.dot(Uqp))*line_normal(c)*phi_r[qp][i];
            
            //FUloc(i*dim + c) += JxW_mid*(-gama_mid*cos_theta0() )*line_normal(c)*phi_r[qp][i];
            //FUloc(i*dim + c) += JxW_mid*zeta(Uqp_norm,0)*Uqp(c)*phi_r[qp][i];

            for (int j = 0; j < n_dofs_u_per_corner/dim; ++j)
              for (int d = 0; d < dim; ++d)
                Aloc_r(i*dim + c, j*dim + d) += JxW_mid* zeta(Uqp_norm,0)*utheta*line_normal(d) *phi_r[qp][j]*line_normal(c)*phi_r[qp][i];
                //if (c==d)
                //  Aloc_r(i*dim + c, j*dim + d) += JxW_mid* zeta(Uqp_norm,0) *utheta*phi_r[qp][j]*phi_r[qp][i];

          }

          //FUloc.segment(i*dim, dim) += JxW_mid* phi_r[qp][i] * (-gama_mid*cos_theta0() + zeta(0,0)*line_normal.dot(Uqp))*line_normal;
        }


      } // end quadratura




      // Projection - to force non-penetrarion bc
      mesh->getCornerNodesId(&*corner, corner_nodes.data());
      getProjectorMatrix(Prj, nodes_per_corner, corner_nodes.data(), Vec_x_1, current_time+dt, *this);

      FUloc = Prj*FUloc;
      Aloc_r = Prj*Aloc_r*Prj;

      VecSetValues(Vec_fun, mapU_r.size(), mapU_r.data(), FUloc.data(), ADD_VALUES);
      MatSetValues(*JJ, mapU_r.size(), mapU_r.data(), mapU_r.size(), mapU_r.data(), Aloc_r.data(),  ADD_VALUES);
      //cout << FUloc.transpose() << endl;

    }



  }  //end LINHA DE CONTATO

  // boundary conditions on global Jacobian
    // solid & triple tags .. force normal
  if (force_dirichlet)
  {
    int      nodeid;
    int      u_dofs[dim];
    Vector   normal(dim);
    Tensor   A(dim,dim);
    Tensor   I(Tensor::Identity(dim,dim));
    int      tag;

    point_iterator point = mesh->pointBegin();
    point_iterator point_end = mesh->pointEnd();
    for ( ; point != point_end; ++point)
    {
      tag = point->getTag();
      if (!(is_in(tag,feature_tags)   ||
            is_in(tag,solid_tags)     ||
            is_in(tag,interface_tags) ||
            is_in(tag,triple_tags)    ||
            is_in(tag,dirichlet_tags) ||
            is_in(tag,neumann_tags)   ||
            is_in(tag,periodic_tags)     )  )
        continue;
      //dof_handler[DH_UNKS].getVariable(VAR_U).getVertexAssociatedDofs(u_dofs, &*point);
      getNodeDofs(&*point, DH_UNKS, VAR_U, u_dofs);

      nodeid = mesh->getPointId(&*point);
      getProjectorMatrix(A, 1, &nodeid, Vec_x_1, current_time+dt, *this);
      A = I - A;
      MatSetValues(*JJ, dim, u_dofs, dim, u_dofs, A.data(), ADD_VALUES);
    }
  }  //end force_dirichlet

  if (force_pressure)
  {
    double const p =1.0;
    MatSetValues(*JJ, 1, &null_space_press_dof, 1, &null_space_press_dof, &p, ADD_VALUES);
  }

  if(print_to_matlab)
  {
    static bool ja_foi=false;
    if (!ja_foi)
    {
      View(Vec_fun, "rhs.m","res");
      View(*JJ,"jacob.m","Jac");
    }
    ja_foi = true;

  }

  Assembly(Vec_fun);
  Assembly(*JJ);

  PetscFunctionReturn(0);

} // END formFunction


PetscErrorCode AppCtx::formJacobian(SNES snes,Vec Vec_up_k, Mat* /*Mat_Jac*/, Mat* /*prejac*/, MatStructure * /*flag*/)
{
  PetscBool          found = PETSC_FALSE;
  char               snes_type[PETSC_MAX_PATH_LEN];

  PetscOptionsGetString(PETSC_NULL,"-snes_type",snes_type,PETSC_MAX_PATH_LEN-1,&found);

  if (found)
    if (string(snes_type) == string("test"))
    {
      cout << "WARNING: TESTING JACOBIAN !!!!! \n";
      this->formFunction(snes, Vec_up_k, Vec_res);
    }

  PetscFunctionReturn(0);
}

// ====================================================================================================

// to apply boundary conditions on linear elasticity problem.
template <typename Derived>
void getProjectorBC(MatrixBase<Derived> & P, int n_nodes, int const* nodes, Vec const& Vec_x_, double t, AppCtx const& app)
{
  int const dim = app.dim;
  Mesh const* mesh = &*app.mesh;
  //DofHandler const* dof_handler = &*app.dof_handler;
  std::vector<int> const& dirichlet_tags  = app.dirichlet_tags;
  std::vector<int> const& neumann_tags    = app.neumann_tags  ;
  std::vector<int> const& interface_tags  = app.interface_tags;
  std::vector<int> const& solid_tags      = app.solid_tags    ;
  std::vector<int> const& triple_tags     = app.triple_tags   ;
  std::vector<int> const& periodic_tags   = app.periodic_tags ;
  std::vector<int> const& feature_tags    = app.feature_tags  ;
  Vec const& Vec_normal = app.Vec_normal;

  P.setIdentity();

  Tensor I(dim,dim);
  Tensor Z(dim,dim);
  Vector X(dim);
  Vector normal(dim);
  int    dofs[dim];
  int    tag;
  Point const* point;

  I.setIdentity();
  Z.setZero();

  bool boundary_smoothing = app.boundary_smoothing;
  //bool boundary_smoothing = false;

  // NODES
  for (int i = 0; i < n_nodes; ++i)
  {
    point = mesh->getNodePtr(nodes[i]);
    tag = point->getTag();
    //m = point->getPosition() - mesh->numVerticesPerCell();
    //cell = mesh->getCellPtr(point->getIncidCell());

    if (is_in(tag,feature_tags))
    {
      if (boundary_smoothing)
      {      
        app.getNodeDofs(&*point, DH_MESH, VAR_M, dofs);
        VecGetValues(Vec_x_, dim, dofs, X.data());
        P.block(i*dim,i*dim,dim,dim)  = feature_proj(X,t,tag);
      }
      else
        P.block(i*dim,i*dim,dim,dim) = Z;
    }
    else
    if (is_in(tag,solid_tags) )
    {
      if (boundary_smoothing)
      {
        app.getNodeDofs(&*point, DH_MESH, VAR_M, dofs);
        VecGetValues(Vec_x_, dim, dofs, X.data());
        normal = -solid_normal(X,t,tag);
        P.block(i*dim,i*dim,dim,dim)  = I - normal*normal.transpose();
      }
      else
      {
        P.block(i*dim,i*dim,dim,dim) = Z;
      }
    }
    else
    if (is_in(tag,interface_tags))
    {
      if (false && boundary_smoothing)
      {
        app.getNodeDofs(&*point, DH_MESH, VAR_M, dofs);
        VecGetValues(Vec_normal, dim, dofs, X.data());
        P.block(i*dim,i*dim,dim,dim) = I - X*X.transpose();
      }
      else
      {
        P.block(i*dim,i*dim,dim,dim) = Z;
      }
    }
    else
    //if (is_in(tag,triple_tags) || is_in(tag,dirichlet_tags) || is_in(tag,neumann_tags) || is_in(tag,periodic_tags) || is_in(tag,feature_tags))
    if (is_in(tag,triple_tags) || is_in(tag,dirichlet_tags) || is_in(tag,neumann_tags) || is_in(tag,periodic_tags))
    {
      P.block(i*dim,i*dim,dim,dim) = Z;
    }


  } // end nodes
}


// function to compute mesh velocity
PetscErrorCode AppCtx::formFunction_mesh(SNES /*snes_m*/, Vec Vec_v, Vec Vec_fun)
{
  Mat *JJ = &Mat_Jac_m;

  double utheta = AppCtx::utheta;
  
  if (is_bdf2)
  {
    if (time_step == 0)
      if (!is_bdf_euler_start)
        utheta = 0.5;
  }
  else
  if (is_bdf3)
  {
    if (time_step <= 1)
      utheta = 0.5;
  }
    
  //SNESGetJacobian(snes_m, JJ, NULL, NULL, NULL);

  // NOTE: solve elasticity problem in the mesh at time step n
  // NOTE: The mesh used is the Vec_x_0
  // WARNING: this function assumes that the boundary conditions was already applied
//PetscReal *vall;  Vec normv;  VecCopy(Vec_x_0,normv);  VecAXPY(normv,-1,Vec_x_1);  VecNorm(normv,NORM_2,vall);  cout << "\n " << vall << endl;

  VecZeroEntries(Vec_fun);
  MatZeroEntries(*JJ);

// LOOP NAS CÉLULAS Parallel (uncomment it)
//#ifdef FEP_HAS_OPENMP
//  FEP_PRAGMA_OMP(parallel default(none) shared(Vec_v, Vec_fun, cout, JJ, utheta))
//#endif
  {
    bool const non_linear = nonlinear_elasticity;

    Tensor            dxV(dim,dim);   // grad u
    Tensor            F_c(dim,dim);
    Tensor            invF_c(dim,dim);
    Tensor            invFT_c(dim,dim);
    Vector            Vqp(dim);
    MatrixXd          v_coefs_c_trans(dim, nodes_per_cell);      // mesh velocity;
    MatrixXd          v_coefs_c(nodes_per_cell, dim);
    MatrixXd          x_coefs_c_trans(dim, nodes_per_cell);
    MatrixXd          x_coefs_c(nodes_per_cell, dim);
    MatrixXd          x_coefs_c_new_trans(dim, nodes_per_cell);
    MatrixXd          x_coefs_c_new(nodes_per_cell, dim);
    MatrixXd          dxqsi_c(nodes_per_cell, dim);
    double            J, weight, JxW;

    VectorXd          Floc(n_dofs_v_per_cell);
    MatrixXd          Aloc(n_dofs_v_per_cell, n_dofs_v_per_cell);

    VectorXi          mapV_c(n_dofs_u_per_cell);  //TODO: i think is n_dofs_v_per_cell

    MatrixXd          Prj(n_dofs_v_per_cell, n_dofs_v_per_cell);
    VectorXi          cell_nodes(nodes_per_cell);

    double            sigma_ck;
    double            dsigma_ckjd;  // dphi, d_compd

    const int tid = omp_get_thread_num();
    const int nthreads = omp_get_num_threads();

    cell_iterator cell = mesh->cellBegin(tid,nthreads);
    cell_iterator cell_end = mesh->cellEnd(tid,nthreads);

    //cell_iterator cell = mesh->cellBegin();
    //cell_iterator cell_end = mesh->cellEnd();
    for (; cell != cell_end; ++cell)
    {

      // mapeamento do local para o global: (ID local para ID global)
      // mapV_c saves global IDs for cell's nodes unknowns enumeration
      dof_handler[DH_MESH].getVariable(VAR_M).getCellDofs(mapV_c.data(), &*cell);  //cout << mapV_c << endl;

      /*  Pega os valores das variáveis nos graus de liberdade */
      VecGetValues(Vec_v ,  mapV_c.size(), mapV_c.data(), v_coefs_c.data());  //cout << v_coefs_c << endl;//VecView(Vec_v,PETSC_VIEWER_STDOUT_WORLD);
      VecGetValues(Vec_x_0, mapV_c.size(), mapV_c.data(), x_coefs_c.data());  //cout << x_coefs_c << endl;
      VecGetValues(Vec_x_1, mapV_c.size(), mapV_c.data(), x_coefs_c_new.data());  //cout << x_coefs_c_new << endl;

      //if (current_time < .5*dt)
      if ((is_bdf2 && time_step > 0) || (is_bdf3 && time_step > 1) )
        x_coefs_c = x_coefs_c_new;
      else
      {
        x_coefs_c = (1.-utheta)*x_coefs_c + utheta*x_coefs_c_new;
        //x_coefs_c += x_coefs_c_new;
        //x_coefs_c /= 2.;
      }

      v_coefs_c_trans = v_coefs_c.transpose();
      x_coefs_c_trans = x_coefs_c.transpose();

      Floc.setZero();
      Aloc.setZero();

      // Quadrature
      for (int qp = 0; qp < n_qpts_cell; ++qp)
      {
        F_c = x_coefs_c_trans * dLqsi_c[qp];  //cout << dLqsi_c[qp] << endl;
        inverseAndDet(F_c,dim,invF_c,J);
        invFT_c= invF_c.transpose();  //usado?

        dxqsi_c = dLqsi_c[qp] * invF_c;

        dxV  = v_coefs_c_trans * dxqsi_c;       // n+utheta
        Vqp  = v_coefs_c_trans * qsi_c[qp];
        //Xqp      = x_coefs_c_trans * qsi_c[qp]; // coordenada espacial (x,y,z) do ponto de quadratura

        weight = quadr_cell->weight(qp);
        JxW = J*weight;  //parece que no es necesario, ver 2141 (JxW/JxW)


        for (int i = 0; i < n_dofs_v_per_cell/dim; ++i)  //sobre cantidad de funciones de forma
        {
          for (int c = 0; c < dim; ++c)  //sobre dimension
          {
            for (int k = 0; k < dim; ++k)  //sobre dimension
            {
              sigma_ck = dxV(c,k) + dxV(k,c);

              if (non_linear)  //is right?
              {
                for (int l = 0; l < dim; ++l)
                {
                  sigma_ck += dxV(l,c)*dxV(l,k);  // i think is dxV(c,l)*dxV(k,l);

                  if (c==k)
                  {
                    sigma_ck -= dxV(l,l);
                    for (int m = 0; m < dim; ++m)
                      sigma_ck -=  dxV(l,m)*dxV(l,m);
                  }

                }
              }  //end non_linear

              Floc(i*dim + c) += sigma_ck*dxqsi_c(i,k) * (JxW/JxW); // (JxW/JxW) is to compiler not complain about unused variables

              for (int j = 0; j < n_dofs_v_per_cell/dim; ++j)
              {
                for (int d = 0; d < dim; ++d)
                {
                  dsigma_ckjd = 0;

                  if (c==d)
                    dsigma_ckjd = dxqsi_c(j,k);

                  if (k==d)
                    dsigma_ckjd += dxqsi_c(j,c);

                  if (non_linear)  //is right?
                  {
                    for (int l = 0; l < dim; ++l)
                    {
                      if (l==d)
                        dsigma_ckjd += dxqsi_c(j,c)*dxV(l,k) + dxV(l,c)*dxqsi_c(j,k);  //is ok?

                      if (c==k)
                      {
                        if (l==d)
                        {
                          dsigma_ckjd -= dxqsi_c(j,l);
                          for (int m = 0; m < dim; ++m)
                            dsigma_ckjd -= 2.*dxqsi_c(j,m)*dxV(l,m);
                        }
                      }
                    }
                  }  //end non_linear

                  Aloc(i*dim + c, j*dim + d) += dsigma_ckjd*dxqsi_c(i,k);

                } // end d

              } // end j

            } // end k

          }// end c
        } // end i


      } // fim quadratura


      // Projection - to force non-penetrarion bc
      mesh->getCellNodesId(&*cell, cell_nodes.data());
      getProjectorBC(Prj, nodes_per_cell, cell_nodes.data(), Vec_x_0, current_time, *this /*AppCtx*/);

      Floc = Prj*Floc;
      Aloc = Prj*Aloc*Prj;

#ifdef FEP_HAS_OPENMP
      FEP_PRAGMA_OMP(critical)
#endif
      {
        VecSetValues(Vec_fun, mapV_c.size(), mapV_c.data(), Floc.data(), ADD_VALUES);
        MatSetValues(*JJ, mapV_c.size(), mapV_c.data(), mapV_c.size(), mapV_c.data(), Aloc.data(),  ADD_VALUES);
      }
    } // end cell loop


  } // end parallel


  // boundary conditions on global Jacobian
    // solid & triple tags .. force normal
  if (force_dirichlet)  //identify the contribution of points in *_tags
  {
    int      nodeid;
    int      v_dofs[dim];
    Vector   normal(dim);
    Tensor   A(dim,dim);
    Tensor   I(Tensor::Identity(dim,dim));
    int      tag;

    point_iterator point = mesh->pointBegin();
    point_iterator point_end = mesh->pointEnd();
    for ( ; point != point_end; ++point)
    {
      tag = point->getTag();
      if (!(is_in(tag,feature_tags)   ||
            is_in(tag,solid_tags)     ||
            is_in(tag,interface_tags) ||
            is_in(tag,triple_tags)    ||
            is_in(tag,dirichlet_tags) ||
            is_in(tag,neumann_tags)   ||
            is_in(tag,periodic_tags)     )  )
        continue;
      //dof_handler[DH_UNKS].getVariable(VAR_U).getVertexAssociatedDofs(v_dofs, &*point);
      getNodeDofs(&*point, DH_MESH, VAR_M, v_dofs);

      nodeid = mesh->getPointId(&*point);
      getProjectorBC(A, 1, &nodeid, Vec_x_0, current_time, *this);
      A = I - A;
      MatSetValues(*JJ, dim, v_dofs, dim, v_dofs, A.data(), ADD_VALUES);
    }
  }

  Assembly(*JJ); //MatView(*JJ,PETSC_VIEWER_STDOUT_WORLD);
  Assembly(Vec_fun);

  //View(*JJ, "ElastOp", "JJ");

  PetscFunctionReturn(0);
}


PetscErrorCode AppCtx::formJacobian_mesh(SNES /*snes*/,Vec /*Vec_up_k*/,Mat* /**Mat_Jac*/, Mat* /*prejac*/, MatStructure * /*flag*/)
{
  // jacobian matrix is done in the formFunction_mesh
  PetscFunctionReturn(0);
}



// ******************************************************************************
//                            FORM FUNCTION_FS
// ******************************************************************************
PetscErrorCode AppCtx::formFunction_fs(SNES /*snes*/, Vec Vec_uzp_k, Vec Vec_fun_fs)
{

  double utheta = AppCtx::utheta;

  if (is_bdf2)
  {
    if (time_step == 0)
      if(!is_bdf_euler_start)
        utheta = 0.5;
  }
  else if (is_bdf3)
  {
    if (time_step <= 1)
      utheta = 0.5;
  }

  bool const compact_bubble = false; // eliminate bubble from convective term


  //PetscErrorCode      ierr;

  int null_space_press_dof=-1;

  int iter;

  SNESGetIterationNumber(snes_fs,&iter);  //cout << iter <<endl;

  if (!iter)
  {
    converged_times=0;
  }

  if (force_pressure && (iter<2))
  {
    Vector X(dim);
    Vector X_new(dim);
    if (behaviors & BH_Press_grad_elim)
    {
      cell_iterator cell = mesh->cellBegin();
      dof_handler[DH_UNKS].getVariable(VAR_P).getCellDofs(&null_space_press_dof, &*cell);
      // fix the initial guess
      VecSetValue(Vec_uzp_k, null_space_press_dof, 0.0, INSERT_VALUES);
    }
    else
    {
      point_iterator point = mesh->pointBegin();
      while (!mesh->isVertex(&*point))
        ++point;
      int x_dofs[3];
      dof_handler[DH_MESH].getVariable(VAR_M).getVertexDofs(x_dofs, &*point);
      VecGetValues(Vec_x_1, dim, x_dofs, X_new.data());
      VecGetValues(Vec_x_0, dim, x_dofs, X.data());
      X = .5*(X+X_new);
      dof_handler[DH_UNKS].getVariable(VAR_P).getVertexDofs(&null_space_press_dof, &*point);
      // fix the initial guess
      VecSetValue(Vec_uzp_k, null_space_press_dof, p_exact(X,current_time+.5*dt,point->getTag()), INSERT_VALUES);
    }

    Assembly(Vec_uzp_k);

  }

  // checking:
  if (null_space_press_dof < 0 && force_pressure==1 && (iter<2))
  {
    cout << "force_pressure: something is wrong ..." << endl;
    throw;
  }


  Mat *JJ = &Mat_Jac_fs;



  //PetscErrorCode      ierr;

  VecZeroEntries(Vec_fun_fs);
  MatZeroEntries(*JJ);

  // LOOP NAS CÉLULAS Parallel (uncomment it)
//#ifdef FEP_HAS_OPENMP
//  FEP_PRAGMA_OMP(parallel default(none) shared(Vec_uzp_k,Vec_fun_fs,cout,null_space_press_dof,JJ,utheta,iter))
//#endif
  {
    VectorXd          FUloc(n_dofs_u_per_cell);  // U subvector part of F
    VectorXd          FPloc(n_dofs_p_per_cell);
    VectorXd          FZloc(nodes_per_cell*3);

    /* local data */
    int                 tag, nod_id;
    MatrixXd            u_coefs_c_mid_trans(dim, n_dofs_u_per_cell/dim);  // n+utheta  // trans = transpost
    MatrixXd            u_coefs_c_old(n_dofs_u_per_cell/dim, dim);        // n
    MatrixXd            u_coefs_c_old_trans(dim,n_dofs_u_per_cell/dim);   // n
    MatrixXd            u_coefs_c_new(n_dofs_u_per_cell/dim, dim);        // n+1
    MatrixXd            u_coefs_c_new_trans(dim,n_dofs_u_per_cell/dim);   // n+1

    MatrixXd            du_coefs_c_old(n_dofs_u_per_cell/dim, dim);        // n
    MatrixXd            du_coefs_c_old_trans(dim,n_dofs_u_per_cell/dim);   // n

    MatrixXd            du_coefs_c_vold(n_dofs_u_per_cell/dim, dim);        // n-1
    MatrixXd            du_coefs_c_vold_trans(dim,n_dofs_u_per_cell/dim);   // n-1

    MatrixXd            v_coefs_c_mid(nodes_per_cell, dim);        // mesh velocity; n
    MatrixXd            v_coefs_c_mid_trans(dim,nodes_per_cell);   // mesh velocity; n

    VectorXd            p_coefs_c_new(n_dofs_p_per_cell);  // n+1
    VectorXd            p_coefs_c_old(n_dofs_p_per_cell);  // n
    VectorXd            p_coefs_c_mid(n_dofs_p_per_cell);  // n

    MatrixXd            z_coefs_c_mid_trans(3, n_dofs_z_per_cell/3);  // n+utheta  // trans = transpost
    MatrixXd            z_coefs_c_old(n_dofs_z_per_cell/3, 3);        // n
    MatrixXd            z_coefs_c_old_trans(3,n_dofs_z_per_cell/3);   // n
    MatrixXd            z_coefs_c_new(n_dofs_z_per_cell/3, 3);        // n+1
    MatrixXd            z_coefs_c_new_trans(3,n_dofs_z_per_cell/3);   // n+1
    MatrixXd            z_coefs_c_auo(dim, nodes_per_cell), z_coefs_c_aun(dim, nodes_per_cell);

    //MatrixXd            x_coefs_c_mid(nodes_per_cell, dim);       // n+utheta
    MatrixXd            x_coefs_c_mid_trans(dim, nodes_per_cell); // n+utheta
    MatrixXd            x_coefs_c_new(nodes_per_cell, dim);       // n+1
    MatrixXd            x_coefs_c_new_trans(dim, nodes_per_cell); // n+1
    MatrixXd            x_coefs_c_old(nodes_per_cell, dim);       // n
    MatrixXd            x_coefs_c_old_trans(dim, nodes_per_cell); // n

    Tensor              F_c_mid(dim,dim);       // n+utheta
    Tensor              invF_c_mid(dim,dim);    // n+utheta
    Tensor              invFT_c_mid(dim,dim);   // n+utheta

    Tensor              F_c_old(dim,dim);       // n
    Tensor              invF_c_old(dim,dim);    // n
    Tensor              invFT_c_old(dim,dim);   // n

    Tensor              F_c_new(dim,dim);       // n+1
    Tensor              invF_c_new(dim,dim);    // n+1
    Tensor              invFT_c_new(dim,dim);   // n+1


    //Tensor              F_c_new(dim,dim);       // n+1
    //Tensor              invF_c_new(dim,dim);    // n+1
    //Tensor              invFT_c_new(dim,dim);   // n+1
    //Tensor              F_c_old(dim,dim);       // n
    //Tensor              invF_c_old(dim,dim);    // n
    //Tensor              invFT_c_old(dim,dim);   // n

    /* All variables are in (n+utheta) by default */

    MatrixXd            dxphi_c(n_dofs_u_per_cell/dim, dim);
    MatrixXd            dxphi_c_new(dxphi_c);
    MatrixXd            dxpsi_c(n_dofs_p_per_cell, dim);
    MatrixXd            dxpsi_c_new(dxpsi_c);
    MatrixXd            dxqsi_c(nodes_per_cell, dim);
    Vector              dxbble(dim);
    Vector              dxbble_new(dim);
    Tensor              dxU(dim,dim);   // grad u
    Tensor              dxU_old(dim,dim);   // grad u
    Tensor              dxU_new(dim,dim);   // grad u
    Tensor              dxUb(dim,dim);  // grad u bble
    Vector              dxP_new(dim);   // grad p
    Vector              Xqp(dim);
    Vector              Xqp_old(dim);
    Vector              Xc(dim);  // cell center; to compute CR element
    Vector              Uqp(dim);
    Vector              Ubqp(dim); // bble
    Vector              Uqp_old(dim), Zqp_old(dim); // n
    Vector              Uqp_new(dim), Zqp_new(dim); // n+1
    Vector              dUqp_old(dim);  // n
    Vector              dUqp_vold(dim);  // n
    Vector              Vqp(dim);
    Vector              Uconv_qp(dim);
    Vector              dUdt(dim);
    double              Pqp_new;
    double              Pqp;
    double              bble_integ=0;
    //VectorXd          FUloc(n_dofs_u_per_cell); // subvetor da função f (parte de U)
    //VectorXd          FPloc(n_dofs_p_per_cell);     // subvetor da função f (parte de P)
    VectorXi            cell_nodes(nodes_per_cell);
    double              J_mid;
    double              J_new, J_old;
    double              JxW_mid;
    double              JxW_new, JxW_old;
    double              weight;
    double              visc=-1; // viscosity
    double              cell_volume;
    double              hk2;
    double              tauk=0;
    double              delk=0;
    double              delta_cd;
    double              rho;
    double ddt_factor;
    if (is_bdf2 && time_step > 0)
      ddt_factor = 1.5;
    else
    if (is_bdf3 && time_step > 1)
      ddt_factor = 11./6.;
    else
      ddt_factor = 1.;


    MatrixXd            Aloc(n_dofs_u_per_cell, n_dofs_u_per_cell);
    MatrixXd            Gloc(n_dofs_u_per_cell, n_dofs_p_per_cell);
    MatrixXd            Dloc(n_dofs_p_per_cell, n_dofs_u_per_cell);
    MatrixXd            Eloc(n_dofs_p_per_cell, n_dofs_p_per_cell);   // GSL, BC
    MatrixXd            Cloc(n_dofs_u_per_cell, n_dofs_p_per_cell);   // GSL
    Tensor              iBbb(dim, dim);                               // BC, i : inverse ..it is not the inverse to CR element
    MatrixXd            Bbn(dim, n_dofs_u_per_cell);                  // BC
    MatrixXd            Bnb(n_dofs_u_per_cell, dim);                  // BC
    MatrixXd            Dpb(n_dofs_p_per_cell, dim);                  // BC
    MatrixXd            Gbp(dim, n_dofs_p_per_cell);                  // BC
    MatrixXd            Gnx(n_dofs_u_per_cell, dim);                  // CR ;; suffix x means p gradient
    Vector              FUb(dim);                                     // BC
    Vector              FPx(dim); // pressure gradient

    MatrixXd            Z1loc(n_dofs_u_per_cell,nodes_per_cell*3);
    MatrixXd            Z2loc(nodes_per_cell*3,n_dofs_u_per_cell);
    MatrixXd            Z3loc(nodes_per_cell*3,nodes_per_cell*3);
    MatrixXd            Z4loc(nodes_per_cell*3,n_dofs_p_per_cell);
    MatrixXd            Z5loc(n_dofs_p_per_cell,nodes_per_cell*3);

    Vector              force_at_mid(dim);
    Vector              Res(dim);                                     // residue
    Tensor              dResdu(dim,dim);                              // residue derivative
    Tensor const        I(Tensor::Identity(dim,dim));
    Vector              vec(dim);     // temp
    Tensor              Ten(dim,dim); // temp

    VectorXi            mapU_c(n_dofs_u_per_cell);
    VectorXi            mapU_r(n_dofs_u_per_corner);
    VectorXi            mapP_c(n_dofs_p_per_cell);
    VectorXi            mapP_r(n_dofs_p_per_corner);
    VectorXi            mapZ_c(nodes_per_cell*3);
    VectorXi            mapZ_f(nodes_per_facet*3);
    // mesh velocity
    VectorXi            mapM_c(dim*nodes_per_cell);
    VectorXi            mapM_f(dim*nodes_per_facet);
    VectorXi            mapM_r(dim*nodes_per_corner);

    MatrixXd            Prj(n_dofs_u_per_cell,n_dofs_u_per_cell); // projector matrix
    //VectorXi            cell_nodes(nodes_per_cell);

    //Q Dofs re-organization
    VectorXi p_id_cor(3); p_id_cor << 1,1,1;
    int n_solid_nodes = dof_handler[DH_UNKM].getVariable(VAR_Z).numPositiveDofs()/3;
    p_id_cor = 3*(n_solid_nodes-N_Solids)*p_id_cor;

    std::vector<bool>   SV(N_Solids,false);  //solid visited history
    std::vector<int>    SV_c(nodes_per_cell,0);   //maximum nodes in solid visited saving the solid tag

    Vector Rotf(dim), Rotv(dim);
    Vector Zw(3); Zw << 0,0,1;
    Vector Xg(dim);
    Vector2d XIp;
    Vector auxRotf(dim), auxRotv(dim);
    Tensor auxTenf(dim,dim), auxTenv(dim,dim);

    const int tid = omp_get_thread_num();
    const int nthreads = omp_get_num_threads();

    cell_iterator cell = mesh->cellBegin(tid,nthreads);
    cell_iterator cell_end = mesh->cellEnd(tid,nthreads);

    //cell_iterator cell = mesh->cellBegin();
    //cell_iterator cell_end = mesh->cellEnd();
    for (; cell != cell_end; ++cell)
    {

      tag = cell->getTag();
      //cout << p_coefs_c_new << endl << endl; cout << z_coefs_c_new << endl << endl;
      // mapeamento do local para o global:
      //
      dof_handler[DH_MESH].getVariable(VAR_M).getCellDofs(mapM_c.data(), &*cell);  //cout << mapM_c << endl;  //unk. global ID's
      dof_handler[DH_UNKM].getVariable(VAR_U).getCellDofs(mapU_c.data(), &*cell);  //cout << mapU_c << endl;
      dof_handler[DH_UNKM].getVariable(VAR_Q).getCellDofs(mapP_c.data(), &*cell);  //cout << mapP_c << endl;
      mapP_c = mapP_c - p_id_cor;

      dof_handler[DH_UNKM].getVariable(VAR_Z).getCellDofs(mapZ_c.data(), &*cell);  // Z global ids for the current cell
      //cout << mapZ_c << endl;
      for (int j = 0; j < cell->numNodes(); ++j){
        tag = mesh->getNodePtr(cell->getNodeId(j))->getTag();
        nod_id = is_in_id(tag,flusoli_tags);
        if (nod_id){
          //if(!SV[nod_id-1]){
            for (int l = 0; l < 3; l++){  // the 3 here is for Z quantity of Dofs for 2D case
              mapZ_c(j*cell->numNodes() + l) = dof_handler[DH_UNKM].getVariable(VAR_U).numPositiveDofs() - 1
          			                         + 3*nod_id - 2 + l;
            }
            //SV[nod_id-1] = true;  //cout << "Solid " << nod_id << " visited." << endl;
          //}
        }
        SV_c[j]=nod_id;
      }
      //for (int j = 0; j < nodes_per_cell; j++) {cout << SV_c[j] << " ";} cout << endl;
      //cout << mapZ_c.transpose() << endl; //VecSetOption(Vec_uzp_0, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
      u_coefs_c_old = MatrixXd::Zero(n_dofs_u_per_cell/dim, dim);
      u_coefs_c_new = MatrixXd::Zero(n_dofs_u_per_cell/dim, dim);
      z_coefs_c_old = MatrixXd::Zero(n_dofs_z_per_cell/3, 3);
      z_coefs_c_new = MatrixXd::Zero(n_dofs_z_per_cell/3, 3);
      z_coefs_c_auo = MatrixXd::Zero(dim,nodes_per_cell);
      z_coefs_c_aun = MatrixXd::Zero(dim,nodes_per_cell);
      VecSetOption(Vec_uzp_0, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
      VecSetOption(Vec_uzp_k, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);

      if ((is_bdf2 && time_step > 0) || (is_bdf3 && time_step > 1))
        VecGetValues(Vec_v_1, mapM_c.size(), mapM_c.data(), v_coefs_c_mid.data());
      else
        VecGetValues(Vec_v_mid, mapM_c.size(), mapM_c.data(), v_coefs_c_mid.data());  //cout << v_coefs_c_mid << endl << endl;//size of vector mapM_c
      VecGetValues(Vec_x_0,     mapM_c.size(), mapM_c.data(), x_coefs_c_old.data());  //cout << x_coefs_c_old << endl << endl;
      VecGetValues(Vec_x_1,     mapM_c.size(), mapM_c.data(), x_coefs_c_new.data());  //cout << x_coefs_c_new << endl << endl;
      VecGetValues(Vec_uzp_0,    mapU_c.size(), mapU_c.data(), u_coefs_c_old.data()); //cout << u_coefs_c_old << endl << endl;
      VecGetValues(Vec_uzp_k,    mapU_c.size(), mapU_c.data(), u_coefs_c_new.data()); //cout << u_coefs_c_new << endl << endl;
      VecGetValues(Vec_uzp_0,    mapP_c.size(), mapP_c.data(), p_coefs_c_old.data()); //cout << p_coefs_c_old << endl << endl;
      VecGetValues(Vec_uzp_k,    mapP_c.size(), mapP_c.data(), p_coefs_c_new.data()); //cout << p_coefs_c_new << endl << endl;
      VecGetValues(Vec_uzp_0,    mapZ_c.size(), mapZ_c.data(), z_coefs_c_old.data()); //cout << z_coefs_c_old << endl << endl;
      VecGetValues(Vec_uzp_k,    mapZ_c.size(), mapZ_c.data(), z_coefs_c_new.data()); //cout << z_coefs_c_new << endl << endl;

      VecGetValues(Vec_duzp,    mapU_c.size(), mapU_c.data(), du_coefs_c_old.data()); // bdf2,bdf3
      if (is_bdf3)
        VecGetValues(Vec_dup_0,  mapU_c.size(), mapU_c.data(), du_coefs_c_vold.data()); // bdf3

      // get nodal coordinates of the old and new cell
      mesh->getCellNodesId(&*cell, cell_nodes.data());
      //mesh->getNodesCoords(cell_nodes.begin(), cell_nodes.end(), x_coefs_c.data());
      //x_coefs_c_trans = x_coefs_c_mid_trans;

      v_coefs_c_mid_trans = v_coefs_c_mid.transpose();
      x_coefs_c_old_trans = x_coefs_c_old.transpose();
      x_coefs_c_new_trans = x_coefs_c_new.transpose();
      u_coefs_c_old_trans = u_coefs_c_old.transpose();
      u_coefs_c_new_trans = u_coefs_c_new.transpose();
      z_coefs_c_old_trans = z_coefs_c_old.transpose();
      z_coefs_c_new_trans = z_coefs_c_new.transpose();

      du_coefs_c_old_trans = du_coefs_c_old.transpose(); // bdf2
      if (is_bdf3)
        du_coefs_c_vold_trans = du_coefs_c_vold.transpose(); // bdf3

      u_coefs_c_mid_trans = utheta*u_coefs_c_new_trans + (1.-utheta)*u_coefs_c_old_trans;
      x_coefs_c_mid_trans = utheta*x_coefs_c_new_trans + (1.-utheta)*x_coefs_c_old_trans;
      p_coefs_c_mid       = utheta*p_coefs_c_new       + (1.-utheta)*p_coefs_c_old;
      z_coefs_c_mid_trans = utheta*z_coefs_c_new_trans + (1.-utheta)*z_coefs_c_old_trans;

      visc = muu(tag);
      rho  = pho(Xqp,tag);
      Aloc.setZero();
      Gloc.setZero();
      Dloc.setZero();
      FUloc.setZero();
      FZloc.setZero();
      FPloc.setZero();
      Eloc.setZero();


      if (behaviors & BH_bble_condens_PnPn) // reset matrices
      {
        iBbb.setZero();
        Bnb.setZero();
        Gbp.setZero();
        FUb.setZero();
        Bbn.setZero();
        Dpb.setZero();
      }

      if(behaviors & BH_GLS)
      {
        cell_volume = 0;
        for (int qp = 0; qp < n_qpts_cell; ++qp) {
          F_c_mid = x_coefs_c_mid_trans * dLqsi_c[qp];
          J_mid = determinant(F_c_mid,dim);
          cell_volume += J_mid * quadr_cell->weight(qp);
        }  //cout << J_mid << " " << cell_volume << endl;

        hk2 = cell_volume / pi; // element size
        double const uconv = (u_coefs_c_old - v_coefs_c_mid).lpNorm<Infinity>();

        tauk = 4.*visc/hk2 + 2.*rho*uconv/sqrt(hk2);
        tauk = 1./tauk;
        if (dim==3)
          tauk *= 0.1;

        delk = 4.*visc + 2.*rho*uconv*sqrt(hk2);
        //delk = 0;

        Eloc.setZero();
        Cloc.setZero();
      }
      if (behaviors & BH_bble_condens_CR)
      {
        bble_integ = 0;
        Gnx.setZero();
        iBbb.setZero();
        Bnb.setZero();
        FUb.setZero();
        FPx.setZero();
        Bbn.setZero();

        cell_volume = 0;
        Xc.setZero();
        for (int qp = 0; qp < n_qpts_cell; ++qp) {
          F_c_mid = x_coefs_c_mid_trans * dLqsi_c[qp];
          J_mid = determinant(F_c_mid,dim);
          Xqp  = x_coefs_c_mid_trans * qsi_c[qp];
          cell_volume += J_mid * quadr_cell->weight(qp);
          Xc += J_mid * quadr_cell->weight(qp) * Xqp;
        }
        Xc /= cell_volume;
      }


      // Quadrature
      for (int qp = 0; qp < n_qpts_cell; ++qp)
      {
        //F_c_new = x_coefs_c_new_trans * dLqsi_c[qp];
        //inverseAndDet(F_c_new,dim,invF_c_new,J_new);
        //invFT_c_new= invF_c_new.transpose();

        //F_c_old = x_coefs_c_old_trans * dLqsi_c[qp];
        //inverseAndDet(F_c_old,dim,invF_c_old,J_old);
        //invFT_c_old= invF_c_old.transpose();

        F_c_mid = x_coefs_c_mid_trans * dLqsi_c[qp];  // (dim x nodes_per_cell) (nodes_per_cell x dim)
        F_c_old = x_coefs_c_old_trans * dLqsi_c[qp];
        F_c_new = x_coefs_c_new_trans * dLqsi_c[qp];
        inverseAndDet(F_c_mid,dim,invF_c_mid,J_mid);
        inverseAndDet(F_c_old,dim,invF_c_old,J_old);
        inverseAndDet(F_c_new,dim,invF_c_new,J_new);
        invFT_c_mid= invF_c_mid.transpose();
        invFT_c_old= invF_c_old.transpose();
        invFT_c_new= invF_c_new.transpose();

        dxphi_c_new = dLphi_c[qp] * invF_c_new;
        dxphi_c     = dLphi_c[qp] * invF_c_mid;
        dxpsi_c_new = dLpsi_c[qp] * invF_c_new;
        dxpsi_c     = dLpsi_c[qp] * invF_c_mid;
        dxqsi_c     = dLqsi_c[qp] * invF_c_mid;

        dxP_new  = dxpsi_c.transpose() * p_coefs_c_new;
        dxU      = u_coefs_c_mid_trans * dLphi_c[qp] * invF_c_mid;       // n+utheta
        dxU_new  = u_coefs_c_new_trans * dLphi_c[qp] * invF_c_new;       // n+1
        dxU_old  = u_coefs_c_old_trans * dLphi_c[qp] * invF_c_old;       // n

        Xqp      = x_coefs_c_mid_trans * qsi_c[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
        Xqp_old  = x_coefs_c_old_trans * qsi_c[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
        Uqp      = u_coefs_c_mid_trans * phi_c[qp]; //n+utheta
        Uqp_new  = u_coefs_c_new_trans * phi_c[qp]; //n+1
        Uqp_old  = u_coefs_c_old_trans * phi_c[qp]; //n
        Pqp_new  = p_coefs_c_new.dot(psi_c[qp]);
        Pqp      = p_coefs_c_mid.dot(psi_c[qp]);
        Vqp      = v_coefs_c_mid_trans * qsi_c[qp];
        //Vqp = v_exact(Xqp_old, current_time, tag);
        //Vqp = v_exact(Xqp, current_time+dt/2., tag);
        Uconv_qp = Uqp - Vqp;
        //Uconv_qp = Uqp_old;
        for (int l = 0; l < nodes_per_cell; l++){
          if (SV_c[l]){
            z_coefs_c_auo.col(l) = SolidVel(Xqp, XG_0[SV_c[l]-1], z_coefs_c_old_trans.col(l));
            z_coefs_c_aun.col(l) = SolidVel(Xqp, XG[SV_c[l]-1], z_coefs_c_new_trans.col(l));
          }
        }
        Zqp_new  = z_coefs_c_aun * phi_c[qp];
        Zqp_old  = z_coefs_c_auo * phi_c[qp];

        dUdt     = (Uqp_new-Uqp_old)/dt + (Zqp_new-Zqp_old)/dt;

        if (is_bdf2 && time_step > 0)
        {
          dUqp_old  = du_coefs_c_old_trans * phi_c[qp]; //n+utheta
          dUdt = 1.5*dUdt - .5*dUqp_old;
        }
        else
        if (is_bdf3 && time_step > 1)
        {
          dUqp_old   = du_coefs_c_old_trans  * phi_c[qp];
          dUqp_vold  = du_coefs_c_vold_trans * phi_c[qp];
          dUdt = 11./6.*dUdt - 7./6.*dUqp_old + 1./3.*dUqp_vold;
        }


        force_at_mid = force(Xqp,current_time+utheta*dt,tag);

        weight = quadr_cell->weight(qp);
        JxW_mid = J_mid*weight;
        JxW_old = J_old*weight;
        JxW_new = J_new*weight;
//cout << "J_mid = " << J_mid << endl;
        //~ if (mesh->getCellId(&*cell) == 0)
        //~ {
          //~ printf("cHEcKKKKKKKKKKK!!\n");
          //~ cout << "x coefs mid:" << endl;
          //~ cout << x_coefs_c_mid_trans.transpose() << endl;
        //~ }
        if (J_mid < 1.e-14)
        {
          FEP_PRAGMA_OMP(critical)
          //if (tid==0)
          {
            printf("in formCellFunction:\n");
            std::cout << "erro: jacobiana da integral não invertível: ";
            std::cout << "J_mid = " << J_mid << endl;
            cout << "trans(f) matrix:\n" << F_c_mid << endl;
            cout << "x coefs mid:" << endl;
            cout << x_coefs_c_mid_trans.transpose() << endl;
            cout << "-----" << endl;
            cout << "cell id: " << mesh->getCellId(&*cell) << endl;
            cout << "cell Contig id: " << mesh->getCellContigId( mesh->getCellId(&*cell) ) << endl;
            cout << "cell nodes:\n" << cell_nodes.transpose() << endl;
            cout << "mapM :\n" << mapM_c.transpose() << endl;
            throw;
          }
        }

        for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
        {
          for (int c = 0; c < dim; ++c)
          {
            FUloc(i*dim + c) += JxW_mid*
                 ( rho*(unsteady*dUdt(c) + has_convec*Uconv_qp.dot(dxU.row(c)))*phi_c[qp][i] + //aceleração
                   visc*dxphi_c.row(i).dot(dxU.row(c) + dxU.col(c).transpose()) - //rigidez  //transpose() here is to add 2 row-vectors
                   force_at_mid(c)*phi_c[qp][i] - //força
                   Pqp_new*dxphi_c(i,c) );        //pressão

            for (int j = 0; j < n_dofs_u_per_cell/dim; ++j)
            {
              for (int d = 0; d < dim; ++d)
              {
                delta_cd = c==d;
                Aloc(i*dim + c, j*dim + d) += JxW_mid*
                    ( has_convec*phi_c[qp][i]*utheta *rho*( delta_cd*Uconv_qp.dot(dxphi_c.row(j)) + dxU(c,d)*phi_c[qp][j] ) + //advecção
                      ddt_factor* unsteady* delta_cd*rho*phi_c[qp][i]*phi_c[qp][j]/dt + //time derivative
                      utheta*visc*( delta_cd * dxphi_c.row(i).dot(dxphi_c.row(j)) + dxphi_c(i,d)*dxphi_c(j,c)) ); //rigidez

              }
            }
            for (int j = 0; j < n_dofs_p_per_cell; ++j)
            {
              Gloc(i*dim + c,j) -= JxW_mid * psi_c[qp][j]* dxphi_c(i,c);
              //Gloc(i*dim + c,j) -= utheta*JxW_mid * psi_c[qp][j]* dxphi_c(i,c);
              //Gloc(i*dim + c,j) -= JxW_new * psi_c[qp][j]* dxphi_c_new(i,c);
              //Dloc(j, i*dim + c) -= JxW_new * psi_c[qp][j]*  dxphi_c_new(i,c);
              Dloc(j, i*dim + c) -= utheta*JxW_mid * psi_c[qp][j]*  dxphi_c(i,c);
            }

          }
        }

        for (int i = 0; i < n_dofs_p_per_cell; ++i)
          FPloc(i) -= JxW_mid* dxU.trace()*psi_c[qp][i];
          //FPloc(i) -= JxW_new* dxU_new.trace() *psi_c[qp][i];

        if (sum_vec(SV_c)){
          // residue
          for (int I = 0; I < nodes_per_cell; I++){
            int K = SV_c[I];
            if (K != 0){
              for (int C = 0; C < 3; C++){
                if (C < 2){
                  FZloc(I*3 + C) += JxW_mid*
                      ( rho*(unsteady*dUdt(C) + has_convec*Uconv_qp.dot(dxU.row(C)))*phi_c[qp][I] + // aceleração
                        visc*dxphi_c.row(I).dot(dxU.row(C) + dxU.col(C).transpose()) - //rigidez  //transpose() here is to add 2 row-vectors
                        force_at_mid(C)*phi_c[qp][I] -// força
                        Pqp_new*dxphi_c(I,C) );
                }
                else{
                  Xg = XG[K-1];
                  mesh->getNodePtr(cell->getNodeId(I))->getCoord(XIp.data(),dim);
                  Rotf = SolidVel(XIp,Xg,Zw);  //fixed evaluation point Ip
                  Rotv = SolidVel(Xqp,Xg,Zw);  //variable evaluation point qp
                  //cout << XIp.transpose() << "   " << Xqp.transpose() << "   " << Xg.transpose()
                  //     << "   " << Rotf.transpose() << "   " << Rotv.transpose() << endl;
                  auxRotf << dxU.col(0).dot(Rotf),dxU.col(1).dot(Rotf);
                  auxRotv << dxU.col(0).dot(Rotv),dxU.col(1).dot(Rotv);
                  auxTenf << auxRotf(0)*dxphi_c(I,0),auxRotf(0)*dxphi_c(I,1),
                             auxRotf(1)*dxphi_c(I,0),auxRotf(1)*dxphi_c(I,1);
                  auxTenv << auxRotf(0)*dxphi_c(I,0),auxRotf(0)*dxphi_c(I,1)-phi_c[qp][I],
                             auxRotf(1)*dxphi_c(I,0)+phi_c[qp][I],auxRotf(1)*dxphi_c(I,1);

                  FZloc(I*3 + C) += JxW_mid*(
                      rho*(unsteady*dUdt.dot(auxRotv) + has_convec*Uconv_qp.dot(auxRotv))*phi_c[qp][I] +
                      visc*(DobCont(dxU,auxTenv)+DobCont(dxU.transpose(),auxTenv)) -
                      force_at_mid.dot(auxRotv)*phi_c[qp][I] -
                      Pqp_new*auxTenv.trace() );
                }
              }
            }
          }
          // jacobian
          for (int J = 0; J < nodes_per_cell; J++){
            int K = SV_c[J];
            if (K != 0){
              for (int c = 0; c < dim; c++){
                for (int i = 0; i < n_dofs_u_per_cell/dim; i++){
                  for (int D = 0; D < 3; D++){
                    //Z1loc(i*dim+c,J*3+D) =
                  }
                }
              }
            }
          }

        }  //endif sum_vec





#if(false)
        // ----------------
        //
        //  STABILIZATION
        //
        //  ----------------
        if (behaviors & (BH_bble_condens_PnPn | BH_bble_condens_CR))
        {
          dxbble = invFT_c_mid * dLbble[qp];

          for (int c = 0; c < dim; c++)
          {
            for (int j = 0; j < n_dofs_u_per_cell/dim; j++)
            {
              for (int d = 0; d < dim; d++)
              {
                delta_cd = c==d;

                if (compact_bubble)
                {
                  Bbn(c, j*dim + d) += JxW_mid*
                                       ( utheta*visc*(delta_cd * dxphi_c.row(j).dot(dxbble) + dxphi_c(j,c)*dxbble(d)) ); // rigidez

                  Bnb(j*dim + d, c) += JxW_mid*
                                       ( utheta*visc*(delta_cd * dxphi_c.row(j).dot(dxbble) + dxphi_c(j,c)*dxbble(d)) ); // rigidez
                }
                else
                {
                  Bbn(c, j*dim + d) += JxW_mid*
                                       ( has_convec*bble[qp]*utheta *rho*( delta_cd*Uconv_qp.dot(dxphi_c.row(j)) + dxU(c,d)*phi_c[qp][j] ) // convective
                                       + unsteady*delta_cd*rho*bble[qp]*phi_c[qp][j]/dt // time derivative
                                       + utheta*visc*(delta_cd * dxphi_c.row(j).dot(dxbble) + dxphi_c(j,c)*dxbble(d)) ); // rigidez

                  Bnb(j*dim + d, c) += JxW_mid*
                                       ( has_convec*phi_c[qp][j]*utheta *rho*( delta_cd*Uconv_qp.dot(dxbble) ) // convective
                                       + delta_cd*rho*phi_c[qp][j]*bble[qp]/dt * unsteady // time derivative
                                       + utheta*visc*(delta_cd * dxphi_c.row(j).dot(dxbble) + dxphi_c(j,c)*dxbble(d)) ); // rigidez
                }
              }
            }
            if (behaviors & BH_bble_condens_PnPn)
              for (int j = 0; j < n_dofs_p_per_cell; ++j)
                Gbp(c, j) -= JxW_mid*psi_c[qp][j]*dxbble(c);
          }

          for (int c = 0; c < dim; c++)
          {
            for (int d = 0; d < dim; d++)
            {
              delta_cd = c==d;

              if (compact_bubble)
              {
                iBbb(c, d) += JxW_mid*
                              ( utheta*visc*(delta_cd* dxbble.dot(dxbble) + dxbble(d)*dxbble(c)) ); // rigidez
              }
              else
              {
                iBbb(c, d) += JxW_mid*
                              ( has_convec*bble[qp]*utheta *rho*( delta_cd*Uconv_qp.dot(dxbble) ) // convective
                              + delta_cd*rho*bble[qp]*bble[qp]/dt * unsteady // time derivative
                              + utheta*visc*(delta_cd* dxbble.dot(dxbble) + dxbble(d)*dxbble(c)) ); // rigidez
              }


            }
            if (compact_bubble)
            {
              FUb(c) += JxW_mid*
                        ( visc*dxbble.dot(dxU.row(c) + dxU.col(c).transpose()) - //rigidez
                          Pqp_new*dxbble(c) - // pressão
                          force_at_mid(c)*bble[qp] ); // força
            }
            else
            {
              FUb(c) += JxW_mid*
                        ( bble[qp]*rho*(dUdt(c)*unsteady + has_convec*Uconv_qp.dot(dxU.row(c))) + // time derivative + convective
                          visc*dxbble.dot(dxU.row(c) + dxU.col(c).transpose()) - //rigidez
                          Pqp_new*dxbble(c) - // pressão
                          force_at_mid(c)*bble[qp] ); // força
            }

          }
        }
        else
        if(behaviors & BH_GLS)
        {
          Res = rho*( dUdt * unsteady + has_convec*dxU*Uconv_qp) + dxP_new - force_at_mid;

          for (int j = 0; j < n_dofs_u_per_cell/dim; ++j)
          {
            dResdu = unsteady*(ddt_factor*rho*phi_c[qp][j]/dt)*I + has_convec*rho*utheta*( phi_c[qp][j]*dxU + Uconv_qp.dot(dxphi_c.row(j))*I );

            for (int i = 0; i < n_dofs_p_per_cell; ++i)
            {
              vec = dxpsi_c.row(i).transpose();
              vec = dResdu.transpose()*vec;
              vec = -JxW_mid*tauk* vec;
              for (int d = 0; d < dim; d++)
                Dloc(i, j*dim + d) += vec(d);

              // atençao nos indices
              vec = JxW_mid*tauk* has_convec*rho*Uconv_qp.dot(dxphi_c.row(j))* dxpsi_c.row(i).transpose();
              for (int d = 0; d < dim; d++)
                Cloc(j*dim + d,i) += vec(d);
            }

            for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
            {
              // supg term
              Ten = JxW_mid*tauk* has_convec*( utheta*rho*phi_c[qp][j]*Res*dxphi_c.row(i) + rho*Uconv_qp.dot(dxphi_c.row(i))*dResdu );
              // divergence term
              Ten+= JxW_mid*delk*utheta*dxphi_c.row(i).transpose()*dxphi_c.row(j);

              for (int c = 0; c < dim; ++c)
                for (int d = 0; d < dim; ++d)
                  Aloc(i*dim + c, j*dim + d) += Ten(c,d);
            }
          }

          for (int i = 0; i < n_dofs_p_per_cell; ++i)
            for (int j = 0; j < n_dofs_p_per_cell; ++j)
              Eloc(i,j) -= tauk*JxW_mid * dxphi_c.row(i).dot(dxphi_c.row(j));


          for (int i = 0; i < n_dofs_u_per_cell/dim; i++)
          {
            vec = JxW_mid*( has_convec*tauk* rho* Uconv_qp.dot(dxphi_c.row(i)) * Res + delk*dxU.trace()*dxphi_c.row(i).transpose() );

            for (int c = 0; c < dim; c++)
              FUloc(i*dim + c) += vec(c);

          }
          for (int i = 0; i < n_dofs_p_per_cell; ++i)
            FPloc(i) -= JxW_mid *tauk* dxpsi_c.row(i).dot(Res);
            //FPloc(i) -= JxW_mid *tauk* dxpsi_c.row(i).dot(dxP_new - force_at_mid); // somente laplaciano da pressao
        }

        if (behaviors & BH_bble_condens_CR)
        {
          bble_integ += JxW_mid*bble[qp];

          for (int c = 0; c < dim; ++c)
          {
            for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
              for (int j = 0; j < dim; ++j) // pressure gradient
                Gnx(i*dim + c,j) -= JxW_mid* (Xqp(j) - Xc(j))*dxphi_c(i,c);

            FPx(c) -= JxW_mid* dxU.trace()*(Xqp(c) - Xc(c));
          }
        }

#endif
      } // fim quadratura
#if(false)
      //Dloc += utheta*Gloc.transpose();

      //
      // stabilization
      //
      if ((behaviors & BH_bble_condens_PnPn) && !compact_bubble)
      {
        //iBbb = iBbb.inverse().eval();
        invert(iBbb,dim);

        FUloc = FUloc - Bnb*iBbb*FUb;
        FPloc = FPloc - utheta*Gbp.transpose()*iBbb*FUb;

        Dpb = utheta*Gbp.transpose();

        // correções com os coeficientes da bolha

        Ubqp = -utheta*iBbb*FUb; // U bolha no tempo n+utheta

        for (int qp = 0; qp < n_qpts_cell; ++qp)
        {
          F_c_mid = x_coefs_c_mid_trans * dLqsi_c[qp];
          inverseAndDet(F_c_mid, dim, invF_c_mid,J_mid);
          invFT_c_mid= invF_c_mid.transpose();

          Uqp = u_coefs_c_mid_trans * phi_c[qp]; //n+utheta
          dxbble = invFT_c_mid * dLbble[qp];
          dxUb = Ubqp*dxbble.transpose();

          weight = quadr_cell->weight(qp);
          JxW_mid = J_mid*weight;

          for (int j = 0; j < n_dofs_u_per_cell/dim; ++j)
          {
            for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
            {
              Ten = has_convec*JxW_mid*rho*utheta*phi_c[qp][i]* phi_c[qp][j] * dxUb; // advecção

              for (int c = 0; c < dim; ++c)
                for (int d = 0; d < dim; ++d)
                  Aloc(i*dim + c, j*dim + d) += Ten(c,d);
            }
            Ten = has_convec*JxW_mid* rho*utheta* bble[qp] * phi_c[qp][j] *dxUb; // advecção
            for (int c = 0; c < dim; ++c)
              for (int d = 0; d < dim; ++d)
                Bbn(c, j*dim + d) += Ten(c,d);
          }
        } // fim quadratura 2 vez

        Aloc -= Bnb*iBbb*Bbn;
        Gloc -= Bnb*iBbb*Gbp;
        Dloc -= Dpb*iBbb*Bbn;
        Eloc = -Dpb*iBbb*Gbp;

      }
      if(behaviors & BH_GLS)
      {
        Gloc += Cloc;
      }
      if(behaviors & BH_bble_condens_CR)
      {
        Ubqp.setZero();
        //for (int i = 0; i < Gnx.cols(); ++i)
        // for (int j = 0; j < Gnx.rows(); ++j)
        // Ubqp(i) += Gnx(j,i)*u_coefs_c_new(j);
        //Ubqp /= -bble_integ;
        //Ubqp *= utheta;

        for (int c = 0; c < dim; ++c)
          for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
            for (int j = 0; j < dim; ++j) // pressure gradient
              //Gnx(i*dim + c,j) -= JxW_mid* (Xqp(j) - Xc(j))*dxphi_c(i,c);
            Ubqp(j) += Gnx(i*dim + c,j) * u_coefs_c_mid_trans(c,i);

        Ubqp /= -bble_integ;
        Ubqp *= utheta;


        //Ubqp = -Gnx.transpose()*u_coefs_c_new; // U bolha no tempo n+utheta

        for (int qp = 0; qp < n_qpts_cell; ++qp)
        {
          F_c_mid = x_coefs_c_mid_trans * dLqsi_c[qp];
          inverseAndDet(F_c_mid, dim, invF_c_mid,J_mid);
          invFT_c_mid= invF_c_mid.transpose();

          Uqp = u_coefs_c_mid_trans * phi_c[qp]; //n+utheta
          dxbble = invFT_c_mid * dLbble[qp];
          dxUb = Ubqp*dxbble.transpose();

          weight = quadr_cell->weight(qp);
          JxW_mid = J_mid*weight;

          for (int j = 0; j < n_dofs_u_per_cell/dim; ++j)
          {
            for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
            {
              Ten = has_convec*JxW_mid*rho*utheta*phi_c[qp][i]* phi_c[qp][j] * dxUb; // advecção

              for (int c = 0; c < dim; ++c)
                for (int d = 0; d < dim; ++d)
                  Aloc(i*dim + c, j*dim + d) += Ten(c,d);
            }
            Ten = has_convec*JxW_mid* rho*utheta* bble[qp] * phi_c[qp][j] *dxUb; // advecção
            for (int c = 0; c < dim; ++c)
              for (int d = 0; d < dim; ++d)
                Bbn(c, j*dim + d) += Ten(c,d);
          }
        } // fim quadratura 2 vez

        double const a = 1./(bble_integ*bble_integ);
        double const b = 1./bble_integ;
        Aloc += utheta*a*Gnx*iBbb*Gnx.transpose() - utheta*b*Bnb*Gnx.transpose() - b*Gnx*Bbn;

        FUloc += a*Gnx*iBbb*FPx - b*Bnb*FPx - b*Gnx*FUb;
      }


//cout << "\n" << FUloc << endl; cout << "\n" << Aloc << endl; cout << "\n" << Gloc << endl; cout << "\n" << Dloc << endl;
      // Projection - to force non-penetrarion bc
      mesh->getCellNodesId(&*cell, cell_nodes.data());
      getProjectorMatrix(Prj, nodes_per_cell, cell_nodes.data(), Vec_x_1, current_time+dt, *this);
//cout << Prj << endl;
      FUloc = Prj*FUloc;
      Aloc = Prj*Aloc*Prj;
      Gloc = Prj*Gloc;
      Dloc = Dloc*Prj;
//cout << "\n" << FUloc << endl; cout << "\n" << Aloc << endl; cout << "\n" << Gloc << endl; cout << "\n" << Dloc << endl;
      if (force_pressure)
      {
        for (int i = 0; i < mapP_c.size(); ++i)
        {
          if (mapP_c(i) == null_space_press_dof)
          {
            Gloc.col(i).setZero();
            Dloc.row(i).setZero();
            FPloc(i) = 0;
            Eloc.col(i).setZero();
            Eloc.row(i).setZero();
            break;
          }
        }
      }

#ifdef FEP_HAS_OPENMP
      FEP_PRAGMA_OMP(critical)
#endif
      {
        VecSetValues(Vec_fun_fs, mapU_c.size(), mapU_c.data(), FUloc.data(), ADD_VALUES);
        VecSetValues(Vec_fun_fs, mapP_c.size(), mapP_c.data(), FPloc.data(), ADD_VALUES);
        MatSetValues(*JJ, mapU_c.size(), mapU_c.data(), mapU_c.size(), mapU_c.data(), Aloc.data(),  ADD_VALUES);
        MatSetValues(*JJ, mapU_c.size(), mapU_c.data(), mapP_c.size(), mapP_c.data(), Gloc.data(),  ADD_VALUES);
        MatSetValues(*JJ, mapP_c.size(), mapP_c.data(), mapU_c.size(), mapU_c.data(), Dloc.data(),  ADD_VALUES);
        MatSetValues(*JJ, mapP_c.size(), mapP_c.data(), mapP_c.size(), mapP_c.data(), Eloc.data(),  ADD_VALUES);
      }
#endif  //for debug
    }  //end for cell


  }  //end LOOP NAS CÉLULAS parallel

  //LOOP FOR SOLID-ONLY CONTRIBUTION
  {
    VectorXd   FZsloc(3);
    VectorXi   mapZ_s(3);
    VectorXd   z_coefs_old(3);
    VectorXd   z_coefs_new(3);
    Vector     dZdt(3);
    Vector     Grav;


    for (int K = 0; K < N_Solids; K++){
      Grav = gravity(XG[K]);
      for (int C = 0; C < 3; C++){
        mapZ_s(C) = dof_handler[DH_UNKM].getVariable(VAR_U).numPositiveDofs() - 1
                                                     + 3*(K+1) - 2 + C;
      }
      VecGetValues(Vec_uzp_0,    mapZ_s.size(), mapZ_s.data(), z_coefs_old.data());
      VecGetValues(Vec_uzp_k ,   mapZ_s.size(), mapZ_s.data(), z_coefs_new.data());
      dZdt = (z_coefs_new - z_coefs_old)/dt;
      for (int C = 0; C < 3; C++){
        if (C < 2){
          FZsloc(C) = MV[K]*dZdt(C) - MV[K]*Grav(C);
        }
        else{
          FZsloc(C) = 0.4*MV[K]*RV[K]*RV[K]*dZdt(C);  // + intertia tensor component instead of dZdt
        }
      }
      VecSetValues(Vec_fun_fs, mapZ_s.size(), mapZ_s.data(), FZsloc.data(), ADD_VALUES);
    }
  }
#if (false)
  // LOOP NAS FACES DO CONTORNO (Neumann)
  //~ FEP_PRAGMA_OMP(parallel default(none) shared(Vec_uzp_k,Vec_fun_fs,cout))
  {
    int                 tag;
    bool                is_neumann;
    bool                is_surface;
    bool                is_solid;
    //MatrixXd           u_coefs_c_new(n_dofs_u_per_facet/dim, dim);
    //VectorXd           p_coefs_f(n_dofs_p_per_facet);
    MatrixXd            u_coefs_f_mid_trans(dim, n_dofs_u_per_facet/dim);  // n+utheta
    MatrixXd            u_coefs_f_old(n_dofs_u_per_facet/dim, dim);        // n
    MatrixXd            u_coefs_f_new(n_dofs_u_per_facet/dim, dim);        // n+1
    MatrixXd            u_coefs_f_old_trans(dim,n_dofs_u_per_facet/dim);   // n
    MatrixXd            u_coefs_f_new_trans(dim,n_dofs_u_per_facet/dim);   // n+1

    MatrixXd            x_coefs_f_mid_trans(dim, n_dofs_v_per_facet/dim); // n+utheta
    MatrixXd            x_coefs_f_new(n_dofs_v_per_facet/dim, dim);       // n+1
    MatrixXd            x_coefs_f_new_trans(dim, n_dofs_v_per_facet/dim); // n+1
    MatrixXd            x_coefs_f_old(n_dofs_v_per_facet/dim, dim);       // n
    MatrixXd            x_coefs_f_old_trans(dim, n_dofs_v_per_facet/dim); // n

    MatrixXd            noi_coefs_f_new(n_dofs_v_per_facet/dim, dim);  // normal interpolada em n+1
    MatrixXd            noi_coefs_f_new_trans(dim, n_dofs_v_per_facet/dim);  // normal interpolada em n+1

    Tensor              F_f_mid(dim,dim-1);       // n+utheta
    Tensor              invF_f_mid(dim-1,dim);    // n+utheta
    Tensor              fff_f_mid(dim-1,dim-1);   // n+utheta; fff = first fundamental form
    //Tensor              invFT_c_mid(dim,dim);   // n+utheta

    MatrixXd            Aloc_f(n_dofs_u_per_facet, n_dofs_u_per_facet);
    VectorXd            FUloc(n_dofs_u_per_facet);

    MatrixXd            tmp(n_dofs_u_per_facet,n_dofs_u_per_facet);

    VectorXi            mapU_f(n_dofs_u_per_facet);
    VectorXi            mapP_f(n_dofs_p_per_facet);
    VectorXi            mapM_f(dim*nodes_per_facet);

    MatrixXd            Prj(n_dofs_u_per_facet,n_dofs_u_per_facet);

    MatrixXd            dxphi_f(n_dofs_u_per_facet/dim, dim);
    Tensor              dxU_f(dim,dim);   // grad u
    Vector              Xqp(dim);
    Vector              Xqp2(dim);
    Vector              Xqp_new(dim);
    Vector              Xqp_old(dim);
    Vector              Uqp(dim);
    Vector              Uqp_new(dim);
    Vector              Uqp_old(dim);
    //VectorXd          FUloc(n_dofs_u_per_facet);
    VectorXi            facet_nodes(nodes_per_facet);
    Vector              normal(dim);
    Vector              noi(dim); // normal interpolada
    Vector              some_vec(dim);
    double              J_mid=0,JxW_mid;
    double              weight=0;
    //double              visc;
    //double              rho;
    Vector              Uqp_solid(dim);

    Vector              traction_(dim);

    //~ const int tid = omp_get_thread_num();
    //~ const int nthreads = omp_get_num_threads();
//~
    //~ facet_iterator facet = mesh->facetBegin(tid,nthreads);
    //~ facet_iterator facet_end = mesh->facetEnd(tid,nthreads);

    // LOOP NAS FACES DO CONTORNO
    facet_iterator facet = mesh->facetBegin();
    facet_iterator facet_end = mesh->facetEnd();  // the next if controls the for that follows
    if (neumann_tags.size() != 0 || interface_tags.size() != 0 || solid_tags.size() != 0)
    for (; facet != facet_end; ++facet)
    {
      tag = facet->getTag();
      is_neumann = is_in(tag, neumann_tags);
      is_surface = is_in(tag, interface_tags);
      is_solid   = is_in(tag, solid_tags);

      //if ((!is_neumann))
      if ((!is_neumann) && (!is_surface) && (!is_solid))
      //PetscFunctionReturn(0);
        continue;

      // mapeamento do local para o global:
      //
      dof_handler[DH_UNKS].getVariable(VAR_U).getFacetDofs(mapU_f.data(), &*facet);  cout << mapU_f << endl << endl;  //unk. global ID's
      dof_handler[DH_UNKS].getVariable(VAR_P).getFacetDofs(mapP_f.data(), &*facet);  //cout << mapP_f << endl;
      dof_handler[DH_MESH].getVariable(VAR_M).getFacetDofs(mapM_f.data(), &*facet);  //cout << mapM_f << endl;

      VecGetValues(Vec_normal,  mapM_f.size(), mapM_f.data(), noi_coefs_f_new.data());
      VecGetValues(Vec_x_0,     mapM_f.size(), mapM_f.data(), x_coefs_f_old.data());
      VecGetValues(Vec_x_1,     mapM_f.size(), mapM_f.data(), x_coefs_f_new.data());
      VecGetValues(Vec_up_0,    mapU_f.size(), mapU_f.data(), u_coefs_f_old.data());
      VecGetValues(Vec_uzp_k ,   mapU_f.size(), mapU_f.data(), u_coefs_f_new.data());

      // get nodal coordinates of the old and new cell
      mesh->getFacetNodesId(&*facet, facet_nodes.data());

      x_coefs_f_old_trans = x_coefs_f_old.transpose();
      x_coefs_f_new_trans = x_coefs_f_new.transpose();
      u_coefs_f_old_trans = u_coefs_f_old.transpose();
      u_coefs_f_new_trans = u_coefs_f_new.transpose();
      noi_coefs_f_new_trans = noi_coefs_f_new.transpose();

      u_coefs_f_mid_trans = utheta*u_coefs_f_new_trans + (1.-utheta)*u_coefs_f_old_trans;
      x_coefs_f_mid_trans = utheta*x_coefs_f_new_trans + (1.-utheta)*x_coefs_f_old_trans;

      FUloc.setZero();
      Aloc_f.setZero();

      //visc = muu(tag);
      //rho  = pho(Xqp,tag);

      //noi_coefs_f_new_trans = x_coefs_f_mid_trans;


      for (int qp = 0; qp < n_qpts_facet; ++qp)
      {

        F_f_mid   = x_coefs_f_mid_trans * dLqsi_f[qp];

        if (dim==2)
        {
          normal(0) = +F_f_mid(1,0);
          normal(1) = -F_f_mid(0,0);
          normal.normalize();
        }
        else
        {
          normal = cross(F_f_mid.col(0), F_f_mid.col(1));
          normal.normalize();
        }

        fff_f_mid.resize(dim-1,dim-1);
        fff_f_mid  = F_f_mid.transpose()*F_f_mid;
        J_mid     = sqrt(fff_f_mid.determinant());
        invF_f_mid = fff_f_mid.inverse()*F_f_mid.transpose();

        weight  = quadr_facet->weight(qp);
        JxW_mid = J_mid*weight;
        Xqp     = x_coefs_f_mid_trans * qsi_f[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
        dxphi_f = dLphi_f[qp] * invF_f_mid;
        dxU_f   = u_coefs_f_mid_trans * dxphi_f; // n+utheta
        Uqp     = u_coefs_f_mid_trans * phi_f[qp];
        noi     = noi_coefs_f_new_trans * qsi_f[qp];

        if (is_neumann)
        {
          //Vector no(Xqp);
          //no.normalize();
          //traction_ = utheta*(traction(Xqp,current_time+dt,tag)) + (1.-utheta)*traction(Xqp,current_time,tag);
          traction_ = traction(Xqp, normal, current_time + dt*utheta,tag);
          //traction_ = (traction(Xqp,current_time,tag) +4.*traction(Xqp,current_time+dt/2.,tag) + traction(Xqp,current_time+dt,tag))/6.;

          for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
          {
            for (int c = 0; c < dim; ++c)
            {
              FUloc(i*dim + c) -= JxW_mid * traction_(c) * phi_f[qp][i] ; // força
            }
          }
        }

        if (is_surface)
        {
          //Vector no(Xqp);
          //no.normalize();
          for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
          {
            for (int c = 0; c < dim; ++c)
            {
//              FUloc(i*dim + c) += JxW_mid *gama(Xqp,current_time,tag)*(dxphi_f(i,c) + (unsteady*dt) *dxU_f.row(c).dot(dxphi_f.row(i))); // correto
              FUloc(i*dim + c) += JxW_mid *gama(Xqp,current_time,tag)*dxphi_f(i,c); //inicialmente descomentado
              //FUloc(i*dim + c) += JxW_mid *gama(Xqp,current_time,tag)*normal(c)* phi_f[qp][i];
              //for (int d = 0; d < dim; ++d)
              //  FUloc(i*dim + c) += JxW_mid * gama(Xqp,current_time,tag)* ( (c==d?1:0) - noi(c)*noi(d) )* dxphi_f(i,d) ;
              //FUloc(i*dim + c) += JxW_mid * gama(Xqp,current_time,tag)* ( unsteady*dt *dxU_f.row(c).dot(dxphi_f.row(i)));
            }
          }

          if (false) // semi-implicit term //inicialmente false
          {
            for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
              for (int j = 0; j < n_dofs_u_per_facet/dim; ++j)
                for (int c = 0; c < dim; ++c)
                  Aloc_f(i*dim + c, j*dim + c) += utheta*JxW_mid* (unsteady*dt) *gama(Xqp,current_time,tag)*dxphi_f.row(i).dot(dxphi_f.row(j));
          }

        }

        if (is_solid)
        {
          Uqp_solid = solid_veloc(Xqp, current_time+utheta*dt, tag);

          for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
          {
            for (int c = 0; c < dim; ++c)
            {
              FUloc(i*dim + c) += JxW_mid *beta_diss()*(Uqp(c)-Uqp_solid(c))*phi_f[qp][i];
              //FUloc(i*dim + c) += x_coefs_f_old_trans.norm()*beta_diss()*Uqp(c)*phi_f[qp][i];

              for (int j = 0; j < n_dofs_u_per_facet/dim; ++j)
                  Aloc_f(i*dim + c, j*dim + c) += utheta*JxW_mid *beta_diss()*phi_f[qp][j]*phi_f[qp][i];

            }
          }
        }

      } // end quadratura


      // Projection - to force non-penetrarion bc
      mesh->getFacetNodesId(&*facet, facet_nodes.data());
      getProjectorMatrix(Prj, nodes_per_facet, facet_nodes.data(), Vec_x_1, current_time+dt, *this);

      FUloc = Prj*FUloc;
      Aloc_f = Prj*Aloc_f*Prj;

      //~ FEP_PRAGMA_OMP(critical)
      {
        VecSetValues(Vec_fun_fs, mapU_f.size(), mapU_f.data(), FUloc.data(), ADD_VALUES);
        MatSetValues(*JJ, mapU_f.size(), mapU_f.data(), mapU_f.size(), mapU_f.data(), Aloc_f.data(),  ADD_VALUES);
      }

    }  //end for facet


  } // end LOOP NAS FACES DO CONTORNO (Neumann)
#endif

  // LINHA DE CONTATO
  //FEP_PRAGMA_OMP(parallel shared(Vec_uzp_k,Vec_fun_fs,cout) default(none))
  {
    // will be useful I hope
    //Real const eps = std::numeric_limits<Real>::epsilon();
    //Real const eps_root = pow(eps,1./3.);
    //double            h;
    //volatile double   hh;

    int              tag;
    bool             is_triple;

    VectorXd         FUloc(n_dofs_u_per_corner);
    MatrixXd         Aloc_r(n_dofs_u_per_corner, n_dofs_u_per_corner);

    VectorXi         mapU_r(n_dofs_u_per_corner);
    VectorXi         mapP_r(n_dofs_p_per_corner);

    MatrixXd         Prj(n_dofs_u_per_corner,n_dofs_u_per_corner);
    VectorXi         corner_nodes(nodes_per_corner);

    bool                gen_error = false;
    //MatrixXd             u_coefs_r_mid(n_dofs_u_per_corner/dim, dim);
    MatrixXd            u_coefs_r_mid_trans(dim, n_dofs_u_per_corner/dim);  // n+utheta
    MatrixXd            u_coefs_r_old(n_dofs_u_per_corner/dim, dim);        // n
    MatrixXd            u_coefs_r_new(n_dofs_u_per_corner/dim, dim);        // n+1
    MatrixXd            x_coefs_r_mid_trans(dim, nodes_per_corner);
    MatrixXd            x_coefs_r_new(nodes_per_corner, dim);
    MatrixXd            x_coefs_r_old(nodes_per_corner, dim);
    Tensor              F_r_mid(dim,dim-2);
    Tensor              invF_r_mid(dim-2,dim);
    MatrixXd            dxphi_r(n_dofs_u_per_corner/dim, dim);
    Tensor              dxU_r(dim,dim);   // grad u
    Vector              Xqp(dim);
    Vector              Uqp(dim);
    //VectorXd          FUloc(n_dofs_u_per_corner);
    //MatrixXd         Aloc_r(n_dofs_u_per_corner, n_dofs_u_per_corner);
    Vector              normal(dim);
    Vector              line_normal(dim);
    Vector              solid_point(dim); // ponto na superfície do sólido .. ele é único
    Vector              point_a(dim); // ponto na linha triplice
    Vector              point_b(dim); // ponto na linha triplice
    Vector              ifacet_normal(dim); // ponto na linha triplice
    double              line_normal_sign = 0; // +1 or -1
    double              J_mid=0, JxW_mid;
    double              weight=0;
    double              gama_mid;
    //double              visc;
    //double              rho;
    int                 iCs[FEPIC_MAX_ICELLS];
    int                 eiCs[FEPIC_MAX_ICELLS];
    int                 *iCs_end;
    int                 *iCs_it;
    Cell                *fluid_cell;

    //VectorXi            mapU_r(n_dofs_u_per_corner);
    //VectorXi            mapP_r(n_dofs_p_per_corner);
    VectorXi            mapM_r(dim*nodes_per_corner);

    //const int tid = omp_get_thread_num();
    //const int nthreads = omp_get_num_threads();
    //const int n_corner_colors = mesh->numCornerColors();


    // LOOP NAS ARESTAS DA LINHA TRIPLICE
    CellElement * corner;

    if (triple_tags.size() != 0)
    for (int _r = 0; _r < n_corners_total; ++_r)
    {
      if (dim==2)
      {
        corner = mesh->getNodePtr(_r);
        if (!mesh->isVertex(corner))
          continue;
      }
      else
        corner = mesh->getCornerPtr(_r);
      if (corner->isDisabled())
        continue;

      tag = corner->getTag();
      is_triple = is_in(tag,triple_tags);
      if (!is_triple)
        continue;

      FUloc.setZero();
      Aloc_r.setZero();

      mesh->getCornerNodesId(&*corner, corner_nodes.data());
     // mesh->getNodesCoords(corner_nodes.begin(), corner_nodes.end(), x_coefs_r_mid.data());

      dof_handler[DH_UNKS].getVariable(VAR_U).getCornerDofs(mapU_r.data(), &*corner);
      dof_handler[DH_UNKS].getVariable(VAR_P).getCornerDofs(mapP_r.data(), &*corner);
      dof_handler[DH_MESH].getVariable(VAR_M).getCornerDofs(mapM_r.data(), &*corner);

      VecGetValues(Vec_x_0,     mapM_r.size(), mapM_r.data(), x_coefs_r_old.data());
      VecGetValues(Vec_x_1,     mapM_r.size(), mapM_r.data(), x_coefs_r_new.data());
      VecGetValues(Vec_up_0,    mapU_r.size(), mapU_r.data(), u_coefs_r_old.data());
      VecGetValues(Vec_uzp_k,    mapU_r.size(), mapU_r.data(), u_coefs_r_new.data());

      u_coefs_r_mid_trans = utheta*u_coefs_r_new.transpose() + (1.-utheta)*u_coefs_r_old.transpose();
      x_coefs_r_mid_trans = utheta*x_coefs_r_new.transpose() + (1.-utheta)*x_coefs_r_old.transpose();

      //visc = muu(tag);
      //rho  = pho(Xqp,tag);

      if (dim==3)
      {
        // encontrando o ponto da superfície sólido. com ele, é calculado uma normal
        // que corrigi o sinal de line_normal
        iCs_end = mesh->edgeStar(corner, iCs, eiCs);
        if (iCs_end == iCs)
        {
          printf("ERROR!: no icell found\n");
          throw;
        }
        gen_error = true;
        for (iCs_it = iCs; iCs_it != iCs_end ; ++iCs_it)
        {
          fluid_cell = mesh->getCellPtr(*iCs_it);

          for (int kk = 0; kk < mesh->numVerticesPerCell(); ++kk)
          {
            int const nodekk_id = fluid_cell->getNodeId(kk);
            Point const* pp      = mesh->getNodePtr(fluid_cell->getNodeId(kk) );
            const int    tag_aux = pp->getTag();

            if ((is_in(tag_aux, solid_tags) || is_in(tag_aux,feature_tags)) && !is_in(tag_aux, triple_tags))
            {
              if (corner_nodes(0) != nodekk_id && corner_nodes(1) != nodekk_id)
              {
                gen_error = false;
                pp->getCoord(solid_point.data(), dim);
                break;
              }
            }
          }
          if (gen_error==false)
            break;

        }
        if (gen_error)
        {
          printf("ERROR!: solid point not found\n");
          cout << "corner id: " << (mesh->getCellPtr(corner->getIncidCell())->getCornerId(corner->getPosition())) << endl;
          cout << "first icell : " << (*iCs) << endl;
          mesh->getNodePtr(corner_nodes(0))->getCoord(point_a.data(),dim);
          mesh->getNodePtr(corner_nodes(1))->getCoord(point_b.data(),dim);
          cout << "point a = " << point_a[0] << " " << point_a[1] << " " << point_a[2] << "\n";
          cout << "point b = " << point_b[0] << " " << point_b[1] << " " << point_b[2] << "\n";
          cout << "number of cells = " << (iCs_end-iCs) << endl;
          throw;
        }


        mesh->getNodePtr(corner_nodes(0))->getCoord(point_a.data(),dim);
        mesh->getNodePtr(corner_nodes(1))->getCoord(point_b.data(),dim);

        // se (a-c) cross (b-c) dot solid_normal > 0, então line_normal_sign = 1, se não, =-1
        Xqp = (point_a+point_b)/2.;
        point_a -= solid_point;
        point_b -= solid_point;
        normal  = solid_normal(Xqp, current_time, tag);
        line_normal_sign = point_a(0)*point_b(1)*normal(2) + point_a(1)*point_b(2)*normal(0) + point_a(2)*point_b(0)*normal(1)
                          -point_a(0)*point_b(2)*normal(1) - point_a(1)*point_b(0)*normal(2) - point_a(2)*point_b(1)*normal(0);
        if (line_normal_sign>0)
          line_normal_sign = 1;
        else
          line_normal_sign = -1;


      }
      else // dim==2
      {
        Point * point = mesh->getNodePtr(corner_nodes[0]);
        Point * sol_point = NULL;
        Point * sol_point_2;
        int iVs[FEPIC_MAX_ICELLS];
        int *iVs_end, *iVs_it;
        Vector aux(dim);

        iVs_end = mesh->connectedVtcs(point, iVs);

        // se esse nó está na linha, então existe um vértice vizinho que está no sólido
        for (iVs_it = iVs; iVs_it != iVs_end ; ++iVs_it)
        {
          sol_point_2 = mesh->getNodePtr(*iVs_it);
          if ( is_in(sol_point_2->getTag(), solid_tags) )
            if( !sol_point || (sol_point_2->getTag() > sol_point->getTag()) )
              sol_point = sol_point_2;
        }
        if (!sol_point)
        {
          //FEP_PRAGMA_OMP(critical)
          {
            printf("ERRO: ponto na linha tríplice não tem um vértice vizinho no sólido");
            throw;
          }
        }
        point->getCoord(Xqp.data(),dim);

        normal = solid_normal(Xqp, current_time, tag);

        sol_point->getCoord(aux.data(),dim);


        // choose a line_normal candidate
        line_normal(0) = -normal(1);
        line_normal(1) =  normal(0);

        // check orientation
        if (line_normal.dot(Xqp-aux) < 0)
          line_normal *= -1;
      }


      for (int qp = 0; qp < n_qpts_corner; ++qp)
      {
        if (dim==3)
        {
          F_r_mid   = x_coefs_r_mid_trans * dLqsi_r[qp];
          J_mid = F_r_mid.norm();
          weight  = quadr_corner->weight(qp);
          JxW_mid = J_mid*weight;
        }
        else
        {
          J_mid = 1;
          weight  = 1;
          JxW_mid = 1;
        }
        //invF_r_mid = F_r_mid.transpose()/(J_mid*J_mid);


        Xqp = x_coefs_r_mid_trans * qsi_r[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
        //dxphi_r = dLphi_r[qp] * invF_r_mid;
        //dxU_r   = u_coefs_r_mid_trans * dxphi_r; // n+utheta
        Uqp  = u_coefs_r_mid_trans * phi_r[qp];

        gama_mid = gama(Xqp, current_time+dt*utheta, tag);

        if (dim==3)
        {
          normal  = solid_normal(Xqp, current_time, tag);
          line_normal(0)= F_r_mid(0,0);
          line_normal(1)= F_r_mid(1,0);
          line_normal(2)= F_r_mid(2,0);
          line_normal *= line_normal_sign;
          line_normal = cross(line_normal, normal);
          line_normal.normalize();
        }

        double const Uqp_norm = abs(line_normal.dot(Uqp));
        for (int i = 0; i < n_dofs_u_per_corner/dim; ++i)
        {
          for (int c = 0; c < dim; ++c)
          {
            FUloc(i*dim + c) += JxW_mid*(-gama_mid*cos_theta0() + zeta(Uqp_norm,0)*line_normal.dot(Uqp))*line_normal(c)*phi_r[qp][i];

            //FUloc(i*dim + c) += JxW_mid*(-gama_mid*cos_theta0() )*line_normal(c)*phi_r[qp][i];
            //FUloc(i*dim + c) += JxW_mid*zeta(Uqp_norm,0)*Uqp(c)*phi_r[qp][i];

            for (int j = 0; j < n_dofs_u_per_corner/dim; ++j)
              for (int d = 0; d < dim; ++d)
                Aloc_r(i*dim + c, j*dim + d) += JxW_mid* zeta(Uqp_norm,0)*utheta*line_normal(d) *phi_r[qp][j]*line_normal(c)*phi_r[qp][i];
                //if (c==d)
                //  Aloc_r(i*dim + c, j*dim + d) += JxW_mid* zeta(Uqp_norm,0) *utheta*phi_r[qp][j]*phi_r[qp][i];

          }

          //FUloc.segment(i*dim, dim) += JxW_mid* phi_r[qp][i] * (-gama_mid*cos_theta0() + zeta(0,0)*line_normal.dot(Uqp))*line_normal;
        }


      } // end quadratura




      // Projection - to force non-penetrarion bc
      mesh->getCornerNodesId(&*corner, corner_nodes.data());
      getProjectorMatrix(Prj, nodes_per_corner, corner_nodes.data(), Vec_x_1, current_time+dt, *this);

      FUloc = Prj*FUloc;
      Aloc_r = Prj*Aloc_r*Prj;

      VecSetValues(Vec_fun_fs, mapU_r.size(), mapU_r.data(), FUloc.data(), ADD_VALUES);
      MatSetValues(*JJ, mapU_r.size(), mapU_r.data(), mapU_r.size(), mapU_r.data(), Aloc_r.data(),  ADD_VALUES);
      //cout << FUloc.transpose() << endl;

    }



  }  //end LINHA DE CONTATO

  // boundary conditions on global Jacobian
    // solid & triple tags .. force normal
  if (force_dirichlet)
  {
    int      nodeid;
    int      u_dofs[dim];
    Vector   normal(dim);
    Tensor   A(dim,dim);
    Tensor   I(Tensor::Identity(dim,dim));
    int      tag;

    point_iterator point = mesh->pointBegin();
    point_iterator point_end = mesh->pointEnd();
    for ( ; point != point_end; ++point)
    {
      tag = point->getTag();
      if (!(is_in(tag,feature_tags)   ||
            is_in(tag,solid_tags)     ||
            is_in(tag,interface_tags) ||
            is_in(tag,triple_tags)    ||
            is_in(tag,dirichlet_tags) ||
            is_in(tag,neumann_tags)   ||
            is_in(tag,periodic_tags)     )  )
        continue;
      //dof_handler[DH_UNKS].getVariable(VAR_U).getVertexAssociatedDofs(u_dofs, &*point);
      getNodeDofs(&*point, DH_UNKS, VAR_U, u_dofs);

      nodeid = mesh->getPointId(&*point);
      getProjectorMatrix(A, 1, &nodeid, Vec_x_1, current_time+dt, *this);
      A = I - A;
      MatSetValues(*JJ, dim, u_dofs, dim, u_dofs, A.data(), ADD_VALUES);
    }
  }  //end force_dirichlet

  if (force_pressure)
  {
    double const p =1.0;
    MatSetValues(*JJ, 1, &null_space_press_dof, 1, &null_space_press_dof, &p, ADD_VALUES);
  }

  if(print_to_matlab)
  {
    static bool ja_foi=false;
    if (!ja_foi)
    {
      View(Vec_fun_fs, "rhs.m","res");
      View(*JJ,"jacob.m","Jac");
    }
    ja_foi = true;

  }

  Assembly(Vec_fun_fs);
  Assembly(*JJ);

  PetscFunctionReturn(0);

} // END formFunction


PetscErrorCode AppCtx::formJacobian_fs(SNES snes_fs,Vec Vec_uzp_k, Mat* /*Mat_Jac*/, Mat* /*prejac*/, MatStructure * /*flag*/)
{
  PetscBool          found = PETSC_FALSE;
  char               snes_type[PETSC_MAX_PATH_LEN];

  PetscOptionsGetString(PETSC_NULL,"-snes_type",snes_type,PETSC_MAX_PATH_LEN-1,&found);

  if (found)
    if (string(snes_type) == string("test"))
    {
      cout << "WARNING: TESTING JACOBIAN !!!!! \n";
      this->formFunction_fs(snes_fs, Vec_uzp_k, Vec_res_fs);
    }

  PetscFunctionReturn(0);
}
