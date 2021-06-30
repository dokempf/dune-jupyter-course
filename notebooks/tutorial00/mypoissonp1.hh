/******************************************************/
/** a local operator for solving the linear convection-
 *  diffusion equation with standard FEM
 *
 * \f{align*}{
 *   \Delta u &=& f \mbox{ in } \Omega,  \\
 *          u &=& g \mbox{ on } \partial\Omega \\
 * \f}
 * Beware of line number changes, they may corrupt docu!
 */
/******************************************************/
template<typename F, typename FiniteElementMap>
class PoissonP1 :
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
private:
  // define useful types
  typedef typename FiniteElementMap::Traits::FiniteElementType
     FiniteElementType;
  typedef typename FiniteElementType::Traits::LocalBasisType
     LocalBasisType;
  typedef typename LocalBasisType::Traits::DomainType
     DomainType;
  typedef typename LocalBasisType::Traits::RangeFieldType
     RF;
  typedef typename LocalBasisType::Traits::RangeType
     RangeType;
  typedef typename LocalBasisType::Traits::JacobianType
     JacobianType;

  // data members
  enum {dim=LocalBasisType::Traits::dimDomain};
  enum {n=dim+1};
  const F f;              // right hand side function
  DomainType qp;          // center of mass of refelem
  double weight;          // quadrature weight on refelem
  double phihat[n];       // basis functions at qp
  double gradhat[dim][n]; // coordinate x #basisfct

public:
  // define flags controlling global assembler
  enum { doPatternVolume = true };
  enum { doAlphaVolume = true };
  enum { doLambdaVolume = true };

  // Constructor precomputes element independent data
  PoissonP1 (const F& f_, const FiniteElementType& fel)
    : f(f_)
  {
    // select quadrature rule
    Dune::GeometryType gt = fel.type();
    const Dune::QuadratureRule<RF,dim>&
      rule = Dune::QuadratureRules<RF,dim>::rule(gt,1);
    if (rule.size()>1) {
      std::cout << "Wrong quadrature rule!" << std::endl;
      exit(1);
    }

    // position and weight of the quadrature point
    weight = rule[0].weight();
    qp = rule[0].position();

    // check size of the basis
    if (fel.localBasis().size()!=n) {
      std::cout << "Wrong basis!" << std::endl;
      exit(1);
    }

    // evaluate basis functions on refelem
    std::vector<RangeType> phi(n);
    fel.localBasis().evaluateFunction(qp,phi);
    for (int i=0; i<n; i++) phihat[i] = phi[i];

    // evaluate gradients of basis functions on refelem
    std::vector<JacobianType> js(n);
    fel.localBasis().evaluateJacobian(qp,js);
    for (int i=0; i<n; i++)
      for (int j=0; j<dim; j++)
        gradhat[j][i] = js[i][0][j];
  }

  // volume integral depending only on test functions
  template<typename EG, typename LFSV, typename R>
  void lambda_volume (const EG& eg, const LFSV& lfsv,
                      R& r) const
  {
    typename F::Traits::RangeType fval;
    f.evaluate(eg.entity(),qp,fval);
    RF factor=fval*weight*eg.geometry().integrationElement(qp);
    for (int i=0; i<n; i++)
      r.accumulate(lfsv,i,-factor*phihat[i]);
  }

  // jacobian of volume term
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename M>
  void jacobian_volume (const EG& eg, const LFSU& lfsu,
                        const X& x, const LFSV& lfsv,
                        M& mat) const
  {
    // get Jacobian and determinant
    // assume the transformation is linear
    const auto geo = eg.geometry();
    const auto S = geo.jacobianInverseTransposed(qp);
    RF factor = weight*geo.integrationElement(qp);

    // compute gradients of basis functions in transformed element
    double grad[dim][n] = {{0.0}}; // coordinate x #basisfct
    for (int i=0; i<dim; i++) // rows of S
      for (int k=0; k<dim; k++) // columns of S
        for (int j=0; j<n; j++) // columns of gradhat
          grad[i][j] += S[i][k] * gradhat[k][j];

    // compute grad^T * grad
    double A[n][n] = {{0.0}};
    for (int i=0; i<n; i++)
      for (int k=0; k<dim; k++)
        for (int j=0; j<n; j++)
          A[i][j] += grad[k][i]*grad[k][j];

    // store in result
    for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
        mat.accumulate(lfsu,i,lfsu,j,A[i][j]*factor);
  }

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu,
                     const X& x, const LFSV& lfsv,
                     R& r) const
  {
    // get Jacobian and determinant
    // assume the transformation is linear
    const auto geo = eg.geometry();
    const auto S = geo.jacobianInverseTransposed(qp);
    RF factor = weight*geo.integrationElement(qp);

    // compute gradients of basis functions in transformed element
    double grad[dim][n] = {{0.0}};  // coordinate x #basisfct
    for (int i=0; i<dim; i++) // rows of S
      for (int k=0; k<dim; k++) // columns of S
        for (int j=0; j<n; j++) // columns of gradhat
          grad[i][j] += S[i][k] * gradhat[k][j];

    // extract coefficients
    double z_T[n];
    for (int j=0; j<n; j++) z_T[j] = x(lfsu,j); // read coeffs

    // compute gradient u_h
    double graduh[dim] = {0.0};
    for (int k=0; k<dim; k++) // rows of grad
      for (int j=0; j<n; j++) // columns of grad
        graduh[k] += grad[k][j]*z_T[j];

    // scalar products
    double a_T[n] = {0.0};
    for (int k=0; k<dim; k++) // rows of grad
      for (int j=0; j<n; j++)
        a_T[j] += grad[k][j]*graduh[k];

    // store in result
    for (int i=0; i<n; i++)
      r.accumulate(lfsv,i,a_T[i]*factor);
  }

  //! apply local jacobian of the volume term
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& z, const LFSV& lfsv,
                              R& r) const
  {
    alpha_volume(eg,lfsu,z,lfsv,r);
  }

  //! apply local jacobian of the volume term
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const X& z, const LFSV& lfsv,
                              R& r) const
  {
    alpha_volume(eg,lfsu,z,lfsv,r);
  }
};
