template<typename Number>
class Problem
{
  Number eta;
  Number t;

public:
  typedef Number value_type;

  //! Constructor without arg sets nonlinear term to zero
  Problem () : eta(0.0), t(0.0) {}

  //! Constructor takes eta parameter
  Problem (const Number& eta_) : eta(eta_), t(0.0) {}

  //! nonlinearity
  Number q (Number u) const
  {
    return eta*u*u;
  }

  //! derivative of nonlinearity
  Number qprime (Number u) const
  {
    return 2*eta*u;
  }

  //! right hand side
  template<typename E, typename X>
  Number f (const E& e, const X& x) const
  {
    return 0.0;
  }

  //! boundary condition type function (true = Dirichlet)
  template<typename I, typename X>
  bool b (const I& i, const X& x) const
  {
    auto global = i.geometry().global(x);
    return (global[0]<=1e-7) ? true : false;
  }

  //! Dirichlet extension
  template<typename E, typename X>
  Number g (const E& e, const X& x) const
  {
    auto global = e.geometry().global(x);
    Number s=sin(2.0*M_PI*t);
    for (std::size_t i=1; i<global.size(); i++)
      s*=sin(global[i]*M_PI)*sin(global[i]*M_PI);
    for (std::size_t i=1; i<global.size(); i++)
      s*=sin(10*global[i]*M_PI)*sin(10*global[i]*M_PI);
    return s;
  }

  //! Neumann boundary condition
  template<typename I, typename X>
  Number j (const I& i, const X& x) const
  {
    return 0.0;
  }

  //! Set time in instationary case
  void setTime (Number t_)
  {
    t = t_;
  }
};
