#include "Solution.h"

registerMooseObject("egretApp", SolutionFunc);

InputParameters
SolutionFunc::validParams()
{
  auto params = Function::validParams();
  params.addRequiredParam<Real>("alpha", "Parameter");
  params.addRequiredParam<Real>("m", "Parameter");
  params.addRequiredParam<Real>("epsilon", "Parameter");
  return params;
}

SolutionFunc::SolutionFunc(InputParameters const & params)
  : Function(params),
    _alpha(getParam<Real>("alpha")),
    _m(getParam<Real>("m")),
    _epsilon(getParam<Real>("epsilon"))
{
}

Real
SolutionFunc::value(Real t, const Point & p) const
{
    auto x = p(0);
    auto y = p(1);
    return sin(M_PI * y + _alpha * (y*y - y)* cos(_m*M_PI*x))
            + _epsilon * cos(2*M_PI*x)*sin(M_PI * y);
}
