#include "AnisoDiffMaterial.h"


registerMooseObject("egretApp", AnisoDiffMaterial);
registerMooseObject("egretApp", SimpleDiffMaterial);


InputParameters
SimpleDiffMaterial::validParams()
{
  auto params = Material::validParams();
  params.addClassDescription("Material object for defining anisotropic diffusion tensor.");
  params.addRequiredParam<MaterialPropertyName>(
      "tensor_name", "Name of the tensor material property to be created");
  return params;
}

SimpleDiffMaterial::SimpleDiffMaterial(InputParameters const & params)
  : Material(params)
  ,  _A(declareProperty<RealTensorValue>(getParam<MaterialPropertyName>("tensor_name")))
{
}

void
SimpleDiffMaterial::computeQpProperties() 
{

   _A[_qp](0,0) = 1.0;
   _A[_qp](0,1) = 0.0;
   _A[_qp](1,0) = 0.0;
   _A[_qp](1,1) = 1.0;
}

InputParameters
AnisoDiffMaterial::validParams()
{
  auto params = Material::validParams();
  params.addClassDescription("Material object for defining anisotropic diffusion tensor.");
  params.addRequiredParam<Real>("alpha","Parameter");
  params.addRequiredParam<Real>("m","Parameter");
  params.addRequiredParam<Real>("epsilon","Parameter");
  params.addRequiredParam<MaterialPropertyName>(
      "tensor_name", "Name of the tensor material property to be created");
  return params;
}

AnisoDiffMaterial::AnisoDiffMaterial(InputParameters const & params)
  : Material(params)
  ,  _A(declareProperty<RealTensorValue>(getParam<MaterialPropertyName>("tensor_name")))
    ,_alpha(getParam<Real>("alpha"))
    ,_m(getParam<Real>("m"))
    ,_epsilon(getParam<Real>("epsilon"))

{
}

void
AnisoDiffMaterial::computeQpProperties() 
{
   auto x = _q_point[_qp](0);
   auto y = _q_point[_qp](1);

   auto b0 = _alpha * (2 * y - 1) * cos(_m*M_PI*x) + M_PI; 
   auto b1 = M_PI * _alpha * _m * (y*y - y) * sin(_m*M_PI*x);
   auto mod = 1.0/(b0*b0 + b1*b1); 


   _A[_qp](0,0) = b0 * b0 * mod + _epsilon*(1.0-b0*b0*mod);  
   _A[_qp](0,1) = b0 * b1 * mod - _epsilon*b0*b1*mod;  
   _A[_qp](1,0) = _A[_qp](0,1);
   _A[_qp](1,1) = b1 * b1 * mod + _epsilon*(1.0-b1*b1*mod);  
}
