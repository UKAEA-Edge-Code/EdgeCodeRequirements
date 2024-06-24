alpha=2.0
m=10.0
epsilon=0.01

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 50
  xmin = 0.0
  xmax = 1.0
  ymin = 0.0
  ymax = 1.0
  elem_type=TRI6
[]

[Variables]
  [./u]
    order=SECOND 
  [../]
[]
[AuxVariables]
  [./err]
  order = CONSTANT
  family = MONOMIAL
  [../]
  [./u_exact]
    order=SECOND 
  [../]
[]
[Functions]
    [./func]
        type = RHSFunc
        alpha = ${alpha}
        m = ${m}
        epsilon = ${epsilon}
    [../]
    [./solution]
        type = SolutionFunc
        alpha = ${alpha}
        m = ${m}
        epsilon = ${epsilon}
    [../]
[]
[AuxKernels]
[./l2_error_aux]
  type = ElementL2ErrorFunctionAux
  variable = err
  function = solution
  coupled_variable = u
[../]
[./exact]
    type = FunctionAux
    function = solution
    variable= u_exact
[../]
[]
[Materials]
  [./D]
    type = AnisoDiffMaterial
    tensor_name =  D
    alpha = ${alpha}
    m = ${m}
    epsilon = ${epsilon}
  [../]
[]
[Kernels]
  [./diff]
    type = MatAnisoDiffusion
    diffusivity = D
    variable = u
  [../]
  [./f]
    type = BodyForce
    function = func
    variable = u
  [../]
[]

[BCs]
    [./top_bottom]
    type = DirichletBC
    variable = u
    boundary = 'top bottom'
    value = 0
    [../]
    [./left_right]
    type = NeumannBC
    variable = u
    boundary = 'right left'
    value = 0
    [../]
[]
[Executioner]
  type = Steady
[]
[Postprocessors]
  [./error]
    type = ElementL2Error
    function = solution
    variable = u
  [../]
[]
[Outputs]
  exodus = true
[]
