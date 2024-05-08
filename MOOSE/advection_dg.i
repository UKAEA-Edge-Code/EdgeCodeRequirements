#Define some variables outside of any block can use later with ${} syntax
sigma=2.0
vx = -1.0
vy = 0.0

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 128
  ny = 16
  xmin = -60
  xmax = 20
  ymin = -5
  ymax = 5
  elem_type=QUAD8
[]

[Variables]
  [u]
    order = SECOND
    family = L2_LAGRANGE
  []
[]

[Functions]
    #Function of initial condition (can be defined in Code)
    [./gaussian]
    type = ParsedFunction
    expression = 'exp(-(x*x + y*y)/(${sigma} * ${sigma}))'
    [../]
    #Exact solution
    [./gaussian_wt]
    type = ParsedFunction
    expression = 'exp(-((x-${vx}*t)*(x-${vx}*t) + (y-${vy}*t)*(y-${vy}*t))/(${sigma} * ${sigma}))'
    [../]
[]

[Kernels]
  [time_u]
    type = TimeDerivative
    variable = u
  []
  [adv_u]
    #Note upwinding not supported for DG
    type = ConservativeAdvection
    variable = u
    velocity = '${vx} ${vy} 0'
  []
[]
[DGKernels]
  [dg_advection_u]
    type = DGConvection
    variable = u
    velocity = '${vx} ${vy} 0'
  []
[]

#add exact solution to output
[AuxVariables]
  [exact_u]
    order = SECOND
    family = L2_LAGRANGE
  []
[]
[AuxKernels]
    [calc_exact_u]
       type = FunctionAux
       variable = exact_u 
       function = gaussian_wt
    []
[]

#Inital conditions
[ICs]
  [u_ic]
    type = FunctionIC
    variable = u
    function = gaussian
  []
[]

[Executioner]
  type = Transient
  [TimeIntegrator]
  #Second order diagonally implicit Runge Kutta method (Dirk) with two stages
  # https://mooseframework.inl.gov/syntax/Executioner/TimeIntegrator/index.html
  # for list of alternatives
    type = LStableDirk2
  []
  solve_type = 'LINEAR'
  end_time=40
  dt = 1e-2
[]

#Calc the error
[Postprocessors]
 [./lt_err]
    type = ElementL2Error
    variable = u
    function = gaussian_wt
[../]
[]

[Outputs]
  [out]
    type = Exodus
    interval=100
  []
  [csv]
    type=CSV
    interval=1000
  []
[]
