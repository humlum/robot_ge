[comment]: # ( <img src="https://render.githubusercontent.com/render/math?math= )
[comment]: # ( "> )

# robot_ge
This Matlab package implements the fixed-point algorithm for solving the dynamic general equilibrium of [Humlum (2020)](https://andershumlum.com/s/humlumJMP.pdf).

## Main Scripts (`~/code/`)

### `solver.m` 
Load input parameter structures `env`, `par`, `init`, `exp`, `sol`, and set counterfactual experiment ``expNo``. Iterate on the equilibrium labor demand mapping until convergence in the path of wages. For a description of the algorithm, see [Online Appendix F.3](https://andershumlum.com/s/humlumJMP.pdf).

### `out_robot.m`
Output figures and tables for [Section 6.2](https://andershumlum.com/s/humlumJMP.pdf#page=29) "The Distributional Impact of Industrial Robots".

### `out_tax.m`
Output figures and tables for [Section 6.3](https://andershumlum.com/s/humlumJMP.pdf#page=29) "Policy Counterfactuals: The Incidence of Robot Taxes".


## Functions (`~/code/auxiliary`)
The functions documented below are located in the scripts ``solve.m`` and ``simulate.m``. The script ``func.m`` contains smaller auxiliary model functions. The script ``output.m`` contains auxiliary functions for the counterfactual analysis in Sections 6.2 and 6.3.

### `solve.geq`
Simulate firm and worker states forward given a path of wages. Given firm and worker states, solve for wage path that clears the static labor demand period-by-period. See [Online Appendix F.3](https://andershumlum.com/s/humlumJMP.pdf) for a further description.
```
[wages1, frmDensity1, wrkDensity1, errWages, errDensityFrm, errSupplySkills] =
  geq(env,par,init,sol,cRobot,wages0,frmDensity0,wrkDensity0,tStart) 
```

#### Description of arguments

`env,par,init,sol`
  : input structures; see [codebook](https://andershumlum.com/s/humlumJMP.pdf)

`cRobot`
  : robot adoption costs wage path, <img src="https://render.githubusercontent.com/render/math?math=c^R_{t}">
   
`wages0`
  : current guess of wage path, <img src="https://render.githubusercontent.com/render/math?math=w^{(j)}_t">
  
`frmDensity0`
  : Initial firm states, <img src="https://render.githubusercontent.com/render/math?math=\mu^{F0}_t">

`wrkDensity0`
  : Initial worker states, <img src="https://render.githubusercontent.com/render/math?math=\mu^{W0}_t">
  
`tStart`
  : initial period, <img src="https://render.githubusercontent.com/render/math?math=t_0">

#### Description of output

  `wages1`
  : wage path that clears the static labor demand period-by-period, <img src="https://render.githubusercontent.com/render/math?math=w^{(e)}_t">
  
`frmDensity1`
  : simulated densities of firm states, <img src="https://render.githubusercontent.com/render/math?math=\mu^{F(j)}_t">

`wrkDensity1`
  : simulated densities of worker states, <img src="https://render.githubusercontent.com/render/math?math=\mu^{W(j)}_t">

`errWages`
  : Distance between current guess and static-equilibrium wage paths, <img src="https://render.githubusercontent.com/render/math?math=||w^{(e)}_{t}-w^{(j)}_{t}||">

`errDensityFrm`
  : Distance between current and new paths of firm densities, <img src="https://render.githubusercontent.com/render/math?math=||\mu^{F(e)}_{t} - \mu^{F(j)}_{t}||">
  
`errSupplySkills`
  : Distance between current and new paths for the labor supply to occupations, <img src="https://render.githubusercontent.com/render/math?math=||L^{S(e)}_{t} - L^{S(j)}_{t}||">


### `solve.frm`
Solve for firm continuation values using dynamic programming. See [Online Appendix F.1](https://andershumlum.com/s/humlumJMP.pdf) for further description.
```
v = frm(env,par,sol,profit,cRobot,tStart)
```
#### Description of arguments

`env,par,sol`
  : input structures; see [codebook](https://andershumlum.com/s/humlumJMP.pdf)

`profit`
  : flow profit function, <img src="https://render.githubusercontent.com/render/math?math=\pi_{t}(R,\varphi)">

`cRobot`
  : robot adoption cost, <img src="https://render.githubusercontent.com/render/math?math=c^R_{t}">
    
`tStart`
  : initial period, <img src="https://render.githubusercontent.com/render/math?math=t_0">
    

#### Description of output

  `v`
  : firm expected value function, <img src="https://render.githubusercontent.com/render/math?math=\mathbb{E} V_t (R,\varphi)">



### `solve.wrk`
Solve for worker continuation values using dynamic programming. See [Online Appendix F.2](https://andershumlum.com/s/humlumJMP.pdf) for further description.
```
[v, policy] = wrk(env,par,wages,tStart)
```
#### Description of arguments

`env,par`
  : input structures; see [codebook](https://andershumlum.com/s/humlumJMP.pdf)

`wages`
  : real wage path, <img src="https://render.githubusercontent.com/render/math?math=w_{t}">
  
`tStart`
  : initial period, <img src="https://render.githubusercontent.com/render/math?math=t_0">
    

#### Description of output

  `v`
  : worker expected value function, <img src="https://render.githubusercontent.com/render/math?math=\mathbb{E} V_t (o,a,\omega)">
  
`policy`
  : policy function for occupational choice, <img src="https://render.githubusercontent.com/render/math?math=o'(o,a,\omega)">

### `simulate.frm`
Simulate firm states forward.
```
[density, adopt0] = frm(env,par,densityInit,cRobot,v,tStart)
```
#### Description of arguments

`env,par`
  : input structures; see [codebook](https://andershumlum.com/s/humlumJMP.pdf)

`densityInit`
  : Initial conditions for firm states
  
`cRobot`
  : Robot adoption cost, <img src="https://render.githubusercontent.com/render/math?math=c^R_t">
    
`v`
  : Firm value function, <img src="https://render.githubusercontent.com/render/math?math=\mathbb{E} V_t (R,\varphi)">
      
`tStart`
  : initial period, <img src="https://render.githubusercontent.com/render/math?math=t_0">
    

#### Description of output

  `density`
  : simulated densities of firm states, <img src="https://render.githubusercontent.com/render/math?math=d\mu^F_t">
  
`adopt0`
  : robot adoption policy function, <img src="https://render.githubusercontent.com/render/math?math=R_t(0,\varphi)">


### `simulate.wrk`
Simulate worker states forward.
```
density = wrk(env,par,init,densityInit,policy,tStart)
```
#### Description of arguments

`env,par,init`
  : input structures; see [codebook](https://andershumlum.com/s/humlumJMP.pdf)

`densityInit`
  : Initial conditions for firm states
  
`policy`
  : policy function for occupational choice, <img src="https://render.githubusercontent.com/render/math?math=o'(o,a,\omega)">
      
`tStart`
  : initial period, <img src="https://render.githubusercontent.com/render/math?math=t_0">
    

#### Description of output

  `density`
  : simulated densities of worker states, <img src="https://render.githubusercontent.com/render/math?math=d\mu^W_t">
