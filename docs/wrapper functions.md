# Mean, Linear, and Equipercentile Equating
## No Smoothing
- Random Groups Design: __Wrapper_RN__
- Single Groups Design: __Wrapper_SN__
- Common-Item Nonequivalent Groups Design: __Wrapper_CN__
## PreSmoothing
- Beta-Binomial
  - Random Groups Design: __Wrapper_RB__
- Log-Linear
  - Random Groups Design: __Wrapper_RL__
  - Single Groups Design: __Wrapper_SL__
  - Common-Item Nonequivalent Groups Design: __Wrapper_CL__

# Kernel Equating
- Random Groups Design: __Wrapper_RK__
- Single Groups Design: __Wrapper_SK__
- Common-Item Nonequivalent Groups Design: __Wrapper_CK__

# Continuized Log-Linear Equating
- Random Groups Design: __Wrapper_RC__
- Single Groups Design: __Wrapper_SC__
- Common-Item Nonequivalent Groups Design: __Wrapper_CC__

# Miscellanous
- Beta Binomial Smoothing: __Wrapper_Smooth_BB__
- Univariate Log Linear Smoothing: __Wrapper_Smooth_ULL__
- Bivariate Log Linear Smoothing: __Wrapper_Smooth_BLL__
- Cubic Spline Postsmoothing: __Wrapper_Smooth_CubSpl__
- Equated Scale Scores: __Wrapper_ESS__
- Raw and Scale Score Bootstrap Standard Errors: __Wrapper_Bootstrap__
- IRT scale transformation: __Wrapper_IRTst__
- IRT equating: __Wrapper_IRTeq__

# Print Functions
- Wrapper_RN(): __Print_RN__
- Wrapper_SN(): __Print_SN__
- Wrapper_CN(): __Print_CN__, __Print_SynDens__
- Wrapper_RB(): __Print_RB__
- Wrapper_RL(): __Print_RL__
- Wrapper_SL(): __Print_SL__ 
- Wrapper_CL(): __Print_CL__, __Print_SynDens__
- Wrapper_ESS(): __Print_ESS__
- Wrapper_Bootstrap(): __Print_Boot_se_eraw__, __Print_Boot_se_ess__
- Wrapper_Smooth BB(): __Print_BB__
- Wrapper_Smooth ULL(): __Print_ULL__
- Wrapper_Smooth BLL(): __Print_BLL__
- Wrapper_Smooth CubSpl(): __Print_CubSpl__
- Wrapper_RK(): __Print_RK__
- Wrapper_SK(): __Print_SK__
- Wrapper_CK(): __Print_CK__
- Wrapper_RC(): __Print_RC__
- Wrapper_SC(): __Print_SC__
- Wrapper_CC(): __Print_CC__
- Wrapper_IRTst(): [wrapper does own printing]
- Wrapper_IRTeq(): __Print_IRTeq__, __Print_ESS_QD__

# Structures
- USTATS: raw-score statistics for a univariate distribution
- BSTATS: raw-score statistics for a bivariate distribution
- PDATA: input for a particular D/M/S schema as well as other data passed among functions
- ERAW_RESULTS: equated raw score results
- ESS_RESULTS: equated scale score results
- BB_SMOOTH: input and output for beta-binomial smoothing
- ULL_SMOOTH: input and output for univariate log-linear smoothing
- BLL_SMOOTH: input and output for bivariate log-linear smoothing
- CS_SMOOTH: input and output for cubic spline postsmoothing
- BOOT_ERAW_RESULTS: equated raw score results for bootstrap
- BOOT_ESS_RESULTS: equated scale score results for bootstrap