# Lake Model Driver (LMD)  

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)  

A unified framework for integrating and running multiple lake models (e.g., FLake, Simstrat) with seamless coupling to the **Common Land Model (CoLM)** or standalone execution.  

---

## Features  
- **Unified Interface**: Standardized input/output formats and configurations for easy model switching.  
- **CoLM Coupling**: Bidirectional data exchange via adapter patterns (supports offline/online modes).  
- **Extensible Design**: Add new lake models with minimal effort.  
- **Toolchain**: Preprocessing, validation, and visualization utilities included.  

## Supported Models

| Model       | Vertical structure               | Turbulent mixing parameterization                                                                 | Website                                  |
|-------------|-----------------------------------|---------------------------------------------------------------------------------------------------|------------------------------------------|
| CoLM-Lake   | Multilayer                        | The water surface temperature is equal to the mixed-layer temperature<br>computed from heat flux | https://github.com/CoLM-SYSU/CoLM202X    |
| FLake       | Two-layer self-similar structure | Hendenson-Sellers thermal diffusion model<br>with wind-driven diffusivity                        | http://www.flake.igb-berlin.de           |
| Simstrat    | Multilayer                        | k-ε turbulence model with<br>buoyancy and internal seiche parameterization                      | https://github.com/Eawag-AppliedSystemAnalysis/Simstrat |
