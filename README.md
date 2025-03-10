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

| Model       | Type               | Notes                     | Website                                  |
|-------------|--------------------|---------------------------|------------------------------------------|
| CoLM-Lake   | 1D Thermodynamic   | Integrated with CoLM      | [CoLM Official](https://github.com/CoLM-SYSU/CoLM202X) |
| FLake       | 1D Thermal         | Built-in configuration    | [FLake Wiki](http://www.flake.igb-berlin.de/) |
| Simstrat    | 1D/3D Hydrodynamic | Requires JSON input       | [Simstrat GitHub](https://github.com/Eawag-AppliedSystemAnalysis/Simstrat) |
