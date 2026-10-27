# Taylor–Green vortex reference data (Re = 1600)

- `TGV_Re1600.dat`: DNS with Incompact3d on 512³ (only the π³ box is computed, using symmetries), from
  T. Dairay, E. Lamballais, S. Laizet & J. C. Vassilicos, "Numerical dissipation vs. subgrid-scale modelling for large
  eddy simulation", *J. Comput. Phys.* 337, 252–274 (2017), https://doi.org/10.1016/j.jcp.2017.02.035.
  Cite this paper when using the data; see the file header for the columns.
- `p7.jld2`: WaterLily.jl run at 128³, with time `t`, kinetic energy `E` and dissipation
  `Z = ν⟨ω²⟩`, from the WaterLily.jl CPC paper (G. D. Weymouth & B. Font).
