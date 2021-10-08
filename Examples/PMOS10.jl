   using MOSLab
   using Plots
   using LaTeXStrings
   Vg = -1.0 # Gate Voltage
   Vd = -1.0 # Dain Voltage
   Nx = 101 # Number of grid points on the x direction
   Ny = 101 # Number of grid points on the y direction
   N_d = 1e20 # Dopants density at the drain and source side  in cm ^-3
   T = 300.0 # Temperature in Kelvin
   Nb = 3e17 # Dopants ensity at the bulk in cm ^-3
   NNB = -3e17
   l = .065e-4 # Channel length
   h = 0.2e-4 # Silicon depth
   t_ox = 1.5e-7 # Oxide thickens in cm
   Ms = MOSFETInputDeck(N_d,Nb,T,l,h,h/2.0,l,t_ox,:PMOS;lscal=1.0,ϕₘ=5.15) # Create the input deck for an NMOS transistor
   #Ms = MOSFETInputDeck(N_d,Nb,T,l,h,h/2.0,l,t_ox,:PMOS;lscal=1.0,Alluminium())
   #Vs = ψs_PSP(Vg,Ms.Gate) # Calculate the surface potential at the gate
   NB = SemiconductorData(T,BoltzmanDist(),NSilicon(Nb,-0.044))
   Nikel() = Metal(5.1*Eunit)
   MOS1 = MOSStructure(Nikel(),SiO2(),NB,t_ox,h)
   Vs = ψs_PSPp(Vg,MOS1) 
   #G,psi,nn,pp = MOSFETSimulation(Ms,Vs,Vd,Nx,Ny;verbose=false,ξ₀=0.01) # run 2D MOSFET simulation
   G,psi,nn,pp = MOSFETSimulation(Ms,Vs,Vd,Nx,Ny;verbose=false,ξ₀=0.01)
   gr()
   XX = unique(G[1,:])
   YY = unique(G[2,:])
   contour(XX,YY,psi,fill=true)
   surface(G[1,:],G[2,:],psi,fill=true)
   xlabel!(L"x \ [\mu m]")
   ylabel!(L"y \ [\mu m]")
   plot!(zlabel = L"\psi \ [V]")
   savefig("PMOS10.png")