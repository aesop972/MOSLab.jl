#include("MOSLab1.jl")
   using MosLab
   using LaTeXStrings
   using Gadfly
   using Cairo
   #Vg = -1.0 # Gate Voltage
   #Vd = -1.00 # Dain Voltage
   Nx = 51 # Number of grid points on the x direction
   Ny = 21 # Number of grid points on the y direction
   N_d = 2e19 # Dopants density at the drain and source side  in cm ^-3
   T = 300.0 # Temperature in Kelvin
   Nb = 1e18 # Dopants ensity at the bulk in cm ^-3
   l = .03e-4 # Channel length
   h = 0.02e-4 # Silicon depth
   t_ox = 1.5e-7 # Oxide thickens in cm
   e_charge = 1.602e-19 # Electron charge
   Nikel() = Metal(5.1*Eunit)
   Ms = MOSFETInputDeck(N_d,Nb,T,l,h,h/2,l,t_ox,:PMOS;lscal=1.0,ϕₘ=5.15) # Create the input deck for an NMOS transistor
   NB = SemiconductorData(T,BoltzmanDist(),NSilicon(Nb,0))
   MOS1 = MOSStructure(Nikel(),SiO2(),NB,t_ox,h)
   
   dRang = [0, 1, 2, 3, 4, 5]
   nCarrier = zeros(length(dRang))
   pCarrier = zeros(length(dRang))
   eleC = zeros(length(dRang),length(dRang))
   holC = zeros(length(dRang),length(dRang))
   Vs0 = ψs_PSPp(0,MOS1) # PMOS
   G0,psi0,nn0,pp0 = MOSFETSimulation(Ms,Vs0,0,Nx,Ny;verbose=false,ξ₀=0.01)
   nCarrier0 = sum(nn0)
   pCarrier0 = sum(nn0)
   for Vg in 0:5
        sVg = string(Vg)
        Vga = Vg * (-0.20)
        Vs = ψs_PSPp(Vga,MOS1) # PMOS
        for Vd in 0:5

            sVd = string(Vd)
            Vda = Vd * (-0.20)
            G,psi,nn,pp = MOSFETSimulation(Ms,Vs,Vda,Nx,Ny;verbose=false,ξ₀=0.01)
            nCarrier[Vd+1] = sum(nn) - nCarrier0
            pCarrier[Vd+1] = sum(pp) - pCarrier0

            mpsi = reshape(reinterpret(Float64, psi),51,21)       
            coord1 = Coord.cartesian(xmin=0, xmax=0.09, ymin=0, ymax=0.02)
            myplot = plot(coord1, z=mpsi, x=0.0:0.0018:0.09, y=0.0:0.001:0.02, 
            Geom.vectorfield(scale=0.00004), Geom.contour,
            Scale.x_continuous(minvalue=0.0, maxvalue=0.09),
            Scale.y_continuous(minvalue=0.0, maxvalue=0.02), 
            Guide.xlabel("x"), Guide.ylabel("y"),
            layer(x=[0.03,0.06,0.06,0.03],y=[0.02,0.02,0.021,0.021],Geom.polygon(preserve_order=true, fill=true)))
            draw(PNG(join([sVg,sVd,".png"]), 5inch, 3.5inch), myplot)

        end 

        aRang = dRang * -0.2
        eleC[Vg+1,:] = nCarrier .* aRang
        eleC[Vg+1,:] = eleC[Vg+1,:] * e_charge * l * h *1e10 # unit uA/um
        holC[Vg+1,:] = pCarrier .* aRang
        holC[Vg+1,:] = holC[Vg+1,:] * e_charge * l * h * 1e10 # unit uA/um
    end
   plot(
        layer(x=aRang, y=holC[1,:]-eleC[1,:]),
        layer(x=aRang, y=holC[2,:]-eleC[2,:]),
        layer(x=aRang, y=holC[3,:]-eleC[3,:]),
        layer(x=aRang, y=holC[4,:]-eleC[4,:]),
        layer(x=aRang, y=holC[5,:]-eleC[5,:]),
        )