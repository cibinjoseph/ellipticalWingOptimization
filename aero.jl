using VortexLattice
using Plots
using Trapz

function getEllipticalChords(n; semispan = 4.0, root_chord = 1.0)
    y = LinRange(0, 1.0, n)
    return root_chord*sqrt.(1 .- y.^2), y*semispan
end

function getGrid(chords, nc = 4; semispan = 4.0)
    ns = length(chords)-1
    xyz = zeros(3, nc+1, ns+1)
    yvec = LinRange(0.0, semispan, ns+1)

    for i in 1:nc+1
        xyz[2, i, :] .= yvec
    end
    for i in 1:ns+1
        xyz[1, :, i] .= collect(LinRange(0.0, chords[i], nc+1))
    end

    # Make the quarter-chord line zero sweep
    for i in 1:ns+1
        xyz[1, :, i] .-= chords[i]*0.25
    end

    return xyz

end

function aerodynamics(inputs; semispan = 4.0, output_plot = false)
    alpha = 5.0*pi/180
    chords = inputs
    ns = length(chords)

    # # geometry (right half of the wing)
    xle = -0.25 .* chords
    yle = LinRange(0, semispan, ns)
    zle = zeros(ns)
    chord = chords
    theta = 0.0*pi/180 * ones(ns)
    phi = zeros(ns)
    #
    # # discretization parameters
    nc = 4
    spacing_s = Uniform()
    spacing_c = Uniform()

    # reference parameters
    cref = sum(chords)/ns
    bref = semispan * 2
    Sref = trapz(yle, chords)*2
    rref = zeros(3)
    rho = 1.2
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    # freestream parameters
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    # construct surface
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
                                           spacing_s=spacing_s, spacing_c=spacing_c)

    # xyz = getGrid(chords)
    # grid, surface = grid_to_surface_panels(xyz)

    # create vector containing all surfaces
    surfaces = [surface]

    # we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
    symmetric = true

    # perform steady state analysis
    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)

    # retrieve near-field forces
    CF, CM = body_forces(system; frame=Wind())

    # perform far-field analysis
    # CDiff = far_field_drag(system)

    # Plot wing
    if output_plot
        properties = get_surface_properties(system)
        write_vtk("wing", surfaces, properties; symmetric)
    end

    CD, CY, CL = CF
    L = CL * 0.5*rho*Vinf^2*Sref
    D = CD * 0.5*rho*Vinf^2*Sref

    return L, D, CL, CD, Sref
end

function elliptical(output_plot=false)
    chords, y = getEllipticalChords(80)
    # plot(y, chords)

    L, D, CL, CD, S = aerodynamics(chords; output_plot=output_plot)
    return L, D, CL, CD, S
end
