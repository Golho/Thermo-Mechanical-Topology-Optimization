include("GmshParser.jl")
include("HeatFEM.jl")

using WriteVTK
using LinearAlgebra

g = parse_mesh("../meshes/square_01x01_005.msh")

time_steps = 100
t_final = 0.001
const DUMP = false

prescribed = BoundaryCondition(
    "Dirichlet", "outlet", 0, 1:time_steps)

flux = BoundaryCondition("Neumann", "inlet",
    0.1/0.001, 1:time_steps)

body = BodyCondition("main", "solid")

material = Material(0.001, 1, 1, 1*Matrix(I, 2, 2))

fem = HeatFEM(g, t_final, time_steps)
fem.material = material
push!(fem.boundary_conditions, prescribed)
push!(fem.boundary_conditions, flux)
push!(fem.body_conditions, body)

# Apply conditions (assemble)

apply_body_conditions!(fem)
@time apply_boundary_conditions!(fem)
solve!(fem)

if DUMP
    pvd = paraview_collection("result/my_pvd_file")

    cell_type = VTKCellTypes.VTK_QUAD
    connectivity = get_elements(g, QUA_4)[:, 2:end]

    cells = [MeshCell(cell_type, connectivity[e, :]) for e = 1:size(connectivity, 1)]

    points = fem.node_coordinates'

    times = LinRange(0, t_final, time_steps)
    for time_step = 1:time_steps
        vtkfile = vtk_grid("result/my_vtk_file_$time_step", points, cells)
        vtkfile["temperature"] = fem.temperatures[:, time_step]
        pvd[times[time_step]] = vtkfile
    end

    outfiles = vtk_save(pvd)
end
