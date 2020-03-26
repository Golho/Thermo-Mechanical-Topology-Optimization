using DelimitedFiles

export ElementType, parse_mesh, Gmsh

@enum ElementType::Int begin
    LIN_2 = 1
    TRI_3 = 2
    QUA_4 = 3
    TET_4 = 4
    HEX_5 = 5
    PRI_6 = 6
    PYR_5 = 7
    LIN_3 = 8
    TRI_6 = 9
    QUA_9 = 10
    TET_10 = 11
    HEX_27 = 12
    PRI_18 = 13
    PYR_14 = 14
    PNT = 15
end

struct MeshFormat
    version::AbstractFloat
    filetype::Integer
    datasize::Integer
end

struct PhysicalName
    dimension::Integer
    physical_tag::Integer
    name::AbstractString
end

struct Point
    tag::Integer
    x::Float64
    y::Float64
    z::Float64
    physical_tags::AbstractArray{Int}
end

struct Curve
    tag::Integer
    min_x::Float64
    min_y::Float64
    min_z::Float64
    max_x::Float64
    max_y::Float64
    max_z::Float64
    physical_tags::AbstractArray{Int}
    point_tags::AbstractArray{Int}
end

struct Surface
    tag::Integer
    min_x::Float64
    min_y::Float64
    min_z::Float64
    max_x::Float64
    max_y::Float64
    max_z::Float64
    physical_tags::AbstractArray{Int}
    curve_tags::AbstractArray{Int}
end

struct Volume
    tag::Integer
    min_x::Float64
    min_y::Float64
    min_z::Float64
    max_x::Float64
    max_y::Float64
    max_z::Float64
    physical_tags::AbstractArray{Int}
    surface_tags::AbstractArray{Int}
end

struct Node
    node_tag::Integer
    x::Float64
    y::Float64
    z::Float64
end

struct NodeBlock
    entity_dim::Integer
    entity_tag::Integer
    parametric::Bool
    nodes::AbstractArray{Node, 1}
end

struct Element
    element_tag::Integer
    node_tags::AbstractArray{Int}
end

struct ElementBlock
    entity_dim::Integer
    entity_tag::Integer
    element_type::ElementType
    elements::AbstractArray{Element, 1}

    ElementBlock(a, b, c, d) = new(a, b, ElementType(c), d)
end

mutable struct Gmsh
    format::MeshFormat
    physical_names::AbstractArray{PhysicalName, 1}

    points::AbstractArray{Point, 1}
    curves::AbstractArray{Curve, 1}
    surfaces::AbstractArray{Surface, 1}
    volumes::AbstractArray{Volume, 1}

    node_blocks::AbstractArray{NodeBlock, 1}
    element_blocks::AbstractArray{ElementBlock, 1}

    Gmsh() = new()
end



function parse_mesh(filename::AbstractString)
    #@info filename
    file = open(filename, "r")
    gmsh = Gmsh()
    while !eof(file)
        blockname = readline(file)
        if blockname == "\$MeshFormat"
            gmsh.format = parse_meshformat(file)
        elseif blockname == "\$PhysicalNames"
            gmsh.physical_names = parse_physical_names(file)
        elseif blockname == "\$Entities"
            entities = parse_entities(file)
            gmsh.points = entities[1]
            gmsh.curves = entities[2]
            gmsh.surfaces = entities[3]
            gmsh.volumes = entities[4]
        elseif blockname == "\$Nodes"
            gmsh.node_blocks = parse_nodes(file)
        elseif blockname == "\$Elements"
            gmsh.element_blocks = parse_elements(file)
        end
    end
    close(file)

    return gmsh
end

function parse_meshformat(file::IOStream)
    meshformat_data = parse_vector(file)
    @info "GMSH-file version $(meshformat_data[1])"
    return MeshFormat(meshformat_data...)
end

function parse_physical_names(file::IOStream)
    num_physical_names = parse_vector(file; type=Int)[1]
    @info "Found $num_physical_names physical names:"
    physical_names = Array{PhysicalName, 1}(undef, num_physical_names)
    for i = 1:num_physical_names
        line = readline(file)
        line_vector = split(line, ' ')
        numericals = map(x -> parse(Float64, x), line_vector[1:2])
        name = trim(line_vector[3])
        @info "\t $name"
        physical_names[i] = PhysicalName(numericals..., name)
    end
    return physical_names
end

function parse_entities(file::IOStream)
    # <num-points num-curves num-surfaces num-volumes>
    num_data = parse_vector(file; type=Int)
    @info "Found $(num_data[1]) points"
    @info "Found $(num_data[2]) curves"
    @info "Found $(num_data[3]) surfaces"
    @info "Found $(num_data[4]) volumes"
    points = Array{Point, 1}(undef, num_data[1])
    curves = Array{Curve, 1}(undef, num_data[2])
    surfaces = Array{Surface, 1}(undef, num_data[3])
    volumes = Array{Volume, 1}(undef, num_data[4])

    for i = 1:num_data[1]
        data = parse_vector(file)
        num_physical_tags = convert(Int, data[5])
        points[i] = Point(data[1:4]..., convert.(Int, data[6:5 + num_physical_tags]))
    end

    function fill_entities!(entities, type::DataType, index::Integer, data)
        num_physical_tags = convert(Int, data[8])
        num_bounding_tags = convert(Int, data[8 + num_physical_tags + 1])

        if num_physical_tags > 0
            # Bounding tags exist
            physical_tags = convert.(Int, data[9:9+num_physical_tags-1])
        else
            physical_tags = Array{Int,1}(undef, 0)
        end
        if num_bounding_tags > 0
            # Bounding tags exist
            bounding_tags = convert.(Int, data[8+num_physical_tags+2:end])
        else
            bounding_tags = Array{Int,1}(undef, 0)
        end
        entities[index] = type(data[1:7]..., physical_tags, bounding_tags)
    end

    for i = 1:num_data[2]
        fill_entities!(curves, Curve, i, parse_vector(file))
    end
    for i = 1:num_data[3]
        fill_entities!(surfaces, Surface, i, parse_vector(file))
    end
    for i = 1:num_data[4]
        fill_entities!(volumes, Volume, i, parse_vector(file))
    end

    return (points, curves, surfaces, volumes)
end

function parse_nodes(file::IOStream)
    # <numEntityBlocks> <numNodes> <minNodeTag> <maxNodeTag>
    num_data = parse_vector(file; type=Int)
    @info "Detected $(num_data[2]) nodes"
    node_blocks = Array{NodeBlock, 1}(undef, num_data[1])

    for i = 1:num_data[1]
        # <entityDim> <entityTag> <parametric> <numNodesInBlock>
        entity_data = parse_vector(file; type=Int)
        if Bool(entity_data[3])
            @error "The parametric option is not yet supported by the parser"
        end
        node_tags = Array{Int, 1}(undef, entity_data[4])
        for j = 1:entity_data[4]
            node_tags[j] = parse_vector(file; type=Int)[1]
        end
        nodes = Array{Node, 1}(undef, entity_data[4])
        for j = 1:entity_data[4]
            pos = parse_vector(file)
            nodes[j] = Node(node_tags[j], pos...)
        end
        node_blocks[i] = NodeBlock(entity_data[1:3]..., nodes)
    end

    return node_blocks
end

function parse_elements(file::IOStream)
    # <numEntityBlocks> <numElements> <minElementTag> <maxElementTag>
    num_data = parse_vector(file; type=Int)
    @info "Detected $(num_data[2]) elements"
    element_blocks = Array{ElementBlock, 1}(undef, num_data[1])

    for i = 1:num_data[1]
        # <entityDim> <entityTag> <elementType> <numElementsInBlock>
        entity_data = parse_vector(file; type=Int)
        elements = Array{Element, 1}(undef, entity_data[4])

        for j = 1:entity_data[4]
            elem_data = parse_vector(file; type=Int)
            elements[j] = Element(elem_data[1], elem_data[2:end])
        end

        element_blocks[i] = ElementBlock(entity_data[1:3]..., elements)
    end

    return element_blocks
end
# Helper functions

function parse_vector(file::IOStream; delim::AbstractChar = ' ',
        type::DataType = Float64)
    line = readline(file)
    line_vector = split(line, delim)
    line_vector = filter(x -> x != "", line_vector)
    float_vector = map(x -> parse(type, x), line_vector)
end

function trim(str)
    str[2:end-1]
end

function get_physical_element_blocks(gmsh::Gmsh, physical_name::AbstractString)
    found_groups = filter(pn -> pn.name == physical_name, gmsh.physical_names)
    if isempty(found_groups)
        @error "No group with the physical name $physical_name found"
    else
        @debug "Found $physical_name"
    end

    blocks = Array{ElementBlock, 1}(undef, 0)
    found_element_blocks = Array{ElementBlock, 1}(undef, 0)
    for group = found_groups
        if group.dimension == 0
            elementary_entities = gmsh.points
        elseif group.dimension == 1
            elementary_entities = gmsh.curves
        elseif group.dimension == 2
            elementary_entities = gmsh.surfaces
        elseif group.dimension == 3
            elementary_entities = gmsh.volumes
        else
            @error "The dimension parameter of the physical group is invalid"
        end

        # Find which entities that are included in the physical group
        found_entities = filter(e -> any(e.physical_tags .== group.physical_tag), elementary_entities)

        # Extract the element blocks which are included in the found entities
        correct_entity(e) = any(e.entity_tag .== (p->p.tag).(found_entities)) &&
            e.entity_dim == group.dimension

        append!(found_element_blocks, filter(correct_entity, gmsh.element_blocks))
    end

    return found_element_blocks
end

function get_global_coordinates(gmsh::Gmsh)
    number_nodes = get_number_nodes(gmsh)
    node_coordinates = zeros(number_nodes, 3)
    for block = gmsh.node_blocks, node = block.nodes
        node_coordinates[node.node_tag, 1] = node.x
        node_coordinates[node.node_tag, 2] = node.y
        node_coordinates[node.node_tag, 3] = node.z
    end

    return node_coordinates
end

function get_number_nodes(gmsh::Gmsh)
    number_nodes = 0
    for block = gmsh.node_blocks
        number_nodes += length(block.nodes)
    end
    return number_nodes
end

function get_elements(gmsh::Gmsh, element_type::ElementType)
    blocks = filter(block -> block.element_type == element_type, gmsh.element_blocks)
    Enod = Array{Array{Float64}, 1}(undef, 0)

    function concat_elements(element_block)
        element_tags = [e.element_tag for e in element_block.elements]
        node_tags = transpose(hcat([e.node_tags for e in element_block.elements]...))
        return hcat(element_tags, node_tags)
    end

    Enod = vcat(map(concat_elements, blocks)...)
end
