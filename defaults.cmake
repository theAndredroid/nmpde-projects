# Default Project Configuration

set(DEFAULT_BASE_MESH_FILE mesh/mesh-square)
set(DEFAULT_MESH_SIZE 0.1)
set(DEFAULT_MODEL "TNNP")
set(DEFAULT_DELTA_T 0.005)
set(DEFAULT_MAX_TIME 60000)
set(DEFAULT_THETA 0.5)

# Derived variables
set(MESH_COMPILING_COMMAND gmsh ${DEFAULT_BASE_MESH_FILE}.geo -save)
get_filename_component(MESH_STEM ${DEFAULT_BASE_MESH_FILE} NAME_WE)
set(DEFAULT_OUTPUT_FILE "output_ActivationTime-${MESH_STEM}")
