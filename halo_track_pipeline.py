
import numpy as np
from settings import *
import functions, plots


from blk import pipeline, cache
cache.set_dir(SCRATCHCACHE)


################################################################################
# First Iteration: Initial Particle Collection
################################################################################

# Grab the halo data for our target from the halo catalog
target_halo_id = 74
#target_halo_id = 6906
#target_halo_id = 1279
target_halo = functions.load_halo(HALO_CATALOGS[2], target_halo_id, 25)

# target_halo_id = 4954
# target_halo = functions.load_halo(ANNAS_CATALOG, target_halo_id, 25)

# set up the arguments for the particle collection stage
dataset = ENZO_OUTPUTS[2]
halo_position = (target_halo.X, target_halo.Y, target_halo.Z)
halo_radius = target_halo.Rvir

iteration_1 = pipeline.Stage(
    functions.get_halo_particles,
    [dataset, halo_position, halo_radius],
    tag=f"z_2_particles_halo_{target_halo_id}"
)


################################################################################
# Second Step: Find Particles at Earlier Times
################################################################################

iteration_2 = [iteration_1]

for i,z in enumerate(REDSHIFTS[1:]):
    
    dataset = ENZO_OUTPUTS[z]
    
    next_stage = pipeline.Stage(
        functions.get_particles_by_id,
        [iteration_1.tag, dataset],
        tag=f"z_{z}_particles_halo_{target_halo_id}",
        depends_on=[iteration_1]
    )

    iteration_2.append( next_stage )

################################################################################
# Third Step: Calculate COM for particles¶
################################################################################
iteration_3 = []

for i,z in enumerate(REDSHIFTS):
        
    next_stage = pipeline.Stage(
        functions.calc_particle_COM,
        [iteration_2[i].tag],
        tag=f"z_{z}_particle_COM_halo_{target_halo_id}",
        depends_on=[iteration_2[i]]
    )

    iteration_3.append( next_stage )

################################################################################
# Fourth Step: Run the Center Finder algorithm
################################################################################
iteration_4 = []

for i,z in enumerate(REDSHIFTS):
    
    dataset = ENZO_OUTPUTS[z]
    initial_radius = 1.5/25 # 1.5 Mpccm/h
    threshold = 1 / 256 / (2**7) # size of smallest refined cell in the simulation
    
    args = [
        iteration_3[i].tag,
        dataset,
        initial_radius,
        threshold
    ]
    
    next_stage = pipeline.Stage(
        functions.find_center,
        args,
        tag=f"z_{z}_halo_{target_halo_id}_center",
        depends_on=[iteration_3[i]]
    )

    iteration_4.append( next_stage )

################################################################################
# Final Step: Interpolate and Extrapolate to create the tracker file
################################################################################

filename = f"halo_{target_halo_id}_refine_region"
box_length = 1.5/25, # 1.5 Mpccm/h
refine_level = 4
redshift_range = [0,1,2,3,4,5,6]

args = [
    filename,
    target_halo_id,
    redshift_range,
    REDSHIFTS,
    box_length,
    refine_level,
    True,
    ENZO_OUTPUTS[2]
]

create_track_file = pipeline.Stage(
    functions.create_tracking_file,
    args,
    tag=filename,
    action="manual",
    depends_on=iteration_4
)

################################################################################
# Run the pipeline
################################################################################
pipeline.execute(create_track_file, parallelism="task")