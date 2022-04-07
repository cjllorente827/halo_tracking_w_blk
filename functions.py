from audioop import reverse
import numpy as np
import yt
from scipy.interpolate import interp1d
yt.set_log_level(40)

class HaloData:

    def __init__(self, id, Mvir, Rvir, X, Y, Z, catalog_fname):
        self.id = int(id)
        self.Mvir = Mvir
        self.Rvir = Rvir
        self.X = X
        self.Y = Y
        self.Z = Z
        self.catalog_fname = catalog_fname

    def __str__(self):
        return f"Halo {self.id}"

##############################################################################
# Returns all halos in the ascii catalog in a dictionary indexed by the 
# halo id.
##############################################################################
def load_catalog(hc_fname, box_length):
    # Units: Masses in Msun / h
    # Units: Radii in kpc / h (comoving)
    # Units: Positions in Mpc / h (comoving)
    #
    # Distances will be returned in code units, i.e. 
    # floating point numbers from 0 to 1, where 1
    # represents the full length of the box

    ids, DescID, Mvir, Vmax, Vrms, Rvir, Rs, Np, X, Y, Z, *other = np.genfromtxt(hc_fname,
        skip_header=1,
        unpack=True)

    result = {}
    for i, halo_id in enumerate(ids):
        result[halo_id] = HaloData(
            ids[i], 
            Mvir[i], 
            Rvir[i]/1000/box_length, 
            X[i]/box_length, 
            Y[i]/box_length, 
            Z[i]/box_length,
            hc_fname)

    return result

def load_halo(hc_fname, target_id, box_length):
    # Units: Masses in Msun / h
    # Units: Radii in kpc / h (comoving)
    # Units: Positions in Mpc / h (comoving)
    #
    # Distances will be returned in code units, i.e. 
    # floating point numbers from 0 to 1, where 1
    # represents the full length of the box

    ids, DescID, Mvir, Vmax, Vrms, Rvir, Rs, Np, X, Y, Z, *other = np.genfromtxt(hc_fname,
        skip_header=1,
        unpack=True)

    result = None
    for i, halo_id in enumerate(ids):
        if halo_id == target_id:
            return HaloData(
                ids[i], 
                Mvir[i], 
                Rvir[i]/1000/box_length, 
                X[i]/box_length, 
                Y[i]/box_length, 
                Z[i]/box_length,
                hc_fname)

    print(f"No halo found with id {target_id}")
    return None

def get_halo_particles(dataset, halo_position, halo_radius):
    
    ds = yt.load(dataset)
    sp = ds.sphere(
        (halo_position[0],halo_position[1],halo_position[2]), 
        halo_radius
    )

    particle_type = np.array(sp[('all', 'particle_type')], dtype=int)

    particle_data = {}
    particle_data["index"] = np.array(sp[('all', 'particle_index')], dtype=int)
    # particle_data["type"]    = np.array(sp[('all', 'particle_type')], dtype=int)
    particle_data["mass"]  = np.array(sp[('all', 'particle_mass')].to('Msun'))
    particle_data["x"]  = np.array(sp[('all', 'particle_position_x')])
    particle_data["y"]  = np.array(sp[('all', 'particle_position_y')])
    particle_data["z"]  = np.array(sp[('all', 'particle_position_z')])
    
    return particle_data

# The unique_to parameter is used to assign some unique parameter that 
# distinguishes this run from other runs, since its likely the user
# will want to pull different data from the same dataset.
# Most likely, you will want to use the halo id as a unique parameter
def get_particles_by_id(data, data_key, enzo_dataset):
    
    ds = yt.load(enzo_dataset)
    ad = ds.all_data()

    
    particles = data[data_key]
    particle_ids = particles["index"]

    all_particles = np.array(ad[('all', 'particle_index')], dtype=int)
    array_mask = [True if particle_id in particle_ids else False for particle_id in all_particles]

    particle_data = {}
    particle_data["index"] = np.array(ad[('all', 'particle_index')], dtype=int)[array_mask]
    # particle_data["type"]    = np.array(ad[('all', 'particle_type')], dtype=int)[array_mask]
    particle_data["mass"]  = np.array(ad[('all', 'particle_mass')].to('Msun'))[array_mask]
    particle_data["x"]  = np.array(ad[('all', 'particle_position_x')])[array_mask]
    particle_data["y"]  = np.array(ad[('all', 'particle_position_y')])[array_mask]
    particle_data["z"]  = np.array(ad[('all', 'particle_position_z')])[array_mask]

    return particle_data


def calc_particle_COM(data, data_key):

    particles = data[data_key]

    M = np.sum(particles["mass"])
    
    x_com = np.sum(particles["x"] * particles["mass"]) / M
    y_com = np.sum(particles["y"] * particles["mass"]) / M
    z_com = np.sum(particles["z"] * particles["mass"]) / M

    return (x_com, y_com, z_com)

# Attempts to determine halo properties from an initial guess location in a 
# simulation dataset
def find_center(
    data,
    data_key,
    dataset, 
    initial_radius,
    threshold, 
    max_iterations = 25,
    reduce_radius_by = 0.75,
    output_all_data=False):

    location = data[data_key]

    ds = yt.load(dataset)

    iterations = 0
    
    eps = 1
    centers = [location]


    while (iterations < max_iterations and eps > threshold):

        radius = initial_radius*(reduce_radius_by**iterations)

        
        sphere = ds.sphere(centers[iterations],(radius, "code_length"))
        
        com = sphere.quantities.center_of_mass(use_particles=True)
        centers.append(
            np.array(com.to('code_length').value)
        )

        eps = np.linalg.norm(centers[iterations+1] - centers[iterations])

        # For debug purposes 
        # plot = yt.SlicePlot(
        #     ds, "z", ("gas", "density"),
        #     center=sphere.center,
        #     width=(2*radius, "code_length"),
        #     data_source=sphere
        # )

        # plot.save(f"iteration_{iterations}")

        iterations += 1

    # if the user wants to see the result of every iteration, e.g. for debug/convergence testing
    # return the result of every iteration
    if output_all_data:
        return centers
    # otherwise, just return the final result
    else:
        return centers[-1]

def create_tracking_file(centers, 
    fname, 
    target_halo_id, 
    redshift_range,
    redshifts,
    box_length, 
    refine_level,
    use_time=False,
    dataset=None):
    
    key = lambda z: f"z_{z}_halo_{target_halo_id}_center"
    redshift_range.sort(reverse=True)
    redshifts.sort(reverse=True)
    
    positions = np.zeros((len(centers.keys()),3))
    for i,z in enumerate(redshifts):
        positions[i] = centers[key(z)][:]

    if use_time:
        ds = yt.load(dataset)
        t_from_z = lambda z : ds.cosmology.t_from_z(z)
        times = t_from_z(np.array(redshifts))
        time_range = t_from_z(np.array(redshift_range))

        interpolated_x = interp1d(times, positions[:,0], fill_value="extrapolate")
        interpolated_y = interp1d(times, positions[:,1], fill_value="extrapolate")
        interpolated_z = interp1d(times, positions[:,2], fill_value="extrapolate")

        extrapolated_positions = np.array([
            interpolated_x(time_range),
            interpolated_y(time_range),
            interpolated_z(time_range)
        ]).T

    else:
        interpolated_x = interp1d(redshifts, positions[:,0], fill_value="extrapolate")
        interpolated_y = interp1d(redshifts, positions[:,1], fill_value="extrapolate")
        interpolated_z = interp1d(redshifts, positions[:,2], fill_value="extrapolate")

        extrapolated_positions = np.array([
            interpolated_x(redshift_range),
            interpolated_y(redshift_range),
            interpolated_z(redshift_range)
        ]).T

    delta = np.ones(3) * box_length
    with open(fname, "w") as outfile:
        for i in range(extrapolated_positions.shape[0]):
            center = extrapolated_positions[i]
            z = redshift_range[i]
            lower_left  = center - delta
            upper_right = center + delta

            next_line = f"{z:.2f} {lower_left[0]:.7f} {lower_left[1]:.7f} {lower_left[2]:.7f} {upper_right[0]:.7f} {upper_right[1]:.7f} {upper_right[2]:.7f} {refine_level:.0f}\n"
            outfile.write(next_line)

        # end for
    # end with

