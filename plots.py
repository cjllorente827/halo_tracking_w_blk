import matplotlib.pyplot as plt
import numpy as np
from blk import cache


colors = [
        "red",
        "orange",
        "yellow",
        "green",
        "blue",
    ]

def particle_plot(stage):
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(projection='3d')

    particles = cache.load(stage.tag)
    show_fraction = 1 # 0.05

    total_particles = len(particles["x"])
    num_particle = int(np.floor(show_fraction * total_particles))
    print(f"Displaying {num_particle} particles out of {total_particles}")
    
    x = 25*particles["x"][:num_particle]
    y = 25*particles["y"][:num_particle]
    z = 25*particles["z"][:num_particle]

    ax.scatter(x,y,z,c="red")
    ax.set_xlabel('x (Mpc/h)')
    ax.set_ylabel('y (Mpc/h)')
    ax.set_zlabel('z (Mpc/h)')

    plt.ioff()

    plt.show()


def multiz_particle_plot(stage_list):
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(projection='3d')

    show_fraction = 1# 0.05

    for i, stage in enumerate(stage_list):
        particles = cache.load(stage.tag)
        
        num_particle = int(np.floor(show_fraction * len(particles["x"])))
        print(f"Displaying {num_particle} particles")
        
        x = 25*particles["x"][:num_particle]
        y = 25*particles["y"][:num_particle]
        z = 25*particles["z"][:num_particle]

        ax.scatter(x,y,z, c=colors[i])
    ax.set_xlabel('x (Mpc/h)')
    ax.set_ylabel('y (Mpc/h)')
    ax.set_zlabel('z (Mpc/h)')

    plt.ioff()
    plt.show()
    
def position_plot(stage_list):
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(projection='3d')

    for i,stage in enumerate(stage_list):
        position = 25 * np.array( cache.load(stage.tag) )
        
    
        ax.scatter(position[0], position[1], position[2], c=colors[i])
    ax.set_xlabel('x (Mpc/h)')
    ax.set_ylabel('y (Mpc/h)')
    ax.set_zlabel('z (Mpc/h)')

    plt.show()
    
def centers_plot(COM_stage_list, center_stage_list):
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(projection='3d')

    for i,com_stage in enumerate(COM_stage_list):
        center_stage = center_stage_list[i]
        position = 25 * np.array( cache.load(com_stage.tag) )
    
        ax.scatter(position[0], position[1], position[2], c=colors[i], marker='x')
        
        position = 25 * np.array( cache.load(center_stage.tag) )
    
        ax.scatter(position[0], position[1], position[2], c=colors[i])
    ax.set_xlabel('x (Mpc/h)')
    ax.set_ylabel('y (Mpc/h)')
    ax.set_zlabel('z (Mpc/h)')

    plt.show()

def plot_tracking_file(fname):

    rshft, Lx, Ly, Lz, Rx, Ry,Rz,ref = np.genfromtxt(fname, unpack=True)

    Cx, Cy, Cz = (Rx+Lx)/2,(Ry+Ly)/2,(Rz+Lz)/2

    more_colors = list(reversed([
        "black",
        "purple",
        "red",
        "orange",
        "yellow",
        "green",
        "blue",
    ]))

    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(projection='3d')

    for i, z in enumerate(rshft):
        
        position = 25 * np.array( [Cx[i], Cy[i],Cz[i]] )
        ax.scatter(position[0], position[1], position[2], c=more_colors[i])

    ax.set_xlabel('x (Mpc/h)')
    ax.set_ylabel('y (Mpc/h)')
    ax.set_zlabel('z (Mpc/h)')

    plt.show()