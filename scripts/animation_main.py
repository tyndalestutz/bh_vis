"""
A script that reads in gravitational wave strain data and black hole positional data,
applies spin-weighted spherical harmonics to the data, and creates a Mayavi animation
of the black holes and their gravitational waves. At each state, the black holes are 
moved to their respective positions and the render is saved as a .png file.
"""

import os
import sys
import csv
import cv2
import time
import quaternionic
import spherical
import vtk  # Even though we don't use it directly, TVTK requires it
import numpy as np
from scipy.interpolate import interp1d
from math import erf
from typing import Tuple, Dict, Any
from numpy.typing import NDArray
from tvtk.api import tvtk
from mayavi import mlab
from mayavi.api import Engine
from mayavi.sources.vtk_data_source import VTKDataSource
from mayavi.sources.parametric_surface import ParametricSurface
from mayavi.modules.surface import Surface

import util_psi4_to_strain_incl_checks as psi4_to_strain


def read_strain_files(
    file_path: str,
) -> Tuple[NDArray[np.float64], Dict[Tuple[int, int], NDArray[np.complex128]]]:
    """
    Read an ASCII file with a header describing the real and imaginary parts of the data
    for a specific l mode. Returns the time states and mode data in an easily retrievable format.

    :param file_path: Path to the file to be read.
    :return: A tuple containing the time (numpy array)
             and a dictionary with keys (l, m) containing the data.
    :raises ValueError: if the length of the time data is inconsistent across different ell values.
    """
    mode_data: Dict[Tuple[int, int], NDArray[np.complex128]] = {}
    time_data_size: int = -1
    for ell in range(2, 9):
        file_name = file_path.replace("[ELLVAL]", str(ell))
        with open(file_name, "r", encoding="utf-8") as file:
            # Read lines that don't start with '#'
            print(f"Reading strain data from {file_name}")
            lines = [line for line in file.readlines() if not line.startswith("#")]

        # Convert lines to arrays and sort by time
        data: NDArray[np.float64] = np.array(
            [list(map(np.float64, line.split())) for line in lines]
        )
        data = data[np.argsort(data[:, 0])]

        # Remove duplicate times
        _, index = np.unique(data[:, 0], return_index=True)
        data = data[index]

        # Store time data and check for discrepancies between modes
        time_data: NDArray[np.float64] = data[:, 0]
        if time_data_size < 0:
            time_data_size = len(time_data)
        elif time_data_size != len(time_data):
            raise ValueError(
                f"""Inconsistent time data size for ell={ell}. 
                Expected {time_data_size}, got {len(time_data)}."""
            )

        # Loop through columns and store real and imaginary parts in dictionary
        for em in range(-ell, ell + 1):
            mode_index: int = 1 + 2 * (em + ell)  # the index for the real valued strain
            mode_data[(ell, em)] = data[:, mode_index] + 1j * data[:, mode_index + 1]

    return time_data, mode_data


def find_swsh_factor(colat: float, azi: float, ell: int, em: int) -> Any:
    """
    Calculates the complex valued spin-weighted spherical harmonic (SWSH) factors to be
    multiplied with strain. Uses spherical package to find the factor and the inputted
    colatitude, azimuthal, l, and m values, and a physically defined spin = -2 value.

    :param colat: colatitude angle for the SWSH factor
    :param azi: azimuthal angle for the SWSH factor
    :param ell: wave mode value l
    :param em: wave mode value m
    :return: a complex valued spin-weighted spherical harmonic factor

    >>> find_swsh_factor(np.pi/2, 0, 5, 0)
    (-0+0j)
    >>> find_swsh_factor(np.pi/2, 5*np.pi/7, 2, 2)
    (-0.035090612830967316-0.15374202011569824j)
    """
    s = -2
    R = quaternionic.array.from_spherical_coordinates(colat, azi)
    winger = spherical.Wigner(8)  # Create a Winger matrix for all modes from l=2 to l=8
    Y = winger.sYlm(s, R)
    swsh_factor = Y[winger.Yindex(ell, em)]

    return swsh_factor


def superimpose_modes_from_angle(
    colat: float,
    azi: float,
    num_time_states: int,
    mode_data: Dict[Tuple[int, int], NDArray[np.complex128]],
):
    """
    Adds up all the strain modes after factoring in corresponding spin-weighted spherical harmonic
    to specified angle in the mesh. Stored as a new array corresponding to time states.

    :param colat: colatitude angle for the SWSH factor
    :param azi: azimuthal angle for the SWSH factor
    :param num_time_states: number of time states in the strain data
    :param mode_data: dictionary containing strain data for all the modes
    :return: a complex valued numpy array of the superimposed wave
    
    >>> mode_data = {}
    >>> for l in range(2, 9):
    ...     for m in range(-l, l+1):
    ...         mode_data[(l, m)] = np.array([1+1j, 2+3j, 4+5j])
    >>> superimpose_modes_from_angle(np.pi/2, 0, 3, mode_data)
    array([ 4.69306041 +4.69306041j,  9.38612083+14.07918124j,
           18.77224166+23.46530207j])
    """
    summation: NDArray[np.int64] = np.zeros(num_time_states, dtype="complex128")
    for ell in range(2, 9):
        for em in range(-ell, ell + 1):
            swsh_factor = find_swsh_factor(colat, azi, ell, em)
            factored_strain = mode_data[(ell, em)] * swsh_factor
            summation += factored_strain
    return summation


def generate_gw_interpolation_points(
    time_array: NDArray[np.float64], radius_values: NDArray[np.float64], r_ext: int
) -> NDArray[np.float64]:
    """
    Fills out a 2D array of adjusted time values for the wave strain to be
    linearly interpolated to. First index of the result represents the simulation
    time state (aka which mesh), and the second index represents radial distance to
    interpolate to.

    :param time_array: numpy array of of strain time states.
    :param radius_values: numpy array of the radial points on the mesh.
    :param r_ext: extraction radius of the original data.
    :return: a 2D numpy array of time values.
    """

    time_0 = np.min(time_array)
    time_f = np.max(time_array)
    target_times = np.zeros((len(time_array), len(radius_values)))
    for state, current_time in enumerate(time_array):
        for j, radius in enumerate(radius_values):
            target_time = current_time - radius + r_ext
            if target_time < time_0:
                target_time = time_0
            elif target_time > time_f:
                target_time = time_f
            target_times[state][j] = target_time
    return target_times


def interpolate_coords_by_time(
    og_time_array: NDArray[np.float64], 
    e1: NDArray[np.float64], 
    e2: NDArray[np.float64],
    e3: NDArray[np.float64],
    new_time_array: float) -> Tuple:
    """
    Interpolates the 3D coordinates to the given time state.
    :param og_time_array: original time array
    :param e1: first coordinate array
    :param e2: second coordinate array
    :param e3: third coordinate array
    :param new_time_array: new time array
    :return: interpolated 3D coordinates
    """

    new_e1 = interp1d(og_time_array, e1, kind="linear", fill_value="extrapolate")(new_time_array)
    new_e2 = interp1d(og_time_array, e2, kind="linear", fill_value="extrapolate")(new_time_array)
    new_e3 = interp1d(og_time_array, e3, kind="linear", fill_value="extrapolate")(new_time_array)
    return new_e1, new_e2, new_e3


def initialize_tvtk_grid(num_azi: int, num_radius: int) -> Tuple:
    """
    Sets initial parameters for the mesh generation module and returns
    a circular, polar mesh with manipulation objects to write and save data.

    :param num_azi: number of azimuthal points on the mesh
    :param num_radius: number of radial points on the mesh
    :returns: tvtk.FloatArray,
             tvtk.UnstructuredGrid,
             tvtk.Points
    """
    # Create tvtk objects
    points = tvtk.Points()
    grid = tvtk.UnstructuredGrid()
    strain_array = tvtk.FloatArray(
        name="Strain", number_of_components=1, number_of_tuples=num_azi * num_radius
    )

    # Create cells
    cell_array = tvtk.CellArray()
    for j in range(num_radius - 1):
        for i in range(num_azi):
            cell = tvtk.Quad()
            point_ids = [
                i + j * num_azi,
                (i + 1) % num_azi + j * num_azi,
                (i + 1) % num_azi + (j + 1) * num_azi,
                i + (j + 1) * num_azi,
            ]
            for idx, pid in enumerate(point_ids):
                cell.point_ids.set_id(idx, pid)
            cell_array.insert_next_cell(cell)

    # Set grid properties
    # grid.points = points
    grid.set_cells(tvtk.Quad().cell_type, cell_array)

    return strain_array, grid, points


def create_gw(
    engine: Engine,
    grid: Any,
    color: Tuple[float, float, float],
    wireframe: bool = False
) -> None:
    """
    Creates and displays a gravitational wave strain from a given grid.
    :param engine: Mayavi engine
    :param grid: tvtk.UnstructuredGrid
    :param color: color of the strain in a tuple ranging from (0, 0, 0) to (1, 1, 1)
    :param wireframe: whether to display the strain as a wireframe or a surface
    """

    scene = engine.scenes[0]
    gw = VTKDataSource(data=grid)
    engine.add_source(gw, scene)
    s = Surface()
    engine.add_filter(s, gw)
    s.actor.mapper.scalar_visibility = False
    s.actor.property.color = color
    if wireframe:
        s.actor.property.representation = "wireframe"
        s.actor.property.color = (0, 0, 0)
        s.actor.property.line_width = 0.005
        s.actor.property.opacity = 0.5


def create_sphere(
    engine: Engine, radius: float = 1, color: tuple[float, float, float] = (1, 0, 0)
) -> Surface:
    """
    Creates and displays a spherical surface with the given parameters.
    :param engine: Mayavi engine
    :param radius: radius of the sphere
    :param color: color of the sphere in a tuple ranging from (0, 0, 0) to (1, 1, 1)
    """
    scene = engine.scenes[0]
    ps = ParametricSurface()
    ps.function = "ellipsoid"
    ps.parametric_function.x_radius = radius
    ps.parametric_function.y_radius = radius
    ps.parametric_function.z_radius = radius

    engine.add_source(ps, scene)
    s = Surface()
    engine.add_filter(s, ps)
    s.actor.mapper.scalar_visibility = False
    s.actor.property.color = color
    return s


def change_object_position(obj: Surface, position: tuple[float, float, float]) -> None:
    """
    Changes the Cartesian position of a Mayavi surface to the given parameters.
    :param engine: Mayavi engine
    :param obj: Mayavi object
    :param position: position of the object
    """
    position = np.array(position)
    obj.actor.actor.position = position


def dhms_time(seconds: float) -> str:
    """
    Converts a given number of seconds into a string indicating the remaining time.
    :param seconds: number of seconds
    :return: a string indicating the remaining time
    """

    days_in_seconds = 24 * 60 * 60
    hours_in_seconds = 60 * 60
    minutes_in_seconds = 60

    # Calculate remaining days, hours, and minutes
    days = int(seconds // days_in_seconds)
    remaining_seconds = seconds % days_in_seconds
    hours = int(remaining_seconds // hours_in_seconds)
    remaining_seconds = remaining_seconds % hours_in_seconds
    minutes = int(remaining_seconds // minutes_in_seconds)

    # Build the output string
    output = ""
    if days > 0:
        output += f"{days} days"
    if hours > 0:
        if days > 0:
            output += " "
        output += f"{hours} hours"
    if minutes > 0:
        if days > 0 or hours > 0:
            output += " "
        output += f"{minutes} minutes"
    return output


def convert_to_movie(input_path: str, movie_name: str, fps: int = 24) -> None:
    """
    Converts a series of .png files into a movie using OpenCV.
    :param input_path: path to the directory containing the .png files
    :param movie_name: name of the movie file
    :param fps: frames per second (24 by default)
    """
    frames = [f for f in os.listdir(input_path) if f.endswith(".png")]
    frames.sort()
    # Create a movie from the frames
    ref = cv2.imread(os.path.join(input_path, frames[0]))
    height, width, layers = ref.shape
    video = cv2.VideoWriter(
        os.path.join(input_path, f"{movie_name}.mp4"),
        cv2.VideoWriter_fourcc(*"mp4v"),
        fps,
        (width, height)
    )
    for frame in frames:
        f = cv2.imread(os.path.join(input_path, frame))
        video.write(f)
    video.release()


def main() -> None:
    """
    Main function that reads the strain data,
    calculates and factors in spin-weighted spherical harmonics,
    linearly interpolates the strain to fit the mesh points,
    and creates .tvtk mesh file for each time state of the simulation.
    The meshes represent the full superimposed waveform at the polar angle pi/2,
    aka the same plane as the binary black hole merger. At each state, the black holes
    are moved to their respective positions and the mesh is saved as a .png file.
    """

    if len(sys.argv) != 3:
        file_name = os.path.basename(__file__)
        raise RuntimeError(
            f"""Please include path to psi4 folder data as well as the extraction radius of that data.
            Usage: python3 {file_name} <path to psi4 folder> <extraction radius (r/M) (4 digits, e.g. 0100)>"""
        )
    psi4_folder_path = str(sys.argv[1])
    radius_of_extraction = int(sys.argv[2])
    #print(f"EXTRACTION RADIUS: {float(radius_of_extraction)}")
    # Convert psi4 data to strain using imported script
    #psi4_to_strain.main()
    
    # File names
    gw_file_name = "Rpsi4_r0100.0_l[ELLVAL]_conv_to_strain.txt"
    gw_file_path = os.path.abspath(
        os.path.join(psi4_folder_path, gw_file_name)
    )
    bh_file_name = "../puncture_posns_vels_regridxyzU.txt"
    bh_file_path = os.path.abspath(
        os.path.join(__file__, "..", bh_file_name)
    )
    movie_file_name = "first_real_movie"
    movie_file_path = os.path.abspath(
        os.path.join(__file__, "..", "..", "..", "movies", movie_file_name)
    )
    if os.path.exists(movie_file_path):  # Check if the directory exists 
        r = input(f"""{movie_file_path} already exists. Would you like to overwrite it? Y or N: """)
        if r.lower() != "y":
            print("Please choose a different directory.")
            exit()
        for file in os.listdir(movie_file_path):
            os.remove(os.path.join(movie_file_path, file))
    else:
        os.makedirs(movie_file_path)
    
    # Mathematical parameters
    num_radius_points = 450
    num_azi_points = 180
    display_radius = 300
    radius_of_extraction = 100
    amplitude_scale_factor = 200
    omitted_radius_length = 10

    colat = np.pi / 2  # colatitude angle representative of the plane of merger

    # Cosmetic parameters
    status_messages = True
    wireframe = True
    frames_per_second = 24
    save_rate = 10  # Saves every Nth frame
    resolution = (1920, 1080)
    gw_color = (0.28, 0.46, 1.0)
    bh_color = (0.1, 0.1, 0.1)
    bh1_mass = 1
    bh2_mass = 1.24 
    bh_scaling_factor = 1

    # Import strain data
    time_array, mode_data = read_strain_files(gw_file_path)
    num_time_states = len(time_array)
    effective_num_time_states = int(num_time_states / save_rate)

    # Pre-compute theta and radius values for the mesh
    radius_values = np.linspace(0, display_radius, num_radius_points)
    azimuth_values = np.linspace(0, 2 * np.pi, num_azi_points, endpoint=False)
    rv, az = np.meshgrid(radius_values, azimuth_values, indexing="ij")
    x_values = rv * np.cos(az)
    y_values = rv * np.sin(az)

    print("""**********************************************************************
Constructing mesh points in 3D...""")
    strain_to_mesh = (
        {}
    )  # Holds the final strain points indexed [azimuthal point (key)][time state][radius]

    # Apply spin-weighted spherical harmonics, superimpose modes, and interpolate to mesh points
    interpolation_times = generate_gw_interpolation_points(
        time_array, radius_values, radius_of_extraction
    )
    for azi_index, azi in enumerate(azimuth_values):
        superimposed_strain = superimpose_modes_from_angle(
            colat, azi, num_time_states, mode_data
        )
        strain_to_mesh[azi_index] = np.interp(
            interpolation_times, time_array, superimposed_strain
        ).real  

    # Import black hole data
    if bh1_mass > bh2_mass:
        (bh1_mass, bh2_mass) = (bh2_mass, bh1_mass) # Swap the masses so that the heavier black hole is bh2
    with open(bh_file_path, "r") as file:
        reader = csv.reader(file, delimiter=" ")
        #_ = next(reader)  # uncomment to skip the header row
        bh_data = np.array(list(reader)).astype(np.float64)
        bh_time = bh_data[:, 0]
        
        bh1_x = -bh_data[:, 5] # x is flipped because the data is in a different coordinate system
        bh1_y = bh_data[:, 7] # z axis in the data is interpreted as y axis in the visualization
        bh1_z = np.zeros(len(bh1_x))
    bh1_x, bh1_y, bh1_z = interpolate_coords_by_time(bh_time, bh1_x, bh1_y, bh1_z, time_array)

    bh_mass_ratio = bh1_mass / bh2_mass
    bh2_x = -bh1_x * bh_mass_ratio
    bh2_y = -bh1_y * bh_mass_ratio
    bh2_z = -bh1_z * bh_mass_ratio

    # Create Mayavi objects
    #mlab.options.offscreen = True
    engine = Engine()
    engine.start()
    #engine.new_scene()
    #engine.scenes[0].scene.jpeg_quality = 100
    mlab.figure(engine=engine, size=resolution)
    strain_array, grid, points = initialize_tvtk_grid(num_azi_points, num_radius_points)
    create_gw(engine, grid, gw_color)
    if wireframe:
        create_gw(engine, grid, gw_color, wireframe=True)
    bh1 = create_sphere(engine, bh1_mass * bh_scaling_factor, bh_color)
    bh2 = create_sphere(engine, bh2_mass * bh_scaling_factor, bh_color)
    mlab.view(
        azimuth=60, elevation=50, distance=80, focalpoint=(0, 0, 0)
    )  # Initialize viewpoint

    start_time = time.time()
    percentage = list(np.round(np.linspace(0, num_time_states, 100)).astype(int))

    @mlab.animate() # ui=False) This doesn't work for some reason?
    def anim():
        for state, t in enumerate(time_array):
            if state % save_rate != 0:
                continue # Skip all but every nth iteration

            # Print status messages
            if state == 10 * save_rate:
                end_time = time.time()
                eta = (end_time - start_time) * effective_num_time_states / 10
                print(
                    f"""Creating {effective_num_time_states} frames and saving them to: 
                    {movie_file_path}\nEstimated time: {dhms_time(eta)}"""
                )
            if status_messages and state != 0 and state > percentage[0]:
                eta = ((time.time() - start_time) / state) * (num_time_states - state)
                print(
                    f"{int(state  * 100 / num_time_states)}% done, ", 
                    f"{dhms_time(eta)} remaining",
                    end="\r",
                )
                percentage.pop(0)

            # Change the position of the black holes
            bh1_xyz = (bh1_x[state], bh1_y[state], bh1_z[state])
            bh2_xyz = (bh2_x[state], bh2_y[state], bh2_z[state])
            change_object_position(bh1, bh1_xyz)
            change_object_position(bh2, bh2_xyz)

            points.reset()
            index = 0
            for j, radius in enumerate(radius_values):
                for i, _ in enumerate(azimuth_values):
                    x = x_values[j, i]
                    y = y_values[j, i]
                    if radius <= omitted_radius_length:
                        z = np.nan
                        strain_value = np.nan
                    else:
                        h_t_real = strain_to_mesh[i][state][j]
                        dropoff_radius = 1.5 * omitted_radius_length
                        width = 0.5 * omitted_radius_length
                        dropoff_factor = 0.5*(erf((radius - dropoff_radius) / width) + 1)
                        strain_value = h_t_real * amplitude_scale_factor * dropoff_factor
                        z = strain_value
                    points.insert_next_point(x, y, z)
                    strain_array.set_tuple1(index, strain_value)
                    index += 1
            grid._set_points(points)
            grid._get_point_data().add_array(strain_array)
            grid.modified()

            mlab.view(
                #azimuth=min(60 + state * 0.018, 78),
                elevation = max(50 - state * 0.016, 34),
                distance = 80 if state < 2000 else min(80 + (state - 2000) * 0.175, 370),
                focalpoint =(0, 0, 0),
            )

            # Save the frame
            frame_num = int(state/save_rate)
            frame_path = os.path.join(
                movie_file_path, f"z_frame_{frame_num:05d}.png"
            )
            mlab.savefig(frame_path, magnification=1)

            if state >= (effective_num_time_states * save_rate) - 1: # Smoothly exit the program
                total_time = time.time() - start_time
                print(f"Done", end="\r")
                print(
                    f"\nSaved {effective_num_time_states} frames to {movie_file_path} ",
                    f"in {dhms_time(total_time)}."
                )
                print(f"Creating movie...")
                convert_to_movie(movie_file_path, movie_file_name, frames_per_second)
                print(f"Movie saved to {movie_file_path}/{movie_file_name}.mp4")
                exit()
            yield

    _ = anim()
    mlab.show()

# This should automatically create the movie file, but if it doesn't work, run this in the movie directory:
# $ffmpeg -framerate 24 -i frame_%05d.png <movie_name>.mp4

if __name__ == "__main__":
    # run doctests first
    import doctest

    results = doctest.testmod()
    psi4_to_strain_results = doctest.testmod(psi4_to_strain)

    if psi4_to_strain_results.failed > 0:
        print(
            f"""Doctest in {psi4_to_strain} failed:
{psi4_to_strain_results.failed} of {psi4_to_strain_results.attempted} test(s) passed"""
        )
        sys.exit(1)
    else:
        print(f"""Doctest in {psi4_to_strain} passed:
All {psi4_to_strain_results.attempted} test(s) passed""")

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s) passed")
        #sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    # run main() after tests
    main()
