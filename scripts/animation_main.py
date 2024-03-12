"""
A script that reads in gravitational wave strain data and black hole positional data,
applies spin-weighted spherical harmonics to the data, and creates a Mayavi animation
of the black holes and their gravitational waves. At each state, the black holes are 
moved to their respective positions and the render is saved as a .png file.
"""

import os
import sys
import csv
import time
from math import erf
from typing import Tuple, Any
import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import interp1d
import quaternionic
import spherical
import cv2
import vtk  # Unused, but Required by TVTK.
from tvtk.api import tvtk
from mayavi import mlab
from mayavi.api import Engine
from mayavi.sources.vtk_data_source import VTKDataSource
from mayavi.sources.parametric_surface import ParametricSurface
from mayavi.modules.surface import Surface
import psi4_FFI_to_strain as psi4strain


BH_DIR = "../BH_VIS/r100"
MOVIE_DIR = "movies"
ELL_MAX = 8
ELL_MIN = 2
S_MODE = -2
EXT_RAD = 100


def swsh_summation_angles(colat: float, azi: NDArray[np.float64], mode_data):
    """
    Adds up all the strain modes after factoring in corresponding spin-weighted spherical harmonic
    to specified angle in the mesh. Stored as an array corresponding to [angle, time] time_idxs.

    :param colat: colatitude angle for the SWSH factor
    :param azi: azimuthal angle for the SWSH factor
    :param mode_data: numpy array containing strain data for all the modes
    :return: a complex valued numpy array of the superimposed wave
    """
    quat_arr = quaternionic.array.from_spherical_coordinates(colat, azi)
    winger = spherical.Wigner(ELL_MAX, ELL_MIN)
    # Create an swsh array shaped like (n_modes, n_quaternions)
    swsh_arr = winger.sYlm(S_MODE, quat_arr).T
    # mode_data has shape (n_modes, n_times), swsh_arr has shape (n_mode, n_pts).
    # Pairwise multiply and sum over modes: the result has shape (n_pts, n_times).
    pairwise_product = mode_data[:, np.newaxis, :] * swsh_arr[:, :, np.newaxis]
    return np.sum(pairwise_product, axis=0)


def generate_interpolation_points(  # Could use some revisiting, currently keeps n_times constant
    time_array: NDArray[np.float64],
    radius_values: NDArray[np.float64],
    r_ext: int = EXT_RAD,
) -> NDArray[np.float64]:
    """
    Fills out a 2D array of adjusted time values for the wave strain to be
    linearly interpolated to. First index of the result represents the simulation
    time time_idx (aka which mesh), and the second index represents radial distance to
    interpolate to.

    :param time_array: numpy array of of strain time time_idxs.
    :param radius_values: numpy array of the radial points on the mesh.
    :param r_ext: extraction radius of the original data.
    :return: a 2D numpy array (n_radius, n_times) of time values.
    """
    # Repeat time_array and radius_values to match the shape of the output array
    time_repeated = np.repeat(time_array[np.newaxis, :], len(radius_values), axis=0)
    radius_repeated = np.repeat(radius_values[:, np.newaxis], len(time_array), axis=1)

    target_times = time_repeated - radius_repeated + r_ext
    # Shape is (n_radius, n_times)
    filtered_target_times = np.clip(target_times, time_array.min(), time_array.max())
    return filtered_target_times


def interpolate_coords_by_time(
    old_times: NDArray[np.float64],
    e1: NDArray[np.float64],
    e2: NDArray[np.float64],
    e3: NDArray[np.float64],
    new_times: float,
) -> Tuple:
    """
    Interpolates the 3D coordinates to the given time state.
    :param old_times: original time array
    :param e1: first coordinate array
    :param e2: second coordinate array
    :param e3: third coordinate array
    :param new_times: new time array
    :return: interpolated 3D coordinates
    """

    new_e1 = interp1d(old_times, e1, fill_value="extrapolate")(new_times)
    new_e2 = interp1d(old_times, e2, fill_value="extrapolate")(new_times)
    new_e3 = interp1d(old_times, e3, fill_value="extrapolate")(new_times)
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
    >>> strain_array, grid, points = initialize_tvtk_grid(3, 4)
    >>> isinstance(strain_array, tvtk.FloatArray)
    True
    >>> isinstance(grid, tvtk.UnstructuredGrid)
    True
    >>> isinstance(points, tvtk.Points)
    True
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
    wireframe: bool = False,
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


def change_view(
    engine: Engine,
    position: tuple[float, float, float] = None,
    focal_point: tuple[float, float, float] = None,
    view_up: tuple[float, float, float] = None,
    view_angle: float = None,
    clipping_range: tuple[float, float] = None,
):
    """
    Changes the view of the Mayavi engine to the given parameters.
    :param engine: Mayavi engine
    :param position: position of the camera, default is current position
    :param focal_point: focal point of the camera, default is current focal point
    :param view_up: view up vector of the camera, default is current view up vector
    :param view_angle: view angle of the camera, default is current view angle
    """

    scene = engine.scenes[0]
    if position is not None:
        scene.scene.camera.position = position
    if focal_point is not None:
        scene.scene.camera.focal_point = focal_point
    if view_up is not None:
        scene.scene.camera.view_up = view_up
    if view_angle is not None:
        scene.scene.camera.view_angle = 30.0
    if clipping_range is not None:
        scene.scene.camera.clipping_range = clipping_range
    scene.scene.camera.compute_view_plane_normal()


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
    height, width, _ = ref.shape
    video = cv2.VideoWriter(
        os.path.join(input_path, f"{movie_name}.mp4"),
        cv2.VideoWriter_fourcc(*"mp4v"),
        fps,
        (width, height),
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

    # Convert psi4 data to strain using imported script
    # psi4_to_strain.main()

    # File names
    bh_file_name = "puncture_posns_vels_regridxyzU.txt"
    bh_file_path = os.path.join(BH_DIR, bh_file_name)

    movie_file_name = "first_real_movie"
    movie_file_path = os.path.join(MOVIE_DIR, movie_file_name)

    if os.path.exists(movie_file_path):
        r = input(
            f"""{movie_file_path} already exists. Would you like to overwrite it? Y/N: """
        )
        if r.lower() != "y":
            print("Please choose a different directory.")
            exit()
        for file in os.listdir(movie_file_path):
            os.remove(os.path.join(movie_file_path, file))
    else:
        os.makedirs(movie_file_path)

    # Mathematical parameters
    n_rad_pts = 450
    n_azi_pts = 180
    display_radius = 300
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

    # ---Preliminary Calculations---

    strain_array, grid, points = initialize_tvtk_grid(n_azi_pts, n_rad_pts)
    width = 0.5 * omitted_radius_length
    dropoff_radius = width + omitted_radius_length

    # Import strain data
    time_array, mode_data = psi4strain.psi4_ffi_to_strain()
    n_times = len(time_array)
    n_frames = int(n_times / save_rate)

    # theta and radius values for the mesh
    radius_values = np.linspace(0, display_radius, n_rad_pts)
    azimuth_values = np.linspace(0, 2 * np.pi, n_azi_pts, endpoint=False)

    rv, az = np.meshgrid(radius_values, azimuth_values, indexing="ij")
    x_values = rv * np.cos(az)
    y_values = rv * np.sin(az)

    # Apply spin-weighted spherical harmonics, superimpose modes, and interpolate to mesh points
    strain_azi = swsh_summation_angles(colat, azimuth_values, mode_data).real
    lerp_times = generate_interpolation_points(time_array, radius_values)

    strain_to_mesh = np.zeros((n_rad_pts, n_azi_pts, n_times))
    for i in range(n_azi_pts):
        # strain_azi, a function of time_array, is evaluated at t = lerp_times.
        strain_to_mesh[:, i, :] = np.interp(lerp_times, time_array, strain_azi[i, :])


    # Import black hole data
    if bh1_mass > bh2_mass:  # then swap
        bh1_mass, bh2_mass = bh2_mass, bh1_mass

    print(
        """**********************************************************************
Constructing mesh points in 3D..."""
    )

    with open(bh_file_path, mode="r", encoding="utf-8") as file:
        reader = csv.reader(file, delimiter=" ")
        # _ = next(reader)  # uncomment to skip the header row
        bh_data = np.array(list(reader)).astype(np.float64)
    bh_time = bh_data[:, 0]
    # x is flipped because the data is in a different coordinate system
    bh1_x = -bh_data[:, 5]
    # z axis in the data is interpreted as y axis in the visualization
    bh1_y = bh_data[:, 7]
    bh1_z = np.zeros(len(bh1_x))
    bh1_x, bh1_y, bh1_z = interpolate_coords_by_time(
        bh_time, bh1_x, bh1_y, bh1_z, time_array
    )

    bh_mass_ratio = bh1_mass / bh2_mass
    bh2_x = -bh1_x * bh_mass_ratio
    bh2_y = -bh1_y * bh_mass_ratio
    bh2_z = -bh1_z * bh_mass_ratio

    # Create Mayavi objects
    # mlab.options.offscreen = True
    engine = Engine()
    engine.start()
    # engine.new_scene()
    # engine.scenes[0].scene.jpeg_quality = 100
    mlab.figure(engine=engine, size=resolution)
    strain_array, grid, points = initialize_tvtk_grid(n_azi_pts, n_rad_pts)

    create_gw(engine, grid, gw_color)
    if wireframe:
        create_gw(engine, grid, gw_color, wireframe=True)
    bh1 = create_sphere(engine, bh1_mass * bh_scaling_factor, bh_color)
    bh2 = create_sphere(engine, bh2_mass * bh_scaling_factor, bh_color)
    mlab.view(
        azimuth=60, elevation=50, distance=80, focalpoint=(0, 0, 0)
    )  # Initialize viewpoint

    start_time = time.time()
    percentage = list(np.round(np.linspace(0, n_times, 100)).astype(int))

    @mlab.animate()  # ui=False) This doesn't work for some reason?
    def anim():
        for time_idx in range(n_times):
            if time_idx % save_rate != 0:
                continue  # Skip all but every nth iteration

            # Print status messages
            if time_idx == 10 * save_rate:
                end_time = time.time()
                eta = (end_time - start_time) * n_frames / 10
                print(
                    f"""Creating {n_frames} frames and saving them to:
{movie_file_path}\nEstimated time: {dhms_time(eta)}"""
                )
            if status_messages and time_idx != 0 and time_idx > percentage[0]:
                eta = ((time.time() - start_time) / time_idx) * (n_times - time_idx)
                print(
                    f"{int(time_idx  * 100 / n_times)}% done, ",
                    f"{dhms_time(eta)} remaining",
                    end="\r",
                )
                percentage.pop(0)

            # Change the position of the black holes
            bh1_xyz = (bh1_x[time_idx], bh1_y[time_idx], bh1_z[time_idx])
            bh2_xyz = (bh2_x[time_idx], bh2_y[time_idx], bh2_z[time_idx])
            change_object_position(bh1, bh1_xyz)
            change_object_position(bh2, bh2_xyz)

            points.reset()
            index = 0
            for rad_idx, radius in enumerate(radius_values):
                dropoff_factor = 0.5 + 0.5 * erf((radius - dropoff_radius) / width)
                for azi_idx, _ in enumerate(azimuth_values):
                    x = x_values[rad_idx, azi_idx]
                    y = y_values[rad_idx, azi_idx]
                    if radius <= omitted_radius_length:
                        strain_value = np.nan
                    else:
                        strain_value = strain_to_mesh[rad_idx][azi_idx][time_idx]
                    z = strain_value * amplitude_scale_factor * dropoff_factor
                    points.insert_next_point(x, y, z)
                    strain_array.set_tuple1(index, strain_value)
                    index += 1

            grid._set_points(points)
            grid._get_point_data().add_array(strain_array)
            grid.modified()
            mlab.view(
                # azimuth=min(60 + time_idx * 0.018, 78),
                elevation=max(50 - time_idx * 0.016, 34),
                distance=(
                    80 if time_idx < 2000 else min(80 + (time_idx - 2000) * 0.175, 370)
                ),
                focalpoint=(0, 0, 0),
            )

            # Save the frame
            frame_num = int(time_idx / save_rate)
            frame_path = os.path.join(movie_file_path, f"z_frame_{frame_num:05d}.png")
            mlab.savefig(frame_path, magnification=1)

            if time_idx >= (n_frames * save_rate) - 1:  # Smoothly exit the program
                total_time = time.time() - start_time
                print("Done", end="\r")
                print(
                    f"\nSaved {n_frames} frames to {movie_file_path} ",
                    f"in {dhms_time(total_time)}.",
                )
                print("Creating movie...")
                convert_to_movie(movie_file_path, movie_file_name, frames_per_second)
                print(f"Movie saved to {movie_file_path}/{movie_file_name}.mp4")
                exit()
            yield

    _ = anim()
    mlab.show()


# This should automatically create the movie file...
# if it doesn't work, run the following in the movie directory:
# $ffmpeg -framerate 24 -i frame_%05d.png <movie_name>.mp4

if __name__ == "__main__":
    # run doctests first
    import doctest

    results = doctest.testmod()
    p4s_results = doctest.testmod(psi4strain)

    if p4s_results.failed > 0:
        print(
            f"""Doctest in {psi4strain} failed:
{p4s_results.failed} of {p4s_results.attempted} test(s) passed"""
        )
        sys.exit(1)
    else:
        print(
            f"""Doctest in {psi4strain} passed:
All {p4s_results.attempted} test(s) passed"""
        )

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
    # run main() after tests
    main()
