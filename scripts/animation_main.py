import os
import csv
import time
import vtk  # Unused, but Required by TVTK.
from tvtk.api import tvtk
import numpy as np
from numpy.typing import NDArray
import quaternionic
import spherical
from mayavi import mlab
from mayavi.api import Engine
from mayavi.sources.vtk_data_source import VTKDataSource
from mayavi.sources.parametric_surface import ParametricSurface
from mayavi.modules.surface import Surface
import psi4_FFI_to_strain

BH_DIR = ""
ELL_MAX = 8
ELL_MIN = 2
S_MODE = -2


def swsh_summation_angles(colat: float, azi: NDArray[np.float64], mode_data):
    """
    Adds up all the strain modes after factoring in corresponding spin-weighted spherical harmonic
    to specified angle in the mesh. Stored as an array corresponding to [angle, time] time_idxs.

    :param colat: colatitude angle for the SWSH factor
    :param azi: azimuthal angle for the SWSH factor
    :param n_times: number of time time_idxs in the strain data
    :param mode_data: dictionary containing strain data for all the modes
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
    time_array: NDArray[np.float64], radius_values: NDArray[np.float64], r_ext: int
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


def initialize_tvtk_grid(num_azi, num_radius):
    """
    Sets initial parameters for the mesh generation module and returns
    mesh manipulation objects to write and save data.

    :param num_azi: number of azimuthal points on the mesh
    :param num_radius: number of radial points on the mesh
    :returns: tvtk.FloatArray(),
              tvtkUnstructuredGrid(),
              tvtkPoints()

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
    grid.set_cells(tvtk.Quad().cell_type, cell_array)

    return strain_array, grid, points


def create_gw(
    engine: Engine, grid: tvtk.UnstructuredGrid, color: tuple[float, float, float]
):
    scene = engine.scenes[0]
    gw = VTKDataSource(data=grid)
    engine.add_source(gw, scene)
    s = Surface()
    engine.add_filter(s, gw)
    s.actor.mapper.scalar_visibility = False
    s.actor.property.color = color
    return


def create_sphere(
    engine: Engine, radius: float = 1, color: tuple[float, float, float] = (1, 0, 0)
):
    """
    Creates and displays a parametric surface with the given parameters.
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


def change_object_position(obj: Surface, position: tuple[float, float, float]):
    """
    Changes the Cartesian position of a Mayavi object to the given parameters.
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


def dhms_time(seconds):
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


def main():
    """
    Main function that reads the strain data,
    calculates and factors in spin-weighted spherical harmonics,
    linearly interpolates the strain to fit the mesh points,
    and creates .vtu mesh file for each time time_idx of the simulation.
    The meshes represent the full superimposed waveform at the polar angle pi/2,
    aka the same plane as the binary black hole merger.
    """

    # File names
    bh_file_name = "bh_synthetic.csv"
    bh_file_path = os.path.join(BH_DIR, bh_file_name)
    movie_file_name = "test"
    movie_file_path = os.path.join(BH_DIR, "movies", movie_file_name)

    if not os.path.exists(movie_file_path):  # Create the directory if it doesn't exist
        os.makedirs(movie_file_path)

    # Mathematical parameters
    n_rad_pts = 450
    n_azi_pts = 180
    display_radius = 300
    ext_rad = 100
    amplitude_scale_factor = 600
    omitted_radius_length = 20
    colat = np.pi / 2  # colatitude angle representative of the plane of merger

    # Cosmetic parameters
    status_messages = True
    save_rate = 10  # Saves every Nth frame
    gw_color = (0.28, 0.46, 1.0)
    bh_color = (0.1, 0.1, 0.1)
    bh1_mass = 1.24
    bh2_mass = 1
    bh_scaling_factor = 1

    # Import strain data
    time_array, mode_data = psi4_FFI_to_strain.psi4_ffi_to_strain()
    n_times = len(time_array)

    # Import black hole data
    with open(bh_file_path, mode="r", encoding="utf-8") as file:
        reader = csv.reader(file)
        _ = next(reader)  # skip the header row
        bh_data = np.array([row for row in reader])

    # Pre-compute theta and radius values for the mesh
    radius_values = np.linspace(0, display_radius, n_rad_pts)
    azimuth_values = np.linspace(0, 2 * np.pi, n_azi_pts, endpoint=False)

    rv, az = np.meshgrid(radius_values, azimuth_values, indexing="ij")
    x_values = rv * np.cos(az)
    y_values = rv * np.sin(az)

    # Apply spin-weighted spherical harmonics, superimpose modes, and interpolate to mesh points
    strain_to_mesh = np.zeros((n_rad_pts, n_azi_pts, n_times))

    strain_azi = swsh_summation_angles(colat, azimuth_values, mode_data).real

    lerp_times = generate_interpolation_points(time_array, radius_values, ext_rad)
    for i in range(n_azi_pts):
        strain_to_mesh[:, i, :] = np.interp(lerp_times, time_array, strain_azi[i, :])

    strain_array, grid, points = initialize_tvtk_grid(n_azi_pts, n_rad_pts)

    # Create Mayavi objects
    engine = Engine()
    engine.start()
    engine.new_scene()
    engine.scenes[0].scene.jpeg_quality = 100
    # mlab.options.offscreen = True
    create_gw(engine, grid, gw_color)
    bh1 = create_sphere(engine, bh1_mass * bh_scaling_factor, bh_color)
    bh2 = create_sphere(engine, bh2_mass * bh_scaling_factor, bh_color)
    mlab.view(
        azimuth=60, elevation=50, distance=80, focalpoint=(0, 0, 0)
    )  # Initialize viewpoint

    start_time = time.time()
    percentage = np.round(np.linspace(0, n_times / save_rate, 101)).astype(int)

    @mlab.animate()  # ui=False) This doesn't work for some reason?
    def anim():
        for time_idx, t, row in zip(range(n_times), time_array, bh_data):
            if time_idx % save_rate != 0:
                continue

            # Print status messages
            if time_idx == 10:
                end_time = time.time()
                eta = (end_time - start_time) * n_times / 10
                print(
                    f"""Creating {n_times} meshes and saving them to: {movie_file_path}\nEstimated time: {dhms_time(eta)} minutes"""
                )

            if status_messages and time_idx != 0 and np.isin(time_idx, percentage):
                print(f" {int(time_idx * 100 / (n_times - 1))}% done", end="\r")

            # Change the position of the black holes
            _ = float(row[0])  # time
            bh1_xyz = (float(row[1]), float(row[2]), float(row[3]))
            bh2_xyz = (float(row[4]), float(row[5]), float(row[6]))
            change_object_position(bh1, bh1_xyz)
            change_object_position(bh2, bh2_xyz)

            points.reset()
            index = 0
            for j, radius in enumerate(radius_values):
                for i, azi in enumerate(azimuth_values):
                    x = x_values[j, i]
                    y = y_values[j, i]
                    # Introduce a discontinuity to make room for the Black Holes
                    if radius <= omitted_radius_length:
                        strain_value = np.nan
                    else:
                        strain_value = strain_to_mesh[i][time_idx][j].real
                    z = strain_value * amplitude_scale_factor
                    points.insert_next_point(x, y, z)
                    strain_array.set_tuple1(index, strain_value)
                    index += 1

            grid._set_points(points)
            grid._get_point_data().add_array(strain_array)
            grid.modified()

            mlab.view(
                azimuth=min(60 + time_idx * 0.018, 78),
                elevation=max(50 - time_idx * 0.016, 34),
                distance=min(80 + time_idx * 0.35, 430),
            )

            # Save the frame
            frame_path = os.path.join(
                movie_file_path, f"{movie_file_name}{time_idx:05d}.png"
            )
            mlab.savefig(frame_path)

            if time_idx == n_times - 1:
                total_time = time.time() - start_time
                print("100% Done")
                print(
                    f"\nSaved {n_times} frames to {movie_file_path} in {dhms_time(total_time)}."
                )
                exit(0)
            yield

    anim()
    mlab.show()


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    main()
