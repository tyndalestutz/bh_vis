import os
import csv
import time
import vtk # Even though we don't use it directly, TVTK requires it
from tvtk.api import tvtk
from typing import Tuple, Dict, Any
import numpy as np
from numpy.typing import NDArray
import quaternionic
import spherical
from mayavi import mlab
from mayavi.api import Engine
from mayavi.sources.vtk_data_source import VTKDataSource
from mayavi.sources.parametric_surface import ParametricSurface
from mayavi.modules.surface import Surface


def read_strain_files(file_path: str) -> Tuple[NDArray[np.float64], Dict[Tuple[int, int], NDArray[np.complex128]]]:
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
            print(file_name)
            lines = [line for line in file.readlines() if not line.startswith("#")]

        # Convert lines to arrays and sort by time
        data: NDArray[np.float64] = np.array([list(map(np.float64, line.split())) for line in lines])
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
    """
    s = -2
    R = quaternionic.array.from_spherical_coordinates(colat, azi)
    winger = spherical.Wigner(8)  # Create a Winger matrix for all modes from l=2 to l=8
    Y = winger.sYlm(s, R)
    swsh_factor = Y[winger.Yindex(ell, em)]
    return swsh_factor


def superimpose_modes_from_angle(colat: float, azi: float, num_time_states: int, mode_data: Dict[Tuple[int, int], NDArray[np.complex128]]):
    """
    Adds up all the strain modes after factoring in corresponding spin-weighted spherical harmonic
    to specified angle in the mesh. Stored as a new array corresponding to time states.

    :param colat: colatitude angle for the SWSH factor
    :param azi: azimuthal angle for the SWSH factor
    :param num_time_states: number of time states in the strain data
    :param mode_data: dictionary containing strain data for all the modes
    :return: a complex valued numpy array of the superimposed wave
    """
    summation: NDArray[np.int64] = np.zeros(num_time_states, dtype="complex128")
    for ell in range(2, 9):
        for em in range(-ell, ell + 1):
            swsh_factor = find_swsh_factor(colat, azi, ell, em)
            factored_strain = mode_data[(ell, em)] * swsh_factor
            summation += factored_strain
    return summation


def generate_interpolation_points(time_array: NDArray[np.float64], radius_values: NDArray[np.float64], r_ext: int
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


def initialize_tvtk_grid(num_azi, num_radius):
   """
   Sets initial parameters for the mesh generation module and returns
   mesh manipulation objects to write and save data.

   :param num_azi: number of azimuthal points on the mesh
   :param num_radius: number of radial points on the mesh
   :returns: tvtk.api.DataArray(),
            tvtkUnstructuredGrid(),
            tvtkPoints()
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
         point_ids = [i + j * num_azi, (i + 1) % num_azi + j * num_azi,
            (i + 1) % num_azi + (j + 1) * num_azi, i + (j + 1) * num_azi]
         for idx, pid in enumerate(point_ids): 
               cell.point_ids.set_id(idx, pid)
         cell_array.insert_next_cell(cell)

   # Set grid properties
   #grid.points = points
   grid.set_cells(tvtk.Quad().cell_type, cell_array) 

   return strain_array, grid, points


def create_gw(engine, grid, color):
    scene = engine.scenes[0]
    gw = VTKDataSource(data = grid)
    engine.add_source(gw, scene)
    s = Surface()
    engine.add_filter(s, gw)
    s.actor.mapper.scalar_visibility = False
    s.actor.property.color = color
    return


def create_sphere(engine: Engine, radius: float = 1, color: tuple[float, float, float] = (1, 0, 0)):
    """
    Creates and displays a parametric surface with the given parameters.
    :param engine: Mayavi engine
    :param radius: radius of the sphere 
    :param color: color of the sphere in a tuple ranging from (0, 0, 0) to (1, 1, 1)
    """
    scene = engine.scenes[0]
    ps = ParametricSurface()
    ps.function = 'ellipsoid'
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


def change_view(engine: Engine, position: tuple[float, float, float] = None, focal_point: tuple[float, float, float] = None, view_up: tuple[float, float, float] = None, view_angle: float = None, clipping_range: tuple[float, float] = None):
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
    and creates .vtu mesh file for each time state of the simulation. 
    The meshes represent the full superimposed waveform at the polar angle pi/2,
    aka the same plane as the binary black hole merger.
    """

    # File names
    gw_file_name = "Rpsi4_r0100.0_l[ELLVAL]_conv_to_strain.txt"
    gw_file_path = os.path.abspath(
        os.path.join(__file__, "..", "..", "r100", gw_file_name)
    )
    bh_file_name = "bh_synthetic.csv"
    bh_file_path = os.path.abspath(
        os.path.join(__file__, "..", bh_file_name)
    )
    movie_file_name = "test"
    movie_file_path = os.path.abspath(
        os.path.join(__file__, "..", "..", "..", "movies", movie_file_name)
    )
    if not os.path.exists(movie_file_path): # Create the directory if it doesn't exist
        os.makedirs(movie_file_path)
    

    # Mathematical parameters
    num_radius_points = 450
    num_azi_points = 180
    display_radius = 300 
    R_extraction = 100
    amplitude_scale_factor = 600 
    omitted_radius_length = 20 
    colat = np.pi / 2  # colatitude angle representative of the plane of merger

    # Cosmetic parameters
    status_messages = True
    save_rate = 10 # Saves every Nth frame
    gw_color = (0.28, 0.46, 1.0)
    bh_color = (0.1, 0.1, 0.1)
    bh1_mass = 1.24 
    bh2_mass = 1
    bh_scaling_factor = 1 

    # Import strain data
    time_array, mode_data = read_strain_files(gw_file_path)
    num_time_states = len(time_array)

    # Import black hole data
    with open(bh_file_path, "r") as file:
        reader = csv.reader(file)
        _ = next(reader)  # skip the header row
        bh_data = np.array([row for row in reader])

    # Pre-compute theta and radius values for the mesh
    radius_values = np.linspace(0, display_radius, num_radius_points)
    azimuth_values = np.linspace(0, 2 * np.pi, num_azi_points, endpoint=False)
    rv, az = np.meshgrid(radius_values, azimuth_values, indexing="ij")
    x_values = rv * np.cos(az)
    y_values = rv * np.sin(az)
    
    strain_to_mesh = (
        {}
    )  # Holds the final strain points indexed [azimuthal point (key)][time state][radius]

    # Apply spin-weighted spherical harmonics, superimpose modes, and interpolate to mesh points
    interpolation_times = generate_interpolation_points(
        time_array, radius_values, R_extraction
    )
    for azi_index, azi in enumerate(azimuth_values):
        superimposed_strain = superimpose_modes_from_angle(
            colat, azi, num_time_states, mode_data
        )
        strain_to_mesh[azi_index] = np.interp(
            interpolation_times, time_array, superimposed_strain
        )

    strain_array, grid, points = initialize_tvtk_grid(
        num_azi_points, num_radius_points
    )

    # Create Mayavi objects
    engine = Engine()
    engine.start()
    engine.new_scene()
    engine.scenes[0].scene.jpeg_quality = 100
    #mlab.options.offscreen = True 
    create_gw(engine, grid, gw_color)
    bh1 = create_sphere(engine, bh1_mass * bh_scaling_factor, bh_color)
    bh2 = create_sphere(engine, bh2_mass * bh_scaling_factor, bh_color)
    mlab.view(azimuth=60, elevation = 50, distance = 80 , focalpoint=(0,0,0)) # Initialize viewpoint
    
    start_time = time.time()
    percentage = np.round(np.linspace(0, num_time_states/save_rate, 101)).astype(int)

    @mlab.animate()#ui=False) This doesn't work for some reason?
    def anim():
        for state, t, row in zip(range(num_time_states), time_array, bh_data):
            if state % save_rate != 0:
                continue

            # Print status messages
            if state == 10:
                end_time = time.time()
                eta = (end_time - start_time) * num_time_states / 10
                print(
                    f"""Creating {num_time_states} frames and saving them to: 
                    {movie_file_path}\nEstimated time: {dhms_time(eta)}"""
                )
            if status_messages and state != 0 and np.isin(state, percentage):
                print(f" {int(state * 100 / (num_time_states - 1))}% done", end="\r")

            # Change the position of the black holes
            _ = float(row[0]) # time
            bh1_xyz = (float(row[1]), float(row[2]), float(row[3]))
            bh2_xyz = (float(row[4]), float(row[5]), float(row[6]))
            change_object_position(bh1, bh1_xyz)
            change_object_position(bh2, bh2_xyz)

            points.reset()
            index = 0
            for j, radius in enumerate(radius_values):
                for i, azi in enumerate(azimuth_values):
                    h_tR = strain_to_mesh[i][state][j].real
                    x = x_values[j, i]
                    y = y_values[j, i]
                    if (
                        radius <= omitted_radius_length
                    ):  # Introduce a discontinuity to make room for the Black Holes
                        z = np.nan
                        strain_value = np.nan
                    else:
                        strain_value = h_tR * amplitude_scale_factor
                        z = strain_value 
                    points.insert_next_point(x, y, z)
                    strain_array.set_tuple1(index, strain_value)
                    index += 1
            grid._set_points(points)
            grid._get_point_data().add_array(strain_array)
            grid.modified()

            mlab.view(
                azimuth= min(60 + state*0.018, 78), 
                elevation = max(50 - state*0.016, 34),
                distance = min(80 + state*0.35, 430) 
            )

            # Save the frame
            frame_path = os.path.join(movie_file_path, f"{movie_file_name}{state:05d}.png")
            mlab.savefig(frame_path)

            if state == num_time_states - 1:
                total_time = time.time() - start_time
                print(f" {int(state * 100 / (num_time_states - 1))}% done", end="\r")
                print(f"\nSaved {num_time_states} frames to {movie_file_path} in {dhms_time(total_time)}.")
                exit(0)
            yield
    
    anim()
    mlab.show()    
                

if __name__ == "__main__":
    main()
