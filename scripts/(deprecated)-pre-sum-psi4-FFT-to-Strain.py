import os, scipy, sys, quaternionic, spherical
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

def read_modes_data(data_dir, ell_max):
    if not os.path.isdir(data_dir):
        print(f"Directory {data_dir} does not exist. Exiting...")
        exit()

    unique_times = []
    mode_data = {}

    for l in range(2, ell_max + 1):
        file_name = f"{data_dir}/Rpsi4_l{l}-r0100.0.txt"

        if not os.path.isfile(file_name):
            print(f"File {file_name} not found. Exiting...")
            exit()
        with open(file_name, "r") as file:
            lines = [line for line in file.readlines() if not line.startswith("#")]

        data = np.array([list(map(np.float64, line.split())) for line in lines])
        #data2 = np.array([np.array(line.split(), dtype=np.float64) for line in lines])
        if len(unique_times) == 0:
            unique_times = np.unique(data[:, 0])

        data = data[np.argsort(data[:, 0])]
        time_data = data[:,0]

        if time_data.all() != unique_times.all():
            print("Error: time_data discrepency. Exiting...")
            exit()

        for m in range(-l, l + 1):
            idx = 1 + 2 * (m + l)  # Calculates the column index of the real-valued data
            mode_data[(l, m)] = data[:, idx], data[:, idx + 1]

    return unique_times, mode_data

def derivative_finite_difference(time, data): #only use on l=m=2
    dt = time[1] - time[0]
    derivative = np.zeros_like(data)
    # Second-order in the interior
    derivative[1:-1] = (data[2:] - data[:-2]) / (2 * dt)
    # Drop to first-order at the endpoints
    derivative[0] = (data[1] - data[0]) / dt
    derivative[-1] = (data[-1] - data[-2]) / dt
    return derivative

def cumulative_phase_data(time, cmplx): #only use on l=m=2 
    phases = np.angle(cmplx)

    cycles = 0
    last_phase = phases[0]
    for index, phase in enumerate(phases):
        #Identify phase wrapping
        if phase - last_phase > np.pi:
            cycles -= 1
        if phase - last_phase < -np.pi:
            cycles += 1
        phases[index] += 2 * np.pi * cycles
        last_phase = phase

    phases_dt = derivative_finite_difference(time, phases)

    return  np.abs(cmplx), phases, phases_dt

def frequency_cutoff_quad_fit(time, omega, t_start=200, t_end=400): # Revisit/Room to investigate and improve. 
    def quadratic(x, a, b, c):
        return a * x**2 + b * x + c

    # Filter the data, note: non-pysical data at t=0 to t=120
    time_filtered = time[(t_start <= time) & (time <= t_end)]
    omega_filtered = omega[(t_start <= time) & (time <= t_end)]

    # Fit a quadratic curve to the omega data using nonlinear least squares
    params, _ = scipy.optimize.curve_fit(quadratic, time_filtered, omega_filtered)
    a, b, c = params

    # Find the extremum value of the quadratic curve and estimate omega(0)
    extremum_x = -b / (2 * a)
    extremum_y = np.abs(quadratic(extremum_x, a, b, c))
    omega_at_time_zero = np.abs(c)

    #print(f"The extremum of the quadratic curve occurs at t = {extremum_x:.15f} with omega = {extremum_y:.15f} . implied omega(t=0) = {omega_at_time_zero:.15f}")
    print(f"c = {omega_at_time_zero:.4g}")
    return omega_at_time_zero

def psi4_fft_to_strain(input_dir="r100", output_dir="new", num_azi_pts=180, colatitude=np.pi/2, s=-2, ell_max=8):
    if len(sys.argv) != 3:
        print("Usage: python3 script.py <path to gw psi dir> <path to output directory>")
        #sys.exit(1)
        input_dir = "r100" # Defaults
        output_dir = "new"
    else:
        input_dir = sys.argv[1]
        output_dir = sys.argv[2]

    azimuth_values = np.linspace(0, 2 * np.pi, num_azi_pts, endpoint = False)

    D = spherical.Wigner(ell_max)
    sY =  D.sYlm(s, quaternionic.array.from_spherical_coordinates(colatitude, azimuth_values))

    time, mode_data = read_modes_data(input_dir,ell_max)
    strain = np.empty((len(time), num_azi_pts),dtype=object)
    weighted = np.empty((len(time), num_azi_pts),dtype=object)
    #_, _, phase_dt = cumulative_phase_data(time, [r*i for r,i in mode_data[(2,2)]])
    omega_0 = 0.05 #frequency_cutoff_quad_fit(time, phase_dt,t_end=400)

    plt.figure(figsize=(10, 6))
    plt.grid(True)

    for idx_azi, azi in enumerate(azimuth_values):
        summation = 0
        for l in range(2, ell_max + 1):
            for m in range(-l, l + 1):
                # if m == 0 :
                #     continue
                real, imag = mode_data[(l, m)]
                cmplx = real + 1j * imag
                swsh = sY[idx_azi][D.Yindex(l,m)]
                summation += cmplx * swsh
        weighted[:,idx_azi] = summation
        # Fixed Frequency Integration
        fft_ang_freqs = np.fft.fftfreq(len(time), d = time[1] - time[0]) * 2 * np.pi
        fft_result = np.fft.fft(summation)

        # Just below Eq. 27 in https://arxiv.org/abs/1006.1632
        for idx_omg, omega in enumerate(fft_ang_freqs):
            if np.abs(omega) <= omega_0:
                fft_result[idx_omg] *= 1 / (1j * omega_0)**2
            else:
                fft_result[idx_omg] *= 1 / (1j * np.abs(omega))**2

        strain[:,idx_azi] = np.fft.ifft(fft_result)

        output_file = f"{output_dir}/Strain_azi_{azi/np.pi:.2f}pi_r0100.txt"
        with open(output_file, "w") as file:
            file.write("# Time    Strain_Real    Strain_Imag\n")
            for t, real, imag in zip(time, np.real(strain[:,idx_azi]), np.imag(strain[:,idx_azi])):
                file.write(f"{t:.16g} {real:.16g} {imag:.16g}\n")
        #print(f"Second time integral data has been saved to {output_file}")

    #idx_azi = 0
    #plt.plot(time, [np.real(x) for x in weighted[:,idx_azi] ])
    #plt.plot(time, [np.real(x) for x in strain[:,idx_azi] ])
    #plt.show()
    return strain

if __name__ == "__main__":
    psi4_fft_to_strain()