import os, math, cmath, time, numpy, scipy
import quaternionic, spherical
def old_sph_harm(s, l, m, azi, colat): #uses scipy.special sph_harm
    sign = (-1) ** s
    root_factor = math.sqrt((2 * l + 1) * math.factorial(l - m) / math.factorial(l + m) / (4 * math.pi)) # why?
    sYlm = sign * root_factor * scipy.special.sph_harm(m, l, azi, colat)
    return sYlm

def spin_weight_sph_harm2(s, l, m, azi, colat):
    sign = (-1)**(m)
    root_factor = math.sqrt(math.factorial(l + m) * math.factorial(l - m) * (2*l + 1)) / (4*math.pi * math.factorial(l + s) * math.factorial(l - s))
    phase = cmath.exp(1j * m * azi)
    legendre_term = scipy.special.lpmv(m, l, math.cos(colat))
    sYlm = sign * root_factor * phase * legendre_term
    return sYlm

def wiki_sph_harm(s, l, m, azi, colat):
    sign = (-1)**(l + m - s)
    root_factor = math.sqrt((math.factorial(l+m) * math.factorial(l-m) * (2*l + 1)) / (4 * math.pi * math.factorial(l+s) * math.factorial(l-s)))
    sine_and_phase = math.sin(colat/2)**(2*l) * cmath.exp(1j * m * azi)
    sum = 0
    for r in range(0, l - s + 1):
        sum += (-1)**r * math.comb(l - s, r) * math.comb(l + s, r + s - m) * ((math.cos(colat/2) / math.sin(colat/2))**(2*r + s - m))
    sYlm = sign * root_factor * sine_and_phase * sum
    return sYlm

def main():
    # Plot parameters
    numRadius = 450  # Number of points in the radial direction
    numAzi = 180  # Number of points per circle
    display_radius = 300  # Mesh radius
    R_ext = 100  # Extraction radius
    s = 0
    l = 4
    m = 1
    scale_factor = 600
    omit_radius = 4

    azi_values = [math.pi/2,] # numpy.linspace(0, 2 * math.pi, numAzi, endpoint=False)
    #old_swsh = old_sph_harm(s, l, m, math.pi/4, math.pi/2)
    new_swsh = wiki_sph_harm(s, l, m, math.pi/4, math.pi/2)
    #print(old_swsh)
    print(new_swsh)

if __name__ == "__main__":
    main()
    