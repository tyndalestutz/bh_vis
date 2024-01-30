import os, math, cmath, time, numpy, scipy
import quaternionic, spherical

def wiki_sph_harm(s, l, m, colat, azi):
    sign = (-1)**(l + m - s)
    root_factor = math.sqrt( (math.factorial(l+m)*math.factorial(l-m)*(2*l + 1))/(4*math.pi*math.factorial(l+s)*math.factorial(l-s)) )
    sine_and_phase = numpy.sin(colat/2)**(2*l) * numpy.exp(1j * m * azi)
    sum = 0
    for r in range(0, l-s+1):
        sum += (-1)**r * scipy.special.comb(l-s, r) * scipy.special.comb(l+s, r+s-m) * (numpy.cos(colat/2)/numpy.sin(colat/2))**(2*r+s-m)
    sYlm = sign * root_factor * sine_and_phase * sum
    return sYlm

def main():
    numAzi = 180
    s = -2 # -l <= s <= l
    l = 4
    m = 4 # -1 <= m <= l

    colat_values = math.pi/6 #numpy.oneslike(azi_values) * math.pi / 2
    azi_values = numpy.linspace(0, 2 * math.pi, numAzi, endpoint = False)

    ell_max = 8
    winger = spherical.Wigner(ell_max)
    
    wiki_swsh = wiki_sph_harm(s, l, m, colat_values, azi_values)
    sY = winger.sYlm(s, quaternionic.array.from_spherical_coordinates(colat_values, azi_values))
    
    for i in range(numAzi):
        if(abs(sY[i][winger.Yindex(l , m)] - wiki_swsh[i])<=0.00000001):
            print(i)
        else:
            print(abs(sY[i][winger.Yindex(l , m)] - wiki_swsh[i]))

if __name__ == "__main__":
    main()