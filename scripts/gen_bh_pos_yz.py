import os
import numpy as np
from numpy.typing import NDArray
file_name="puncture_posns_vels_regridxyzU.txt"
file_path = os.path.abspath(os.path.join(__file__, '..', '..', "r100", file_name))
with open(file_path, mode="r", encoding="utf-8") as file:
    lines = [line for line in file.readlines() if not line.startswith("#")]
    data: NDArray[np.float64] = np.array(
            [list(map(np.float64, line.split())) for line in lines]
        )
    data = data[np.argsort(data[:, 0])]
print(data)
