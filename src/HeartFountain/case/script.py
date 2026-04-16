# Copyright (c) 2026 JackLee
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT


from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def generateFileHeader():
    return """
FoamFile
{
    version     2.0;
    format      ascii;
    class       Tuple2VectorVector;
    object      positions;
}
"""


def writePositions(filename: str, points: np.ndarray):
    with open(filename, "w") as f:
        f.write(generateFileHeader())
        f.write(f"{len(points)}\n(\n")
        for point in points:
            f.write(
                f"(({point[0]} {point[1]} {point[2]}) ({point[3]} {point[4]} {point[5]}))\n"
            )
        f.write(")\n")


def heartCurve(t):
    x = 16 * np.sin(t) ** 3
    y = 13 * np.cos(t) - 5 * np.cos(2 * t) - 2 * np.cos(3 * t) - np.cos(4 * t)
    return x, y


def dHeartCurve(t):
    dx = 48 * np.sin(t) ** 2 * np.cos(t)
    dy = -13 * np.sin(t) + 10 * np.sin(2 * t) + 6 * np.sin(3 * t) + 4 * np.sin(4 * t)
    return dx, dy


if __name__ == "__main__":
    target_file = Path("./constant/positionsFile")
    nPoints = 128
    scale = 0.001
    offset_x = 0.05
    offset_y = 0.05
    z = 0.005

    t = np.linspace(0, 2 * np.pi, nPoints, endpoint=False)
    x, y = heartCurve(t)
    points = []

    with np.errstate(divide="ignore", invalid="ignore"):
        for i in range(len(x)):
            dx, dy = dHeartCurve(t[i])
            dx /= np.linalg.norm([dx, dy])  # Normalize the tangent vector
            dy /= np.linalg.norm([dx, dy])
            direction = np.array([-dy, dx, 2])
            direction /= np.linalg.norm(direction)  # Normalize the direction vector
            points.append(
                (
                    x[i] * scale + offset_x,
                    y[i] * scale + offset_y,
                    z,
                    direction[0],
                    direction[1],
                    direction[2],
                )
            )

    points = np.array(points)
    writePositions(target_file, points[~np.isnan(points).any(axis=1)])

    plt.plot(points[:, 0], points[:, 1], "ro")
    plt.axis("equal")
    plt.show()
