#!/usr/bin/env python
import json, sys
import numpy as np

sys.exit("Do not use this script since it will remove the hand-tuned starting values for fitting")

def process_data(filename_prefix, filename_prefix_esp, skiprows):
    grid_data = np.loadtxt(f'ESP_Alkali_Halides/{filename_prefix}000.dat', skiprows=skiprows)
    grid_esp_data = np.loadtxt(f'ESP_Alkali_Halides/grid_esp-{filename_prefix_esp}000.dat', skiprows=skiprows)
    distances = np.sqrt(np.sum(grid_data ** 2, axis=1))
    esp_potentials = grid_esp_data * 2625.4996
    return distances, esp_potentials

bounds = {
    'F': [
        ([0., .1], [26, 30]),
        ([.01, .1], [26, 4]),
        ([.01, -10, .1, .1], [10, -0.01, 30, 30]),
        ([.01, -10, 1.0, 1], [26, -0.01, 30, 40]),
    ],
    'Cl': [
        ([1, .1], [16, 30]),
        ([.01, .1], [30, 4]),
        ([.1, -10, .10, .1], [30, -0.01, 30, 10]),
        ([.1, -10, 1., 1], [50, -0.01, 30, 40])
    ],
    'Br': [
        ([1, .1], [50, 10]),
        ([.01, .001], [100, 10]),
        ([.1, -10, 1, .1], [40, -0.01, 40, 40]),
        ([.1, -20, .1, .1], [50, -0.01, 30, 40])
    ],
    'Li': [
        ([.1, .1], [26, 5]),
        ([.1, .1], [26, 8]),
        ([.01, -10, .1, .1], [56, -0.01, 100, 100]),
        ([.01, -10, .1, .1], [26, -0.01, 40, 40])
    ],
    'Na': [
        ([.1, .1], [26, 10]),
        ([.1, .1], [26, 40]),
        ([.1, -10, .1, .1], [26, -0.01, 60, 60]),
        ([.1, -10, .1, .1], [26, -0.01, 100, 100])
    ],
    'K': [
        ([.1, .1], [26, 10]),
        ([.1, .1], [26, 40]),
        ([.1, -10, .1, .1], [60, -0.01, 40, 40]),
        ([.1, -10, 1.0, 1], [60, -0.01, 40, 40]),
    ]
}

initial_guesses = {
    'F': [
        [1, .1],
        [.01, .1],
        [.1, -1.10, .1, .1],
        [.1, -.10, 1, 1]
    ],
    'Cl': [
        [1, .1],
        [.01, .1],
        [.1, -1.10, .1, .1],
        [.1, -.10, 1, 1]
    ],
    'Br': [
        [1, .1],
        [.01, .001],
        [.1, -1.10, 1, .1],
        [.1, -10.1, 1, .1]
    ],
    'Li': [
        [.1, .1],
        [.1, .1],
        [.1, -.10, .1, .1],
        [.1, -.10, .1, .1]
    ],
    'Na': [
        [.1, .1],
        [.1, .1],
        [.1, -.10, 1, .1],
        [.1, -.10, .1, .1]
    ],
    'K': [
        [.1, .1],
        [.1, .1],
        [.1, -.10, .1, .1],
        [.1, -.10, 1, 1]
    ]
}

for T, N in [(100, 0), (10, 200)]:
    data = {}
    for elem in ["Br", "F", "Cl", "Li", "Na", "K"]:
        dist, pot = process_data('grid', f'hf-{elem.lower()}', N)
        data[elem] = [dist.tolist(), pot.tolist()]

    output_data = {
        "charge": -1,
        "data": data,
        "initial_guesses": initial_guesses,
        "bounds": bounds
    }

    with open(f'output_4_{T}.json', 'w') as json_file:
        json.dump(output_data, json_file, indent=4)
