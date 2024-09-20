import numpy as np
from psychopy import core, visual
from psychopy.visual import filters


def update_stim_flicker(stim, flicker_function, flicker_frequency=15):
    animation_distance = core.monotonicClock.getTime()
    stim.color = flicker_function(animation_distance, flicker_frequency)
    stim.draw()


def create_checkerboard_texture(width, height, square_size):
    result = np.zeros((width, height))
    for y in range(height):
        for x in range(width):
            result[x][y] = (((x // square_size) + (y // square_size)) % 2 - 0.5) * 2
    return result


def create_plaid_texture(width, height, square_size):
    result = np.zeros((width, height))
    for y in range(height):
        for x in range(width):
            result[x][y] = (((x // square_size) % 2) + ((y // square_size) % 2)) - 1
    return result


def create_textured_dot_stim(win, radius, texture_function, square_size):
    # convert radius to pixels if win units is height
    diameter = radius * 2
    if win.units == "height":
        pixel_width = int(win.size[1] * radius * 2)
        square_size = int(win.size[1] * square_size)
    else:
        pixel_width = diameter
    texture_matrix = texture_function(
        width=pixel_width, height=pixel_width, square_size=square_size
    )
    dot_mask = filters.makeMask(pixel_width)
    return visual.GratingStim(
        win, tex=texture_matrix, mask=dot_mask, size=diameter, sf=1, interpolate=True
    )


def generate_grid_positions(cols, rows, padding):
    """
    :return: A list of (x, y) normalised coordinate tuples
    """
    # the 1e-10 allows us to get an inclusive range
    return [
        (x, y)
        for x in np.arange(
            -0.5 + padding[0],
            0.5 - padding[0] + 1e-10,
            (1 - (2 * padding[0])) / (rows - 1),
        )
        for y in np.arange(
            -0.5 + padding[1],
            0.5 - padding[1] + 1e-10,
            (1 - (2 * padding[1])) / (cols - 1),
        )
    ]
    

def flicker_sin(distance, frequency):
    return np.sin(distance * 2 * np.pi * frequency)


def flicker_square(distance, frequency):
    return (int(distance * 2 * frequency) % 2) * 2 - 1
