from eyelink_wrapper import EyelinkWrapper
import numpy as np
import os
from psychopy import visual, event, core
from psychopy.visual import filters
import random
from experiment_utils import ExperimentLogCSV, prompt_for_params
from stim_utils import (
    create_textured_dot_stim,
    create_checkerboard_texture,
    flicker_sin,
    generate_grid_positions,
    update_stim_flicker,
)

experiment_params = {
    "flicker_frequency": 8,
    "flicker_duration": 10,  # 10s ON and the last 20 seconds OFF
    "dot_radius": 0.14,  # 0.13 correspons to a diameter of 9 cm [measured]; 0.14 should correspond to 10 cm, VA = 5° 
    "dot_checker_size": 0.015,
    "grid_rows": 3,
    "grid_cols": 3,
    "grid_padding": (0, 0.2),
    "repetitions": 10,
    '"foo': '"bar;";',
}

EXPERIMENT_NAME = "Bold"
SCREEN_SIZE = (800, 600)
PROGRESS_SIGNAL = "b"
SYNCHRONIZE_SIGNAL = "s"
RESULTS_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), "results"))


def wait_for_key(key):
    while not event.getKeys(keyList=[key]):
        core.wait(0.1)
        if event.getKeys(keyList=["escape"]):
            raise Exception("Experiment sequence terminated with esc key")


def run_fixation_sequence(win, logger):
    logger.log("Start fixation")
    visual.TextStim(win, text="+", height=0.04).draw()
    win.flip()
    wait_for_key(PROGRESS_SIGNAL)
    logger.log("End fixation")
    win.flip()


def update_dot_stim_flicker(stim, flicker_function, flicker_frequency):
    animation_distance = core.monotonicClock.getTime()
    stim.color = flicker_function(animation_distance, flicker_frequency)
    stim.draw()


def run_stimulus_sequence(
    win,
    dot_radius,
    dot_checker_size,
    grid_rows,
    grid_cols,
    grid_padding,
    flicker_frequency,
    flicker_duration,
    repetitions,
    logger,
    **params
):
    # Initialize stimulus
    checkered_dot = create_textured_dot_stim(
        win, dot_radius, create_checkerboard_texture, dot_checker_size
    )
    grid_positions = generate_grid_positions(
        rows=grid_rows, cols=grid_cols, padding=grid_padding
    )
    #print(grid_positions)
    [grid_positions.pop(i) for i in [8,6,5,4,3,2,0]]
    print(grid_positions)
    #print(grid_positions)
    #print(jfk)
    random.shuffle(grid_positions)
    flicker_duration_timer = core.Clock()
    wait_for_key(SYNCHRONIZE_SIGNAL)
    # Initialize a generator of ordered updates to stim params
    for position in grid_positions * repetitions:
        logger.log("Present stimulus", *position)
        checkered_dot.pos = position
        # Play flciker animation
        flicker_duration_timer.reset()
        while flicker_duration_timer.getTime() < flicker_duration:
            update_stim_flicker(checkered_dot, flicker_sin, flicker_frequency)
            win.flip()
            core.wait(0.017)
            if event.getKeys(keyList=["escape"]):
                raise Exception("Experiment sequence terminated with esc key")
        # Set stim to middle grey
        checkered_dot.color = 0
        checkered_dot.draw()
        visual.TextStim(win, text="+",height=0.04, pos=(position[0],position[1]), alignVert='center',color=[-1,-1,-1]).draw()
        win.flip()
        # Wait for control signal
        wait_for_key(SYNCHRONIZE_SIGNAL)


def run_experiment(use_eyelink):
    params = {
        **experiment_params,
        **prompt_for_params(("subject_id", "Subject id:"), title=EXPERIMENT_NAME),
    }
    logger = ExperimentLogCSV(EXPERIMENT_NAME, params)
    logger.start()
    output_file_name = "%s_%s.edf" % (EXPERIMENT_NAME, params["subject_id"])
    if use_eyelink:
        eyelink = EyelinkWrapper(
            screen_size=SCREEN_SIZE, output_file_name=output_file_name
        )
    print('#############################')
    print(output_file_name)
    print('#############################')
    if use_eyelink:
        eyelink.start_recording()
    win = visual.Window(
        size=SCREEN_SIZE, fullscr=True, color=[-1, -1, -1], units="height"
    )
    try:
        logger.log("Setup complete")
        wait_for_key(PROGRESS_SIGNAL)
        run_fixation_sequence(win, logger)
        run_stimulus_sequence(win, logger=logger, **params)
        logger.log("Experiment complete")
    except Exception as error:
        print("EXITING WITH ERROR:", error)
        logger.log(repr(error))
    finally:
        logger.finish()
        if use_eyelink:
            eyelink.finish_recording(RESULTS_DIR)


run_experiment(use_eyelink=True)
