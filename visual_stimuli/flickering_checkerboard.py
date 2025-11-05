#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment display a parameterizable flickering checkerboard stimuli
and track eye gaze.
"""

# --- Import packages ---
from random import randint, choice, uniform
import numpy as np
import os
import socket
from psychopy import visual, hardware, gui, data, core, __version__
from psychopy.iohub.client.connect import launchHubServer
from psychopy.constants import NOT_STARTED, STARTED, FINISHED
from psychopy import prefs, plugins, logging
from psychopy.hardware import keyboard
from psychopy.monitors import Monitor
from typing import Tuple as tuple

plugins.activatePlugins()
prefs.hardware['audioLib'] = 'ptb'
prefs.hardware['audioLatencyMode'] = '3'


# --- Parameters ---
DUMMY_MODE: bool = True  # TODO: Set False when using real eye-tracker
FULL_SCREEN: bool = True
MOUSE_VISIBLE: bool = False
FRAME_TOLERANCE: float = 0.001  # s
PROGRESS_SIGNAL: str = "s"
# SYNCHRONIZE_SIGNAL: str = "b"
WIN_SIZE: tuple[int, int] = (800, 600)  # is probably the scanner room's screen size
ON_DURATION: float = 1.0  # Duration in seconds of the checkerboard flickering, the OFF time is determined by the scanner's trigger box
OFF_DURATION: float = 4.0
NUM_REPETITIONS: int = int((16*60 + 35)/(ON_DURATION + OFF_DURATION))  # 252
FREQUENCIES: tuple[float, ...] = (8,)  # , 8, 10, 12, 14)  # Hz
FIXATION_SHAPE: str = "circle_cross_circle"
assert FIXATION_SHAPE in ("circle", "cross", "circle_cross_circle", "peripheral_circles")
FIXATION_COLORS: tuple[str, ...] = ("red", "green", "blue")
FIXATION_CHANGE_DELAY: tuple[float, float] = (1.0, 2.0)

FLICKER_SHAPE: str = "square"
assert FLICKER_SHAPE in ("circle", "ring", "square", "wedge", "windmill")  # Add rotating wedge ? Expanding/Contracting ring ?

num_textures: int = 0
if FLICKER_SHAPE == "square":
    NUM_SQUARES: tuple[int, ...] = (8, )
    num_textures: int = len(NUM_SQUARES)

elif FLICKER_SHAPE == "circle":
    OUTER_FLICKER_RADIUS: tuple[int, ...] = (200, 400)  # in pixels
    INNER_FLICKER_RADIUS: tuple[int, ...] = (0, 200)  # in pixels
    NUM_RINGS: tuple[int, ...] =            (6,   6, )
    NUM_SECTORS: tuple[int, ...] =          (16,  16,)
    assert len(set(map(len, (NUM_RINGS, NUM_SECTORS, OUTER_FLICKER_RADIUS, INNER_FLICKER_RADIUS)))) == 1  # All tuple of parameter must have the same length
    num_textures: int = len(NUM_RINGS)

elif FLICKER_SHAPE == "ring":
    OUTER_FLICKER_RADIUS: tuple[int, ...] = (600, 600)  # in pixels
    RING_THICKNESS: tuple[int, ...] = (200, 300)  # in pixels
    EXPANSION_SPEED: tuple[float, ...] = (10.0, 20.0)
    NUM_RINGS: tuple[int, ...] = (10, 10,)
    NUM_SECTORS: tuple[int, ...] = (16, 16,)
    assert len(set(map(len, (NUM_RINGS, NUM_SECTORS, OUTER_FLICKER_RADIUS, RING_THICKNESS, EXPANSION_SPEED)))) == 1  # All tuple of parameter must have the same length
    num_textures: int = len(NUM_RINGS)

elif FLICKER_SHAPE == "wedge":
    OUTER_FLICKER_RADIUS: tuple[int, ...] = (400, 400)  # in pixels
    INNER_FLICKER_RADIUS: tuple[int, ...] = (0, 0)  # in pixels
    NUM_RINGS: tuple[int, ...] =            (6, 6)
    NUM_SECTORS: tuple[int, ...] =          (6, 6)
    START_ANGLE: tuple[float, ...] = (90, 270)
    DELTA_ANGLE: tuple[float, ...] = (180, 180)
    ROTATION_SPEED: tuple[float, ...] = (0, 0)  # °/s
    assert len(set(map(len, (NUM_RINGS, NUM_SECTORS, OUTER_FLICKER_RADIUS, INNER_FLICKER_RADIUS, START_ANGLE, DELTA_ANGLE, ROTATION_SPEED)))) == 1  # All tuple of parameter must have the same length
    num_textures: int = len(NUM_RINGS)

elif FLICKER_SHAPE == "windmill":
    OUTER_FLICKER_RADIUS: tuple[int, ...] = (400, 400)  # in pixels
    INNER_FLICKER_RADIUS: tuple[int, ...] = (0, 0)  # in pixels
    NUM_RINGS: tuple[int, ...] = (6, 6)
    NUM_SECTORS: tuple[int, ...] = (6, 6)
    NUM_WEDGES: tuple[int, ...] = (4, 8)
    ROTATION_SPEED: tuple[float, ...] = tuple((360/(n_wedge))/ON_DURATION for n_wedge in NUM_WEDGES)  # °/s
    assert len(set(map(len, (NUM_RINGS, NUM_SECTORS, OUTER_FLICKER_RADIUS, INNER_FLICKER_RADIUS, NUM_WEDGES, ROTATION_SPEED)))) == 1  # All tuple of parameter must have the same length
    num_textures: int = len(NUM_RINGS)


routine_duration: float = (ON_DURATION + OFF_DURATION) * num_textures * len(FREQUENCIES) * NUM_REPETITIONS

# --- Global variables ---
end_exp_now: bool = False

print(f"Routine duration = {int(routine_duration // 60)}:{int(routine_duration % 60):02d}")
print(f"With {NUM_REPETITIONS} stimuli periods.")


def send_message(message, addr: str = "localhost", port: int = 2023) -> None:
    client_socket: socket.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    client_socket.connect((addr, port))
    client_socket.sendall(message)
    client_socket.close()


def make_ring_checkerboard(size: int,
                           n_radial: int,
                           n_angular: int,
                           outer_radius: float,
                           inner_radius: float = 0) -> np.ndarray:
    """
    Create a ring (annular) checkerboard image.

    :param size: Size in pixels of the square image (size x size)
    :param outer_radius: Outer radius of the ring (in pixels)
    :param n_radial: Number of radial divisions (rings)
    :param n_angular: Number of angular divisions (wedges)
    :param inner_radius: Inner radius of the ring (in pixels)
    :return: 2D numpy array with values in [-1, 1], 0 outside the ring
    """
    if inner_radius >= outer_radius:
        raise ValueError("inner_radius must be smaller than outer_radius")

    y, x = np.indices((size, size))
    cx, cy = size // 2, size // 2
    r = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
    theta = np.arctan2(y - cy, x - cx)
    theta[theta < 0] += 2 * np.pi  # Map to [0, 2π]

    # Normalize radius to [0, 1] in the annulus region
    r_norm = (r - inner_radius) / (outer_radius - inner_radius)

    # Checkerboard pattern only in the ring
    radial_idx = np.floor(r_norm * n_radial).astype(int)
    angular_idx = np.floor(theta / (2 * np.pi / n_angular)).astype(int)

    checker = ((radial_idx + angular_idx) % 2) * 2 - 1  # values in {-1, 1}

    # Mask outside the ring
    mask = (r >= inner_radius) & (r <= outer_radius)
    checker[~mask] = 0

    return checker


def annulus_mask(size: int,
                 outer_radius: float,
                 ring_thickness: float,
                 expansion_speed: float,
                 t: float):
    """Return annulus mask in [-1, 1] for continuous expanding ring."""
    # Moving inner/outer edges
    inner_r = (expansion_speed * t) % outer_radius
    outer_r = inner_r + ring_thickness

    y, x = np.indices((size, size))
    cx, cy = size // 2, size // 2
    r = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)

    if outer_r <= outer_radius:
        # Normal case: simple annulus
        mask = (r >= inner_r) & (r <= outer_r)
    else:
        # Wrap case: ring is split across center and edge
        wrap_outer = outer_r - outer_radius
        mask = ((r >= inner_r) & (r <= outer_radius)) | (r <= wrap_outer)

    return mask.astype(float) * 2 - 1  # PsychoPy expects [-1, 1]


def make_wedge_checkerboard(size: int,
                            n_radial: int,
                            n_angular: int,
                            outer_radius: int,
                            inner_radius: int = 0,
                            start_angle: float = 0.0,
                            end_angle: float = 90.0) -> np.ndarray:
    """
    Create a wedge-shaped checkerboard (section of a disc) with correct number of angular sectors.
    """
    if inner_radius >= outer_radius:
        raise ValueError("inner_radius must be smaller than outer_radius")

    y, x = np.indices((size, size))
    cx, cy = size // 2, size // 2
    r = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
    theta = np.arctan2(y - cy, x - cx)
    theta[theta < 0] += 2 * np.pi  # Map to [0, 2π]

    # Convert angles to radians
    start_rad = np.deg2rad(start_angle)
    end_rad = np.deg2rad(end_angle)
    wedge_rad = end_rad - start_rad
    if wedge_rad <= 0:
        wedge_rad += 2 * np.pi  # wrap-around

    # Radial index
    r_norm = (r - inner_radius) / (outer_radius - inner_radius)
    radial_idx = np.floor(r_norm * n_radial).astype(int)

    # Angular index scaled to the wedge
    angular_idx = np.floor((theta - start_rad) / wedge_rad * n_angular).astype(int)

    # Checker pattern
    checker = ((radial_idx + angular_idx) % 2) * 2 - 1

    # Mask outside the wedge
    mask_radius = (r >= inner_radius) & (r <= outer_radius)
    mask_angle = (theta >= start_rad) & (theta <= end_rad) if start_rad < end_rad else (theta >= start_rad) | (
                theta <= end_rad)
    mask = mask_radius & mask_angle
    checker[~mask] = 0

    return checker


def make_windmill_checkerboard(size: int,
                               n_radial: int,
                               n_angular: int,
                               n_wedges: int,
                               outer_radius: int,
                               inner_radius: int = 0) -> np.ndarray:
    """
    Create a windmill-shaped checkerboard pattern with both radial and angular alternation.

    Parameters
    ----------
    size : int
        Image size (square, pixels).
    n_radial : int
        Number of radial divisions in the checkerboard.
    n_angular : int
        Number of angular divisions per wedge (checkerboard alternation around circle).
    n_wedges : int
        Total number of wedges (blades in the windmill).
        Half will be checkerboard, half empty (gray=0).
    outer_radius : int
        Outer radius of the windmill (pixels).
    inner_radius : int, optional
        Inner radius (default = 0).
    """
    if inner_radius >= outer_radius:
        raise ValueError("inner_radius must be smaller than outer_radius")

    n_wedges *= 2

    # Coordinates
    y, x = np.indices((size, size))
    cx, cy = size // 2, size // 2
    r = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
    theta = np.arctan2(y - cy, x - cx)
    theta[theta < 0] += 2 * np.pi  # [0, 2π]

    # Angular wedge index
    angular_wedge_idx = np.floor(theta / (2 * np.pi) * n_wedges).astype(int)

    # Normalized radial index
    r_norm = (r - inner_radius) / (outer_radius - inner_radius)
    radial_idx = np.floor(r_norm * n_radial).astype(int)

    # Within each wedge, create angular checker index
    wedge_angle = 2 * np.pi / n_wedges
    angle_within_wedge = theta % wedge_angle
    angular_idx = np.floor(angle_within_wedge / wedge_angle * n_angular).astype(int)

    # Checkerboard (radial + angular alternation)
    checker = ((radial_idx + angular_idx) % 2) * 2 - 1

    # Mask: only within radius
    mask_radius = (r >= inner_radius) & (r <= outer_radius)

    # Windmill selection: only even wedges are checkerboard
    windmill_mask = (angular_wedge_idx % 2 == 0)

    # Apply masks
    result = np.zeros_like(checker, dtype=float)
    result[mask_radius & windmill_mask] = checker[mask_radius & windmill_mask]

    return result


def show_exp_info_dlg(exp_info: dict) -> dict:
    # Show participant info dialog
    dlg = gui.DlgFromDict(
        dictionary=exp_info, title=f"{exp_info['expName|hid']} \n8 characters or less, only letter and number\n(no special char)",
        sortKeys=False, alwaysOnTop=True
    )
    if not dlg.OK:
        core.quit()  # If user pressed cancel

    return exp_info
    

def setup_data(exp_info: dict) -> data.ExperimentHandler:
    # data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
    data_dir: str = os.path.dirname(os.path.abspath(__file__))
    os.chdir(data_dir)
    filename: str = f"data/{exp_info['participant']}_{exp_info['expName|hid']}_{exp_info['date|hid']}"
    
    # make sure filename is relative to data_dir
    if os.path.isabs(filename):
        data_dir = os.path.commonprefix([data_dir, filename])
        filename = os.path.relpath(filename, data_dir)
    
    # an ExperimentHandler isn't essential but helps with data saving
    this_exp = data.ExperimentHandler(
        name=exp_info['expName|hid'],
        version='',
        extraInfo=exp_info,
        runtimeInfo=None,
        originPath=data_dir,
        savePickle=True,
        saveWideText=True,
        dataFileName=data_dir + os.sep + filename,
        # sortColumns='time',
    )
    # this_exp.setPriority('thisRow.t', priority.CRITICAL)
    # this_exp.setPriority('expName', priority.LOW)

    return this_exp


def setup_logging(filename: str) -> logging.LogFile:
    # set how much information should be printed to the console / app
    logging.console.setLevel(logging.WARNING)
    
    # save a log file for detail verbose info
    log_file: logging.LogFile = logging.LogFile(f"{filename}.log")
    log_file.setLevel(logging.WARNING)
    return log_file


def setup_window(exp_info: dict, full_scr: bool = True) -> visual.Window:
    monitor = Monitor("expMonitor")
    monitor.setWidth(369.54e-3)  # screen width in meters
    monitor.setDistance(1020e-3)  # viewing distance in meters
    monitor.setSizePix(WIN_SIZE)  # screen resolution
    monitor.saveMon()

    win: visual.Window = visual.Window(
        size=WIN_SIZE,
        fullscr=full_scr,
        screen=0,
        winType='pyglet',
        allowStencil=False,
        monitor="expMonitor",
        color=(0, 0, 0),
        colorSpace='rgb',
        backgroundImage='',
        backgroundFit='none',
        blendMode='avg',
        useFBO=True,
        units='deg',
        checkTiming=True,
    )

    win.mouseVisible = MOUSE_VISIBLE
    # win.hideMessage()

    # get/measure frame rate if not already in expInfo
    exp_info['frameRate'] = win.getActualFrameRate()

    # get frame duration from frame rate in expInfo
    if 'frameRate' in exp_info and exp_info['frameRate'] is not None:
        frame_dur: float = 1.0 / round(exp_info['frameRate'])
    else:
        frame_dur: float = 1.0 / 60.0  # could not measure, so guess
    
    return win


def setup_devices(exp_info: dict, win):
    # --- Setup input devices ---
    io_config: dict = {
        'eyetracker.hw.sr_research.eyelink.EyeTracker':
        {
            'name': 'tracker',
            'model_name': 'EYELINK 1000 DESKTOP',
            'simulation_mode': DUMMY_MODE,
            'network_settings': '100.1.1.1',
            'default_native_data_file_name': 'EXPFILE',
            'runtime_settings': {
                'sampling_rate': 1000.0,
                'track_eyes': 'RIGHT_EYE',
                'sample_filtering': {
                    'sample_filtering': 'FILTER_LEVEL_2',
                    'elLiveFiltering': 'FILTER_LEVEL_OFF',
                },
                'vog_settings': {
                    'pupil_measure_types': 'PUPIL_AREA',
                    'tracking_mode': 'PUPIL_CR_TRACKING',
                    'pupil_center_algorithm': 'ELLIPSE_FIT',
                }
            }
        },
        'Keyboard': dict(use_keymap='psychopy'),
    }
    
    io_session: str = "1"
    if "session" in exp_info:
        io_session: str = str(exp_info['session'])

    # Start ioHub server
    io_server = launchHubServer(window=win, **io_config)

    # Get eye tracker (from psychopy.iohub.devices.eyetracker.hw.sr_research.eyelink import EyeTracker)
    eye_tracker = io_server.getDevice("tracker")

    # Get default keyboard
    default_keyboard = keyboard.Keyboard(backend="iohub")

    return eye_tracker, default_keyboard


def run_calibration(default_keyboard, routine_timer, win, eye_tracker) -> None:
    calibration_target = visual.TargetStim(win,
                                            name='calibration_2Target',
                                            radius=0.15, fillColor=[0.5, 0.5, 0.5], borderColor=[0.5, 0.5, 0.5],
                                            lineWidth=2.0,
                                            innerRadius=0.07, innerFillColor=[0.5, 0.5, 0.5],
                                            innerBorderColor=[0.5, 0.5, 0.5], innerLineWidth=2.0,
                                            colorSpace='rgb', units=None
                                            )
    # define parameters for calibration
    calibration = hardware.eyetracker.EyetrackerCalibration(win,
                                                            eye_tracker, calibration_target,
                                                            units=None, colorSpace='rgb',
                                                            progressMode='time', targetDur=1.5, expandScale=1.5,
                                                            targetLayout='FIVE_POINTS', randomisePos=True,
                                                            textColor='white',
                                                            movementAnimation=True, targetDelay=1.0
                                                            )
    # run calibration
    calibration.run()
    # clear any keypresses from during calibration_2 so they don't interfere with the experiment
    default_keyboard.clearEvents()
    # the Routine "calibration_2" was not non-slip safe, so reset the non-slip timer
    routine_timer.reset()


def start_eye_tracker(eye_tracker, win, routine_timer, this_exp, default_keyboard) -> None:
    eye_tracker_record = hardware.eyetracker.EyetrackerControl(
        tracker=eye_tracker,
        actionType='Start Only'
    )

    # --- Prepare to start Routine "start_ET" ---
    continueRoutine = True
    # update component parameters for each repeat
    # keep track of which components have finished
    start_ETComponents = [eye_tracker_record]
    for thisComponent in start_ETComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1

    # --- Run Routine "start_ET" ---
    routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routine_timer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routine_timer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        # *etRecord* updates

        # if etRecord is starting this frame...
        if eye_tracker_record.status == NOT_STARTED and tThisFlip >= 0.0 - FRAME_TOLERANCE:
            # keep track of start time/frame for later
            eye_tracker_record.frameNStart = frameN  # exact frame index
            eye_tracker_record.tStart = t  # local t and not account for scr refresh
            eye_tracker_record.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(eye_tracker_record, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            this_exp.timestampOnFlip(win, 'etRecord.started')
            # update status
            eye_tracker_record.status = STARTED
            # Run 'Begin Routine' code from code_channel2
            eye_tracker.sendMessage("Hello tracker record")

        # if etRecord is stopping this frame...
        if eye_tracker_record.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > eye_tracker_record.tStartRefresh + 0 - FRAME_TOLERANCE:
                # keep track of stop time/frame for later
                eye_tracker_record.tStop = t  # not accounting for scr refresh
                eye_tracker_record.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                this_exp.timestampOnFlip(win, 'etRecord.stopped')
                # update status
                eye_tracker_record.status = FINISHED
                eye_tracker.sendMessage("Bye tracker record")

        # check for quit (typically the Esc key)
        if end_exp_now or default_keyboard.getKeys(keyList=["escape"]):
            eye_tracker.sendMessage("key board escape")
            core.quit()

        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in start_ETComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished

        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over, or we'll get a blank screen
            win.flip()


def stop_eye_tracker(eye_tracker, win, default_keyboard, routine_timer, this_exp) -> None:
    text = visual.TextStim(win=win, name='text',
        text="End press 't'",
        height=0.12, color='white', units="norm",
    )

    eye_tracker_stop = hardware.eyetracker.EyetrackerControl(
        tracker=eye_tracker,
        actionType='Stop Only'
    )

    # --- Prepare to start Routine "end" ---
    continueRoutine = True
    # update component parameters for each repeat
    default_keyboard.keys = []
    default_keyboard.rt = []
    _default_keyboard_allKeys = []
    # Run 'Begin Routine' code from code_3
    eye_tracker.sendMessage("ET: Prepare to start routine 'end'")
    # keep track of which components have finished
    endComponents = [text, eye_tracker_stop, default_keyboard]
    for thisComponent in endComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1

    # --- Run Routine "end" ---
    routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routine_timer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routine_timer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame

        # *text* updates

        # if text is starting this frame...
        if text.status == NOT_STARTED and tThisFlip >= 0.0 - FRAME_TOLERANCE:
            # keep track of start time/frame for later
            text.frameNStart = frameN  # exact frame index
            text.tStart = t  # local t and not account for scr refresh
            text.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            this_exp.timestampOnFlip(win, 'text.started')
            # update status
            text.status = STARTED
            text.setAutoDraw(True)

        # if text is active this frame...
        if text.status == STARTED:
            # update params
            pass

        # *ET_stop* updates
        if eye_tracker_stop.status == NOT_STARTED:
            eye_tracker_stop.frameNStart = frameN  # exact frame index
            eye_tracker_stop.tStart = t  # local t and not account for scr refresh
            eye_tracker_stop.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(eye_tracker_stop, 'tStartRefresh')  # time at next scr refresh
            eye_tracker_stop.status = STARTED

        # if ET_stop is stopping this frame...
        if eye_tracker_stop.status == STARTED:
            # is it time to stop? (based on local clock)
            if tThisFlip > 1.0 - FRAME_TOLERANCE:
                # keep track of stop time/frame for later
                eye_tracker_stop.tStop = t  # not accounting for scr refresh
                eye_tracker_stop.tStopRefresh = tThisFlipGlobal  # on global time
                eye_tracker_stop.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                this_exp.timestampOnFlip(win, 'ET_stop.stopped')
                # update status
                eye_tracker_stop.status = FINISHED
                eye_tracker_stop.stop()

        # *default_keyboard* updates
        waitOnFlip = False

        # if default_keyboard is starting this frame...
        if default_keyboard.status == NOT_STARTED and tThisFlip >= 0.0 - FRAME_TOLERANCE:
            # keep track of start time/frame for later
            default_keyboard.frameNStart = frameN  # exact frame index
            default_keyboard.tStart = t  # local t and not account for scr refresh
            default_keyboard.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(default_keyboard, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            this_exp.timestampOnFlip(win, 'default_keyboard.started')
            # update status
            default_keyboard.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(default_keyboard.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(default_keyboard.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if default_keyboard.status == STARTED and not waitOnFlip:
            theseKeys = default_keyboard.getKeys(keyList=['t'], waitRelease=False)
            _default_keyboard_allKeys.extend(theseKeys)
            if len(_default_keyboard_allKeys):
                default_keyboard.keys = _default_keyboard_allKeys[-1].name  # just the last key pressed
                default_keyboard.rt = _default_keyboard_allKeys[-1].rt
                # a response ends the routine
                continueRoutine = False

        # check for quit (typically the Esc key)
        if end_exp_now or default_keyboard.getKeys(keyList=["escape"]):
            eye_tracker.sendMessage("key board escape")
            core.quit()

        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in endComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished

        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()

    # --- Ending Routine "end" ---
    for thisComponent in endComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # make sure the eyetracker recording stops
    if eye_tracker_stop.status != FINISHED:
        eye_tracker_stop.status = FINISHED
    eye_tracker.sendMessage("ET: eye-tracker stopped")
    # check responses
    if default_keyboard.keys in ['', [], None]:  # No response was made
        default_keyboard.keys = None
    this_exp.addData('default_keyboard.keys', default_keyboard.keys)
    if default_keyboard.keys != None:  # we had a response
        this_exp.addData('default_keyboard.rt', default_keyboard.rt)
    this_exp.nextEntry()
    # the Routine "end" was not non-slip safe, so reset the non-slip timer
    routine_timer.reset()


def flickering_checkerboard(eye_tracker, win, default_keyboard, routine_timer, this_exp) -> None:
    # current_fixation_color: str = choice(FIXATION_COLORS)
    current_fixation_color: int = randint(0, len(FIXATION_COLORS) - 1)
    next_color_change_time: float = uniform(*FIXATION_CHANGE_DELAY)

    checkerboard = visual.GratingStim(
        win=win,
        tex=None,  # np.array(((1, -1) * (NUM_SQUARES[0] // 2), (-1, 1) * (NUM_SQUARES[0] // 2)) * (NUM_SQUARES[0] // 2)),
        mask="none",
        units="pix",
        autoLog=False,
        interpolate=False,
        autoDraw=False,
    )

    precomputed_textures: list = [0] * num_textures
    for i in range(num_textures):
        if FLICKER_SHAPE == "square":
            num_square: int = NUM_SQUARES[i]
            precomputed_textures[i] = np.array(((1, -1) * (num_square // 2), (-1, 1) * (num_square // 2)) * (num_square // 2))

        elif FLICKER_SHAPE == "circle":
            num_rings: int = NUM_RINGS[i]
            num_sectors: int = NUM_SECTORS[i]
            outer_flicker_radius: int = OUTER_FLICKER_RADIUS[i]
            inner_flicker_radius: int = INNER_FLICKER_RADIUS[i]
            precomputed_textures[i] = make_ring_checkerboard(size=min(win.size),
                                                      outer_radius=outer_flicker_radius,
                                                      inner_radius=inner_flicker_radius,
                                                      n_radial=num_rings,
                                                      n_angular=num_sectors)

        elif FLICKER_SHAPE == "ring":
            num_rings: int = NUM_RINGS[i]
            num_sectors: int = NUM_SECTORS[i]
            outer_flicker_radius: int = OUTER_FLICKER_RADIUS[i]
            precomputed_textures[i] = make_ring_checkerboard(size=min(win.size),
                                                             outer_radius=outer_flicker_radius,
                                                             inner_radius=0,
                                                             n_radial=num_rings,
                                                             n_angular=num_sectors)


        elif FLICKER_SHAPE == "wedge":
            num_rings: int = NUM_RINGS[i]
            num_sectors: int = NUM_SECTORS[i]
            outer_flicker_radius: int = OUTER_FLICKER_RADIUS[i]
            inner_flicker_radius: int = INNER_FLICKER_RADIUS[i]
            start_angle: float = START_ANGLE[i]
            delta_angle: float = DELTA_ANGLE[i]
            precomputed_textures[i] = make_wedge_checkerboard(size=min(win.size),
                                                       outer_radius=outer_flicker_radius,
                                                       inner_radius=inner_flicker_radius,
                                                       n_radial=num_rings,
                                                       n_angular=num_sectors,
                                                       start_angle=start_angle,
                                                       end_angle=start_angle + delta_angle)

        elif FLICKER_SHAPE == "windmill":
            num_rings: int = NUM_RINGS[i]
            num_sectors: int = NUM_SECTORS[i]
            outer_flicker_radius: int = OUTER_FLICKER_RADIUS[i]
            inner_flicker_radius: int = INNER_FLICKER_RADIUS[i]
            num_wedges: int = NUM_WEDGES[i]
            precomputed_textures[i] = make_windmill_checkerboard(size=min(win.size),
                                                          n_radial=num_rings,
                                                          n_angular=num_sectors,
                                                          n_wedges=num_wedges,
                                                          outer_radius=outer_flicker_radius,
                                                          inner_radius=inner_flicker_radius)

    fixation_points = []
    if FIXATION_SHAPE == "cross":
        fixation_points.append(
            visual.ShapeStim(
                win=win,
                units="deg",
                vertices="cross",
                size=1.0/2,
                lineWidth=0.0,
                fillColor=FIXATION_COLORS[current_fixation_color],
            )
        )
    elif FIXATION_SHAPE == "circle":
        fixation_points.append(
            visual.Circle(
                win=win,
                units="deg",
                edges=64,
                radius=0.6/2,
                lineWidth=0.0,
                fillColor=FIXATION_COLORS[current_fixation_color],
                interpolate=True,
            )
        )
    elif FIXATION_SHAPE == "circle_cross_circle":
        fixation_points.extend([
            visual.Circle(
                win=win,
                units="deg",
                edges=64,
                radius=0.4/2,  # 0.01,
                lineWidth=0.0,
                fillColor=FIXATION_COLORS[current_fixation_color],
                interpolate=True,
            ),
            visual.ShapeStim(
                win=win,
                units="deg",
                vertices="cross",
                size=0.8/2,
                lineWidth=0.0,
                fillColor=FIXATION_COLORS[(current_fixation_color + 1) % len(FIXATION_COLORS)],
            ),
            visual.Circle(
                win=win,
                units="deg",
                edges=64,
                radius=0.1/2,  # 0.01,
                lineWidth=0.0,
                fillColor=FIXATION_COLORS[current_fixation_color],
                interpolate=True,
            )
        ])
    elif FIXATION_SHAPE == "peripheral_circles":
        outer_radius = 2 / 2  # 2° diameter -> 1° radius

        for i in range(8):
            angle = i * 360 / 8  # equally spaced around the circle
            x = outer_radius * np.cos(np.deg2rad(angle))
            y = outer_radius * np.sin(np.deg2rad(angle))

            circle = visual.Circle(
                win=win,
                units="deg",
                edges=64,
                radius=0.2/2,
                lineWidth=0.0,
                fillColor=FIXATION_COLORS[current_fixation_color],
                interpolate=True,
                pos=(x, y),  # place on the circle
            )
            fixation_points.append(circle)

    instruction = visual.TextStim(
        win=win, text=f"Waiting for scanner trigger ('{PROGRESS_SIGNAL}')...",
        height=0.08, color="white", units="norm"
    )

    instruction.draw()
    win.flip()

    while not default_keyboard.getKeys(keyList=[PROGRESS_SIGNAL], ignoreKeys=["escape"], waitRelease=False):
        core.wait(0.01)
        if default_keyboard.getKeys(keyList=["escape"]):
            eye_tracker.sendMessage("key board escape")
            core.quit()
    routine_timer.reset()  # FIXME: Certainly not the correct way to do this

    # --- Prepare to start Routine "checkerboard" ---
    continueRoutine = True
    # update component parameters for each repeat
    # keep track of which components have finished
    checkerboard_Components = [checkerboard, *fixation_points]

    for thisComponent in checkerboard_Components:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    t_last_flip: float = 0.0  # To time the flickering

    # --- Run Routine "checkerboard" ---
    routineForceEnded = not continueRoutine
    while continueRoutine and routine_timer.getTime() < routine_duration:
        # get current time
        t = routine_timer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routine_timer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame

        # *checkerboard* updates

        # if checkerboard is starting this frame...
        if checkerboard.status == NOT_STARTED and tThisFlip >= 0.0 - FRAME_TOLERANCE:
            # keep track of start time/frame for later
            checkerboard.frameNStart = frameN  # exact frame index
            checkerboard.tStart = t  # local t and not account for scr refresh
            checkerboard.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(checkerboard, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            this_exp.timestampOnFlip(win, 'checkerboard.started')
            # update status
            checkerboard.status = STARTED
            eye_tracker.sendMessage("ET: Start routine 'flickering_checkerboard'")
            t_last_flip = t

        # Get trial ID and determine flicker frequency and num squares from it
        trial_id: int = int(t / (ON_DURATION + OFF_DURATION))
        trial_id //= NUM_REPETITIONS
        frequency: float = FREQUENCIES[trial_id // num_textures]
        flicker_interval: float = 1.0 / (2.0 * frequency)  # Hz

        # Set correct texture:
        checkerboard.tex = precomputed_textures[trial_id % num_textures]
        # eye_tracker.sendMessage(f"ET: frequency = {frequency:.1f}, trial_id = {trial_id}")

        # Flicker logic
        if t - t_last_flip >= flicker_interval:
            checkerboard.contrast *= -1
            t_last_flip = t

        # Expansion
        if FLICKER_SHAPE == "ring":
            checkerboard.mask = annulus_mask(min(win.size),
                                             outer_radius=OUTER_FLICKER_RADIUS[trial_id % num_textures],
                                             ring_thickness=RING_THICKNESS[trial_id % num_textures],
                                             expansion_speed=EXPANSION_SPEED[trial_id % num_textures],
                                             t=t)

        # Rotation
        if FLICKER_SHAPE in ("wedge", "windmill"):
            checkerboard.ori = (ROTATION_SPEED[trial_id % num_textures] * t) % 360

        # Activate or deactivate checkerboard
        if (t % (ON_DURATION + OFF_DURATION)) < ON_DURATION:
            checkerboard.draw()

        # Fixation color change
        if t >= next_color_change_time:
            new_colors = [c for c in list(range(len(FIXATION_COLORS))) if c != current_fixation_color]
            current_fixation_color = choice(new_colors)

            if FIXATION_SHAPE in ("cross", "circle", "peripheral_circles"):
                for f_point in fixation_points:
                    f_point.fillColor = FIXATION_COLORS[current_fixation_color]

            elif FIXATION_SHAPE == "circle_cross_circle":
                fixation_points[0].fillColor = FIXATION_COLORS[current_fixation_color]
                fixation_points[1].fillColor = FIXATION_COLORS[(current_fixation_color + 1) % len(FIXATION_COLORS)]
                fixation_points[2].fillColor = FIXATION_COLORS[current_fixation_color]

            next_color_change_time += uniform(*FIXATION_CHANGE_DELAY)

        # Update screen at each frame
        for f_point in fixation_points:
            f_point.draw()

        # is it time to stop? (based on global clock, using actual start)
        if tThisFlipGlobal > checkerboard.tStartRefresh + routine_duration - FRAME_TOLERANCE:
            # keep track of stop time/frame for later
            checkerboard.tStop = t  # not accounting for scr refresh
            checkerboard.frameNStop = frameN  # exact frame index
            # add timestamp to datafile
            this_exp.timestampOnFlip(win, 'checkerboard.stopped')
            # update status
            checkerboard.status = FINISHED

        # check for quit (typically the Esc key)
        if end_exp_now or default_keyboard.getKeys(keyList=["escape"]):
            eye_tracker.sendMessage("key board escape")
            core.quit()

        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in checkerboard_Components:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished

        # refresh the screen
        # if continueRoutine:  # don't flip if this routine is over, or we'll get a blank screen
        win.flip()

    # --- Ending Routine "centered_dot" ---
    for thisComponent in checkerboard_Components:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
    if routineForceEnded:
        routine_timer.reset()
    else:
        routine_timer.addTime(-routine_duration)
    this_exp.nextEntry()


def save_data(this_exp: data.ExperimentHandler) -> None:
    filename = this_exp.dataFileName
    # these shouldn't be strictly necessary (should auto-save)
    this_exp.saveAsWideText(filename + '.csv', delim='auto')
    this_exp.saveAsPickle(filename)


def terminate(this_exp, win: visual.Window, eye_tracker = None) -> None:
    this_exp.abort()  # or data files will save again on exit
    
    # Flip one final time so any remaining win.callOnFlip() 
    # and win.timeOnFlip() tasks get executed before quitting
    win.flip()
    win.close()
    
    logging.flush()
    
    if eye_tracker:
        eye_tracker.setConnectionState(False)

    # terminate Python process
    core.quit()


def main():
    exp_name: str = "flickering_checkerboard"
    exp_info: dict = {
        'participant': f"{randint(0, 999999):06d}",
        'session': '001',
        'date|hid': data.getDateStr(),
        'expName|hid': exp_name,
        'psychopyVersion|hid': __version__,
    }

    # Startup functions
    exp_info: dict = show_exp_info_dlg(exp_info=exp_info)
    this_exp: data.ExperimentHandler = setup_data(exp_info=exp_info)
    _: logging.LogFile = setup_logging(filename=this_exp.dataFileName)
    win: visual.Window = setup_window(exp_info=exp_info, full_scr=FULL_SCREEN)
    eye_tracker, default_keyboard = setup_devices(exp_info=exp_info, win=win)

    # Clocks to sync with MRI scanner
    routine_timer: core.Clock = core.Clock()  # to track time remaining of each (possibly non-slip) routine

    # Routines (calibration, start ET, ...., stop ET)
    run_calibration(default_keyboard, routine_timer, win, eye_tracker)
    start_eye_tracker(eye_tracker, win, routine_timer, this_exp, default_keyboard)
    flickering_checkerboard(eye_tracker, win, default_keyboard, routine_timer, this_exp)
    stop_eye_tracker(eye_tracker, win, default_keyboard, routine_timer, this_exp)

    # Shutdown function
    save_data(this_exp=this_exp)
    terminate(this_exp=this_exp, win=win, eye_tracker=eye_tracker)


if __name__ == "__main__":
    print(f"Psychopy version: {__version__}")
    main()