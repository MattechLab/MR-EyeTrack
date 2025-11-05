#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A copy of `fixed_dot-16_grid_T1w.py` rewritten to avoid `iohub` and hardware eyetracker APIs.
It uses PsychoPy visual, core and event modules only, and a safe socket sender for external messages.

This script aims to preserve the experiment flow: intro, calibration pause (simulated), centered dots, grid dots,
and final message saving data via ExperimentHandler. It intentionally does not use `iohub`,
`hardware.eyetracker` or any tracker-specific calls.

Run: python fixed_dot-16_grid_T1w_noiohub.py
"""

from psychopy import prefs
prefs.hardware['audioLib'] = 'ptb'
prefs.hardware['audioLatencyMode'] = '3'
from psychopy import sound, gui, visual, core, data, event, logging
from psychopy.constants import (NOT_STARTED, STARTED, FINISHED)
import socket
import os
import sys
import math
from numpy.random import randint, shuffle

# Safe send_message that won't crash if no server is listening
def send_message(message, addr="localhost", port=2023):
    try:
        client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        client_socket.settimeout(0.5)
        client_socket.connect((addr, port))
        # ensure bytes
        if isinstance(message, str):
            message = message.encode('utf-8')
        client_socket.sendall(message)
        client_socket.close()
    except Exception:
        # silently ignore networking errors so the experiment can run without a tracker
        pass

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)

psychopyVersion = '2024.2.1'
expName = 'fixed_dot-16_grid_T1w_noiohub'
expInfo = {
    'participant': f"{randint(0, 999999):06.0f}",
    'session': '001',
}
# dialog
try:
    dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)
    if not dlg.OK:
        core.quit()
except Exception:
    # Running in non-GUI environment - proceed with defaults
    pass

expInfo['date'] = data.getDateStr()
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])
thisExp = data.ExperimentHandler(name=expName, extraInfo=expInfo, dataFileName=filename)
logFile = logging.LogFile(filename + '.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)

# Window
win = visual.Window(size=[800, 600], fullscr=True, screen=0, winType='pyglet',
                    monitor='testMonitor', color=[0, 0, 0], units='norm')
win.mouseVisible = False

frameRate = win.getActualFrameRate()
if frameRate is not None:
    frameDur = 1.0 / round(frameRate)
else:
    frameDur = 1.0 / 60.0

# Components
waiting_trigger = visual.TextStim(win=win, name='waiting_trigger',
                                  text="The program is ready for the scanner trigger. Press 's' to proceed manually.",
                                  pos=(0, -0.4), height=0.12, wrapWidth=1.7, color='white')
fix_desc = visual.TextStim(win=win, name='fix_desc',
                           text='In this task you will see a dot moving randomly within different positions on the screen. You have to follow the dot :)',
                           pos=(0, 0.25), height=0.12, wrapWidth=1.0, color='white')

# Dot stimuli: composite dot made of components (outer circle, cross as two rectangles, inner circle)
dot = []

# Geometry (units = 'norm') -- use original sizes but compensate for aspect when drawing
outer_radius = 0.1 / 2       # original outer radius in norm units
inner_radius = 0.02 / 2       # original inner radius in norm units
cross_length = outer_radius * 2.1
cross_thickness = inner_radius * 2.1

# Aspect and x-axis scale for norm units: x should be scaled by 1/aspect to counter horizontal stretching
aspect = float(win.size[0]) / float(win.size[1]) if hasattr(win, 'size') else 1.0
x_scale = 1.0 / aspect

# helper: build polygon vertices approximating a circle, scaled in x by x_scale
def _circle_vertices(radius, n=64, x_scale=1.0):
    verts = []
    for i in range(n):
        theta = (2.0 * math.pi * i) / n
        x = math.cos(theta) * radius * x_scale
        y = math.sin(theta) * radius
        verts.append((x, y))
    return verts

# outer circle (as ShapeStim polygon so we can scale x separately)
outer_verts = _circle_vertices(outer_radius, n=64, x_scale=x_scale)
# inner circle verts
inner_verts = _circle_vertices(inner_radius, n=64, x_scale=x_scale)

dot.extend([
    visual.ShapeStim(
        win=win,
        units='norm',
        vertices=outer_verts,
        fillColor=(-1, -1, -1),
        lineColor=(-1, -1, -1),
        interpolate=True,
    ),
    # cross as rects but scale width by x_scale so they match norm coordinate system
    visual.Rect(
        win=win,
        units='norm',
        width=cross_length * x_scale,
        height=cross_thickness,
        fillColor=(0, 0, 0),
        lineColor=(0, 0, 0),
    ),
    visual.Rect(
        win=win,
        units='norm',
        width=cross_thickness * x_scale,
        height=cross_length,
        fillColor=(0, 0, 0),
        lineColor=(0, 0, 0),
    ),
    visual.ShapeStim(
        win=win,
        units='norm',
        vertices=inner_verts,
        fillColor=(-1, -1, -1),
        lineColor=(-1, -1, -1),
        interpolate=True,
    ),
])

# Helper utilities to manage the composite dot
def set_dot_autodraw(on=True):
    for comp in dot:
        comp.setAutoDraw(on)

def set_dot_pos(pos):
    for comp in dot:
        comp.pos = pos

# Experiment timing / positions
grid_size = 4
dot_size = 0.05
# Preserve same t_dot as original
t_dot = 5 * 6.2 / 8.01

# Compute reasonable normalized offsets (keep in norm units)
x_offset_norm = 1.33 * 2 / 3
y_offset_norm = 1 / 2
positions = {
    (0, y_offset_norm): "up",
    (0, -y_offset_norm): "down",
    (-x_offset_norm, 0): "left",
    (x_offset_norm, 0): "right"
}

# Ensure data directory exists
os.makedirs(os.path.join(_thisDir, 'data'), exist_ok=True)

# Run experiment inside a try/finally so we always save and close
try:
    # Start of experiment
    send_message("ET: Start experiment 'dots'")

    # Intro / waiting for trigger (use event for keypress)
    event.clearEvents()
    waiting_trigger.setAutoDraw(True)
    fix_desc.setAutoDraw(True)
    win.flip()
    send_message("ET: waiting for scanner trigger or manual 's' press")
    keys = event.waitKeys(keyList=['s', 'escape', 'q'])
    if 'escape' in keys or 'q' in keys:
        raise KeyboardInterrupt
    waiting_trigger.setAutoDraw(False)
    fix_desc.setAutoDraw(False)
    thisExp.addData('waiting_trigger', 'triggered')

    # Simulated calibration (we'll just wait and show a message)
    cal_msg = visual.TextStim(win=win, text='Calibration (simulated). Please fixate on the dot when it appears.', color='white')
    cal_msg.draw()
    win.flip()
    send_message("ET: calibration run (simulated)")
    core.wait(1.5)

    # Repeat centered dot 6 times (use composite dot at center)
    for i in range(6):
        set_dot_pos((0, 0))
        set_dot_autodraw(True)
        send_message("ET: Start routine 'centered_dot'")
        thisExp.addData(f'centered_dot_{i}_start', core.getTime())
        win.flip()
        core.wait(5.0)
        thisExp.addData(f'centered_dot_{i}_end', core.getTime())
        set_dot_autodraw(False)

    # Dots routine repeated 30 times as in original file
    for rep in range(30):
        # Shuffle positions
        shuffled_positions = list(positions.items())
        shuffle(shuffled_positions)
        current_position_index = 0
        total_positions = len(shuffled_positions)
        pos, direction = shuffled_positions[current_position_index]
        set_dot_pos(pos)
        time_of_last_change = core.getTime()
        send_message("ET: Start routine 'dots'")
        set_dot_autodraw(True)
        win.flip()
        continueRoutine = True
        while continueRoutine:
            t = core.getTime()
            # change position if enough time passed
            if t - time_of_last_change >= t_dot:
                current_position_index += 1
                if current_position_index < total_positions:
                    pos, direction = shuffled_positions[current_position_index]
                    set_dot_pos(pos)
                    send_message(f"ET: dot moved {direction}!")
                    # log the movement
                    thisExp.addData('rep', rep)
                    thisExp.addData('position_label', direction)
                    thisExp.addData('position', pos)
                    thisExp.addData('time', t)
                    thisExp.nextEntry()
                    time_of_last_change = t
                else:
                    continueRoutine = False
            # check for escape or quit key
            if event.getKeys(keyList=['escape', 'q']):
                raise KeyboardInterrupt
            win.flip()
        set_dot_autodraw(False)

    # Final centered dots 5 times (use composite dot at center)
    for i in range(5):
        set_dot_pos((0, 0))
        set_dot_autodraw(True)
        send_message("ET: Start routine 'centered_dot'")
        win.flip()
        core.wait(5.0)
        set_dot_autodraw(False)

    # End routine
    text = visual.TextStim(win=win, text="End press 't'", color='white')
    text.setAutoDraw(True)
    send_message("ET: Prepare to start routine 'end'")
    win.flip()
    keys = event.waitKeys(keyList=['t', 'escape', 'q'])
    if 'escape' in keys or 'q' in keys:
        raise KeyboardInterrupt
    text.setAutoDraw(False)
    thisExp.addData('end_key', keys)

except KeyboardInterrupt:
    # Graceful early exit triggered by escape
    send_message("key board escape")
    thisExp.addData('early_exit', True)
except Exception:
    logging.error("Unexpected error during experiment", exc_info=True)
    raise
finally:
    # Save and close (always run)
    try:
        thisExp.saveAsWideText(filename + '.csv', delim='auto')
        thisExp.saveAsPickle(filename)
    except Exception:
        logging.error('Failed saving experiment data', exc_info=True)
    logging.flush()
    win.flip()
    try:
        win.close()
    except Exception:
        pass
    core.quit()
