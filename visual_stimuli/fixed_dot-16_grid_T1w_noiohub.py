#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Variant of fixed_dot-16_grid_T1w.py with all iohub-related code removed.
Automatically generated to run without psychopy.iohub and Eyelink hardware calls.
"""

from psychopy import prefs
from psychopy import plugins
plugins.activatePlugins()
prefs.hardware['audioLib'] = 'ptb'
prefs.hardware['audioLatencyMode'] = '3'
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors, layout, hardware
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED, STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)
import socket
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average, sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle, choice as randchoice
import os  # handy system and path functions
import sys  # to get file system encoding
from psychopy.monitors import Monitor

from psychopy.hardware import keyboard

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)
# Store info about the experiment session
psychopyVersion = '2024.2.1'
expName = 'fixed_dot-16_grid_T1w_noiohub'  # from the Builder filename that created this script
expInfo = {
    'participant': f"{randint(0, 999):03.0f}",
    'session': '001',
}
# --- Show participant info dialog --
dlg = gui.DlgFromDict(dictionary=expInfo, sortKeys=False, title=expName)
if dlg.OK == False:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath=None,
    savePickle=True, saveWideText=True,
    dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp
frameTolerance = 0.001  # how close to onset before 'same' frame

# Start Code - component code to be run after the window creation

# ET?
DUMMY_MODE = True  # Set to False when running with the real eyetracker (not used in this variant)

# Monitor setup
WIN_SIZE = (800, 600)
monitor = Monitor("expMonitor")
monitor.setWidth(369.54e-3)  # screen width in meters
monitor.setDistance(1020e-3)  # viewing distance in meters
monitor.setSizePix(WIN_SIZE)  # screen resolution
monitor.saveMon()

# --- Setup the Window ---
win = visual.Window(
    size=WIN_SIZE,
    fullscr=True,
    screen=0, 
    winType='pyglet',
    allowStencil=False,
    monitor='expMonitor',
    color=(0, 0, 0),
    colorSpace='rgb',
    backgroundImage='',
    backgroundFit='none',
    blendMode='avg',
    useFBO=True, 
    units='deg')
win.mouseVisible = False
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess

# Compute a dynamic text height based on the monitor/window size when units are 'deg'.
try:
    screen_width_deg = monitor.getWidth() / monitor.getDistance() * (180.0 / np.pi)
    screen_height_deg = screen_width_deg * (win.size[1] / win.size[0])
    _text_height = max(0.02, min(0.18, screen_height_deg * 0.06))
except Exception:
    _text_height = 0.12

# Setup eyetracking placeholders removed in this variant

# Setup keyboard (no iohub)
defaultKeyboard = keyboard.Keyboard()

# --- Initialize components for Routine "trail" ---
waiting_trigger = visual.TextStim(win=win, name='waiting_trigger',
    text="The program is ready for the scanner trigger. Press 's' to proceed manually.",
    font='Open Sans',
    pos=(0, -0.4), height=_text_height, wrapWidth=1.7, ori=0.0, 
    color='white', colorSpace='rgb', opacity=1.0, 
    languageStyle='LTR',
    depth=0.0);
key_resp = keyboard.Keyboard()
fix_desc = visual.TextStim(win=win, name='fix_desc',
    text='In this task you will see a dot moving randomly within different positions on the screen. You have to follow the dot :)',
    font='Open Sans',
    pos=(0, -0.4), height=_text_height, wrapWidth=1.7, ori=0.0, 
    color='white', colorSpace='rgb', opacity=1.0, 
    languageStyle='LTR',
    depth=0.0);

# --- Initialize components for Routine "start_ET" ---
# ET-control removed; we'll keep a simple placeholder function if needed

dot = []

dot.append([
    visual.Circle(
        win=win,
        units="deg",
        edges=64,
        radius=0.4/2,
        lineWidth=0.0,
        fillColor=[0.5, 0.5, 0.5],
        interpolate=True,
    ),
    visual.ShapeStim(
        win=win,
        units="deg",
        vertices="cross",
        size=0.8/2,
        lineWidth=0.0,
        fillColor=[0.5, 0.5, 0.5],
    ),
    visual.Circle(
        win=win,
        units="deg",
        edges=64,
        radius=0.1/2,  # 0.01,
        lineWidth=0.0,
        fillColor=[0.5, 0.5, 0.5],
        interpolate=True,
    )
])

# Keep the raw configuration (list of stimuli) available and expose a small
# wrapper so the rest of the code can treat `dot` like a single stimulus
# while still allowing multiple stimulus variants.
dot_config = dot[0]

class DotGroup:
    """A thin wrapper around a list of visual stimuli.

    - Delegates attribute access to the currently selected stimulus.
    - Assigning to the attribute `pos` updates all configured stimuli so they
      stay aligned; other attribute sets are forwarded to the active stimulus.
    - Use `select(index)` to choose which stimulus is active.
    """
    def __init__(self, items, default_index=None):
        self._items = list(items)
        if default_index is None:
            default_index = len(self._items) - 1
        self._idx = int(default_index)

    def select(self, idx: int):
        self._idx = int(idx)

    def current(self):
        return self._items[self._idx]

    def __getattr__(self, name):
        # Delegate missing attributes to the currently selected stimulus
        return getattr(self.current(), name)

    def __setattr__(self, name, value):
        if name in ('_items', '_idx'):
            object.__setattr__(self, name, value)
            return
        if name == 'pos':
            # Keep all configured stimuli at the same position
            for it in self._items:
                setattr(it, 'pos', value)
            return
        # Forward other sets to the active stimulus
        setattr(self.current(), name, value)

    def __repr__(self):
        return f"DotGroup(active={self._idx}, items={self._items})"

# Replace `dot` with the wrapper while keeping `dot_config` intact
dot = DotGroup(dot_config)

# Begin Experiment
t_dot = 5*6.2/8.01  # Seconds of showing the dot per position. TR is decreased from 8.1 to 6.2ms

# Define offsets (in degrees of visual angle)
x_offset_deg = 5
y_offset_deg = 5

# Define positions in 'deg' units
positions = {
    (0, y_offset_deg): "up",
    (0, -y_offset_deg): "down",
    (-x_offset_deg, 0): "left",
    (x_offset_deg, 0): "right"
}

# --- Initialize components for Routine "end" ---
text = visual.TextStim(win=win, name='text',
    text="End press 't'",
    font='Open Sans',
    pos=(0, 0), height=_text_height, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
key_resp_3 = keyboard.Keyboard()

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.Clock()  # to track time remaining of each (possibly non-slip) routine 

# define target for calibration_2 (removed eyetracker dependency)
calibration_2Target = visual.TargetStim(win, 
    name='calibration_2Target',
    radius=0.15, fillColor=[0.5, 0.5, 0.5], borderColor=[0.5, 0.5, 0.5], lineWidth=2.0,
    innerRadius=0.07, innerFillColor=[0.5, 0.5, 0.5], innerBorderColor=[0.5, 0.5, 0.5], innerLineWidth=2.0,
    colorSpace='rgb', units=None
)

# reset timers
continueRoutine = True
# update component parameters for each repeat
key_resp.keys = []
key_resp.rt = []
_key_resp_allKeys = []
# keep track of which components have finished
trailComponents = [waiting_trigger, key_resp, fix_desc]
for thisComponent in trailComponents:
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

# --- Run Routine "trail" ---
routineForceEnded = not continueRoutine
while continueRoutine:
    # get current time
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *waiting_trigger* updates
    if waiting_trigger.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        waiting_trigger.frameNStart = frameN
        waiting_trigger.tStart = t
        waiting_trigger.tStartRefresh = tThisFlipGlobal
        win.timeOnFlip(waiting_trigger, 'tStartRefresh')
        waiting_trigger.status = STARTED
        waiting_trigger.setAutoDraw(True)

    if waiting_trigger.status == STARTED:
        pass

    # *key_resp* updates
    waitOnFlip = False
    if key_resp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        key_resp.frameNStart = frameN
        key_resp.tStart = t
        key_resp.tStartRefresh = tThisFlipGlobal
        win.timeOnFlip(key_resp, 'tStartRefresh')
        key_resp.status = STARTED
        waitOnFlip = True
        win.callOnFlip(key_resp.clock.reset)
        win.callOnFlip(key_resp.clearEvents, eventType='keyboard')
    if key_resp.status == STARTED and not waitOnFlip:
        theseKeys = key_resp.getKeys(keyList=['s'], waitRelease=False)
        _key_resp_allKeys.extend(theseKeys)
        if len(_key_resp_allKeys):
            key_resp.keys = _key_resp_allKeys[-1].name
            key_resp.rt = _key_resp_allKeys[-1].rt
            continueRoutine = False

    # *fix_desc* updates
    if fix_desc.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        fix_desc.frameNStart = frameN
        fix_desc.tStart = t
        fix_desc.tStartRefresh = tThisFlipGlobal
        win.timeOnFlip(fix_desc, 'tStartRefresh')
        fix_desc.status = STARTED
        fix_desc.setAutoDraw(True)

    if fix_desc.status == STARTED:
        pass

    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()

    # check if all components have finished
    if not continueRoutine:
        routineForceEnded = True
        break
    continueRoutine = False
    for thisComponent in trailComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break

    if continueRoutine:
        win.flip()

# --- Ending Routine "trail" ---
for thisComponent in trailComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# check responses
if key_resp.keys in ['', [], None]:
    key_resp.keys = None
thisExp.addData('key_resp.keys',key_resp.keys)
if key_resp.keys != None:
    thisExp.addData('key_resp.rt', key_resp.rt)
thisExp.nextEntry()
routineTimer.reset()

# --- Prepare to start Routine "start_ET" ---
continueRoutine = True
start_ETComponents = []  # no eyetracker components
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

# run a short pause instead of eyetracker calibration
core.wait(0.2)

# Repeat the centered_dot routine 6 times
for _ in range(6):
    continueRoutine = True
    dotComponents = [dot]
    for thisComponent in dotComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1

    routineForceEnded = not continueRoutine
    while continueRoutine and routineTimer.getTime() < 5.0:
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1
        if dot.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            dot.frameNStart = frameN
            dot.tStart = t
            dot.tStartRefresh = tThisFlipGlobal
            win.timeOnFlip(dot, 'tStartRefresh')
            thisExp.timestampOnFlip(win, 'dot.started')
            dot.status = STARTED
            dot.setAutoDraw(True)
        if dot.status == STARTED:
            pass
        if dot.status == STARTED:
            if tThisFlipGlobal > dot.tStartRefresh + 5.0-frameTolerance:
                dot.tStop = t
                dot.frameNStop = frameN
                thisExp.timestampOnFlip(win, 'dot.stopped')
                dot.status = FINISHED
                dot.setAutoDraw(False)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        if not continueRoutine:
            routineForceEnded = True
            break
        continueRoutine = False
        for thisComponent in dotComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break
        if continueRoutine:
            win.flip()
    for thisComponent in dotComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    if routineForceEnded:
        routineTimer.reset()
    else:
        routineTimer.addTime(-5.000000)
    thisExp.nextEntry()

# Repeat the dots routine 30 times
for _ in range(30):
    continueRoutine = True
    shuffled_positions = list(positions.items())
    shuffle(shuffled_positions)
    current_position_index = 0
    total_positions = len(shuffled_positions)
    dot.pos, direction = shuffled_positions[current_position_index]
    time_of_last_change = 0
    continueRoutine = True
    # update component parameters for each repeat
    dotsComponents = [dot]
    for thisComponent in dotsComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1

    routineForceEnded = not continueRoutine
    while continueRoutine:
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1
        if dot.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            dot.frameNStart = frameN
            dot.tStart = t
            dot.tStartRefresh = tThisFlipGlobal
            win.timeOnFlip(dot, 'tStartRefresh')
            thisExp.timestampOnFlip(win, 'dot.started')
            dot.status = STARTED
            dot.setAutoDraw(True)
        if dot.status == STARTED:
            pass

        t = routineTimer.getTime()
        if t - time_of_last_change >= t_dot:
            current_position_index += 1
            if current_position_index < total_positions:
                dot.pos, direction = shuffled_positions[current_position_index]
                # optionally log to experiment file instead of eyetracker
                thisExp.addData('dot_moved', direction)
                time_of_last_change = t
            else:
                continueRoutine = False

        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()

        if not continueRoutine:
            routineForceEnded = True
            break
        continueRoutine = False
        for thisComponent in dotsComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break
        if continueRoutine:
            win.flip()

    for thisComponent in dotsComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    routineTimer.reset()
    thisExp.nextEntry()

# Repeat the centered_dot routine 5 times
for _ in range(5):
    continueRoutine = True
    dotComponents = [dot]
    for thisComponent in dotComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1

    routineForceEnded = not continueRoutine
    while continueRoutine and routineTimer.getTime() < 5.0:
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1
        if dot.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            dot.frameNStart = frameN
            dot.tStart = t
            dot.tStartRefresh = tThisFlipGlobal
            win.timeOnFlip(dot, 'tStartRefresh')
            thisExp.timestampOnFlip(win, 'dot.started')
            dot.status = STARTED
            dot.setAutoDraw(True)
        if dot.status == STARTED:
            pass
        if dot.status == STARTED:
            if tThisFlipGlobal > dot.tStartRefresh + 5.0-frameTolerance:
                dot.tStop = t
                dot.frameNStop = frameN
                thisExp.timestampOnFlip(win, 'dot.stopped')
                dot.status = FINISHED
                dot.setAutoDraw(False)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        if not continueRoutine:
            routineForceEnded = True
            break
        continueRoutine = False
        for thisComponent in dotComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break
        if continueRoutine:
            win.flip()
    for thisComponent in dotComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    if routineForceEnded:
        routineTimer.reset()
    else:
        routineTimer.addTime(-5.000000)
    thisExp.nextEntry()

# --- Prepare to start Routine "end" ---
continueRoutine = True
key_resp_3.keys = []
key_resp_3.rt = []
_key_resp_3_allKeys = []
# optionally record that experiment is finishing
thisExp.addData('experiment_end', True)
endComponents = [text, key_resp_3]
for thisComponent in endComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
frameN = -1

routineForceEnded = not continueRoutine
while continueRoutine:
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1
    if text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        text.frameNStart = frameN
        text.tStart = t
        text.tStartRefresh = tThisFlipGlobal
        win.timeOnFlip(text, 'tStartRefresh')
        thisExp.timestampOnFlip(win, 'text.started')
        text.status = STARTED
        text.setAutoDraw(True)
    if text.status == STARTED:
        pass

    # *key_resp_3* updates
    waitOnFlip = False
    if key_resp_3.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        key_resp_3.frameNStart = frameN
        key_resp_3.tStart = t
        key_resp_3.tStartRefresh = tThisFlipGlobal
        win.timeOnFlip(key_resp_3, 'tStartRefresh')
        thisExp.timestampOnFlip(win, 'key_resp_3.started')
        key_resp_3.status = STARTED
        waitOnFlip = True
        win.callOnFlip(key_resp_3.clock.reset)
        win.callOnFlip(key_resp_3.clearEvents, eventType='keyboard')
    if key_resp_3.status == STARTED and not waitOnFlip:
        theseKeys = key_resp_3.getKeys(keyList=['t'], waitRelease=False)
        _key_resp_3_allKeys.extend(theseKeys)
        if len(_key_resp_3_allKeys):
            key_resp_3.keys = _key_resp_3_allKeys[-1].name
            key_resp_3.rt = _key_resp_3_allKeys[-1].rt
            continueRoutine = False

    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()

    if not continueRoutine:
        routineForceEnded = True
        break
    continueRoutine = False
    for thisComponent in endComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break
    if continueRoutine:
        win.flip()

for thisComponent in endComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)

# final save and cleanup
thisExp.saveAsWideText(filename + '.csv', delim='auto')
thisExp.saveAsPickle(filename)
logging.flush()
thisExp.abort()
win.close()
core.quit()
