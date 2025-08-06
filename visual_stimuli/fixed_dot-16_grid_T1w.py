#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2024.2.1),
    on January 27, 2025, at 18:21
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

# --- Import packages ---
from psychopy import prefs
from psychopy import plugins
plugins.activatePlugins()
prefs.hardware['audioLib'] = 'ptb'
prefs.hardware['audioLatencyMode'] = '3'
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors, layout, iohub, hardware
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED, STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)
import socket
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average, sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle, choice as randchoice
import os  # handy system and path functions
import sys  # to get file system encoding

import psychopy.iohub as io
from psychopy.hardware import keyboard

def send_message(message, addr="localhost", port=2023):
    client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    client_socket.connect((addr, port))
    client_socket.sendall(message)
    client_socket.close()

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)
# Store info about the experiment session
psychopyVersion = '2024.2.1'
expName = 'fixed_dot-16_grid_T1w'  # from the Builder filename that created this script
expInfo = {
    'participant': f"{randint(0, 999999):06.0f}",
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
    originPath='/c:/Users/jaime.barranco/Desktop/repos/mattechlab/MR-EyeTrack/visual_stimuli/fixed_dot-16_grid_T1w.py',
    savePickle=True, saveWideText=True,
    dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp
frameTolerance = 0.001  # how close to onset before 'same' frame

# Start Code - component code to be run after the window creation

# --- Setup the Window ---
win = visual.Window(
    size=[800, 600], fullscr=True, screen=0, 
    winType='pyglet', allowStencil=False,
    monitor='testMonitor', color=[-1,-1,-1], colorSpace='rgb',
    backgroundImage='', backgroundFit='none',
    blendMode='avg', useFBO=True, 
    units='norm')
win.mouseVisible = False
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess
# --- Setup input devices ---
ioConfig = {}

# Setup eyetracking
ioConfig['eyetracker.hw.sr_research.eyelink.EyeTracker'] = {
    'name': 'tracker',
    'model_name': 'EYELINK 1000 DESKTOP',
    'simulation_mode': False,
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
}

# Setup iohub keyboard
ioConfig['Keyboard'] = dict(use_keymap='psychopy')

ioSession = '1'
if 'session' in expInfo:
    ioSession = str(expInfo['session'])
ioServer = io.launchHubServer(window=win, **ioConfig)
eyetracker = ioServer.getDevice('tracker')

# create a default keyboard (e.g. to check for escape)
defaultKeyboard = keyboard.Keyboard(backend='iohub')

# --- Initialize components for Routine "trail" ---
waiting_trigger = visual.TextStim(win=win, name='waiting_trigger',
    text="The program is ready for the scanner trigger. Press 's' to proceed manually.",
    font='Open Sans',
    pos=(0, -0.4), height=0.12, wrapWidth=1.7, ori=0.0, 
    color='white', colorSpace='rgb', opacity=1.0, 
    languageStyle='LTR',
    depth=0.0);
key_resp = keyboard.Keyboard()
fix_desc = visual.TextStim(win=win, name='fix_desc',
    text='In this task you will see a dot moving randomly within different positions on the screen. You have to follow the dot :)',
    font='Open Sans',
    pos=(0, 0.25), height=0.12, wrapWidth=1.0, ori=0.0, 
    color='white', colorSpace='rgb', opacity=1.0, 
    languageStyle='LTR',
    depth=-2.0);

# --- Initialize components for Routine "start_ET" ---
etRecord = hardware.eyetracker.EyetrackerControl(
    tracker=eyetracker,
    actionType='Start Only'
)
# Send logs to the ET
def send_message(message, addr="localhost", port=2023):
    client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    client_socket.connect((addr, port))
    client_socket.sendall(message)
    client_socket.close()

# --- Initialize components for Routine "centered_dot" ---
dot_centered = visual.ShapeStim(
    win=win, name='dot_centered',
    size=(0.1, 0.1), vertices='circle',
    ori=0.0, pos=(0, 0), anchor='center',
    lineWidth=1.0, colorSpace='rgb', lineColor=[0.5000, 0.5000, 0.5000], fillColor=[0.5000, 0.5000, 0.5000],
    opacity=1.0, depth=0.0, interpolate=True)

# --- Initialize components for Routine "dots" ---
dot = visual.ShapeStim(
    win=win, name='dot',
    size=(0.1, 0.1), vertices='circle',
    ori=0.0, pos=(0, 0), anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=[0.5, 0.5, 0.5], fillColor=[0.5, 0.5, 0.5],
    opacity=1.0, depth=0.0, interpolate=True)
# Begin Experiment
grid_size = 4  # 4x4 grid
dot_size = 0.05  # Size of the grey dot
t_dot = 3  # Seconds of showing the dot per position

# Get the screen dimensions
# In norm units, screen goes from -1 to +1 vertically, and aspect-ratio-scaled horizontally.
# So on an 800×600 display, full x-range is from -800/600 = -1.333 to +1.333
# We normalize the pixel offsets to norm coordinates.
screen_width, screen_height = win.size
x_offset_norm = 1.33*2/3  # of 1.33
y_offset_norm = 1/2  # of 1.00

# Define the positions in 'norm' units
positions = {
    (0, y_offset_norm): "up",
    (0, -y_offset_norm): "down",
    (-x_offset_norm, 0): "left",
    (x_offset_norm, 0): "right"
}

ioServer.getDevice('tracker').sendMessage("ET: Start experiment 'dots'")

# --- Initialize components for Routine "end" ---
text = visual.TextStim(win=win, name='text',
    text="End press 't'",
    font='Open Sans',
    pos=(0, 0), height=0.12, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
ET_stop = hardware.eyetracker.EyetrackerControl(
    tracker=eyetracker,
    actionType='Stop Only'
)
key_resp_3 = keyboard.Keyboard()

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.Clock()  # to track time remaining of each (possibly non-slip) routine 
# define target for calibration_2
calibration_2Target = visual.TargetStim(win, 
    name='calibration_2Target',
    radius=0.15, fillColor=[0.5, 0.5, 0.5], borderColor=[0.5, 0.5, 0.5], lineWidth=2.0,
    innerRadius=0.07, innerFillColor=[0.5, 0.5, 0.5], innerBorderColor=[0.5, 0.5, 0.5], innerLineWidth=2.0,
    colorSpace='rgb', units=None
)
# define parameters for calibration_2
calibration_2 = hardware.eyetracker.EyetrackerCalibration(win, 
    eyetracker, calibration_2Target,
    units=None, colorSpace='rgb',
    progressMode='time', targetDur=1.5, expandScale=1.5,
    targetLayout='FIVE_POINTS', randomisePos=True, textColor='white',
    movementAnimation=True, targetDelay=1.0
)
# run calibration
calibration_2.run()
# clear any keypresses from during calibration_2 so they don't interfere with the experiment
defaultKeyboard.clearEvents()
# the Routine "calibration_2" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# --- Prepare to start Routine "trail" ---
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
    
    # if waiting_trigger is starting this frame...
    if waiting_trigger.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        waiting_trigger.frameNStart = frameN  # exact frame index
        waiting_trigger.tStart = t  # local t and not account for scr refresh
        waiting_trigger.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(waiting_trigger, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'waiting_trigger.started')
        # update status
        waiting_trigger.status = STARTED
        waiting_trigger.setAutoDraw(True)
    
    # if waiting_trigger is active this frame...
    if waiting_trigger.status == STARTED:
        # update params
        pass
    
    # *key_resp* updates
    waitOnFlip = False
    
    # if key_resp is starting this frame...
    if key_resp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        key_resp.frameNStart = frameN  # exact frame index
        key_resp.tStart = t  # local t and not account for scr refresh
        key_resp.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(key_resp, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'key_resp.started')
        # update status
        key_resp.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(key_resp.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(key_resp.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if key_resp.status == STARTED and not waitOnFlip:
        theseKeys = key_resp.getKeys(keyList=['s'], waitRelease=False)
        _key_resp_allKeys.extend(theseKeys)
        if len(_key_resp_allKeys):
            key_resp.keys = _key_resp_allKeys[-1].name  # just the last key pressed
            key_resp.rt = _key_resp_allKeys[-1].rt
            # a response ends the routine
            continueRoutine = False
    
    # *fix_desc* updates
    
    # if fix_desc is starting this frame...
    if fix_desc.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        fix_desc.frameNStart = frameN  # exact frame index
        fix_desc.tStart = t  # local t and not account for scr refresh
        fix_desc.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(fix_desc, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'fix_desc.started')
        # update status
        fix_desc.status = STARTED
        fix_desc.setAutoDraw(True)
    
    # if fix_desc is active this frame...
    if fix_desc.status == STARTED:
        # update params
        pass
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        ioServer.getDevice('tracker').sendMessage("key board escape")
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        routineForceEnded = True
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in trailComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# --- Ending Routine "trail" ---
for thisComponent in trailComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# check responses
if key_resp.keys in ['', [], None]:  # No response was made
    key_resp.keys = None
thisExp.addData('key_resp.keys',key_resp.keys)
if key_resp.keys != None:  # we had a response
    thisExp.addData('key_resp.rt', key_resp.rt)
thisExp.nextEntry()
# the Routine "trail" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# --- Prepare to start Routine "start_ET" ---
continueRoutine = True
# update component parameters for each repeat
# keep track of which components have finished
start_ETComponents = [etRecord]
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
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    # *etRecord* updates
    
    # if etRecord is starting this frame...
    if etRecord.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        etRecord.frameNStart = frameN  # exact frame index
        etRecord.tStart = t  # local t and not account for scr refresh
        etRecord.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(etRecord, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'etRecord.started')
        # update status
        etRecord.status = STARTED
        # Run 'Begin Routine' code from code_channel2
        ioServer.getDevice('tracker').sendMessage("Hello tracker record")
    
    
    # if etRecord is stopping this frame...
    if etRecord.status == STARTED:
        # is it time to stop? (based on global clock, using actual start)
        if tThisFlipGlobal > etRecord.tStartRefresh + 0-frameTolerance:
            # keep track of stop time/frame for later
            etRecord.tStop = t  # not accounting for scr refresh
            etRecord.frameNStop = frameN  # exact frame index
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'etRecord.stopped')
            # update status
            etRecord.status = FINISHED
            ioServer.getDevice('tracker').sendMessage("Bye tracker record")
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        ioServer.getDevice('tracker').sendMessage("key board escape")
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
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# --- Ending Routine "start_ET" ---
for thisComponent in start_ETComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# the Routine "start_ET" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# Repeat the centered_dot routine 6 times
for _ in range(6):  # Change to 6 if you want to repeat it 6 times
    # --- Prepare to start Routine "centered_dot" ---
    continueRoutine = True
    # update component parameters for each repeat
    # keep track of which components have finished
    centered_dotComponents = [dot_centered]
    for thisComponent in centered_dotComponents:
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

    # --- Run Routine "centered_dot" ---
    routineForceEnded = not continueRoutine
    while continueRoutine and routineTimer.getTime() < 5.0:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *dot_centered* updates
        
        # if dot_centered is starting this frame...
        if dot_centered.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            dot_centered.frameNStart = frameN  # exact frame index
            dot_centered.tStart = t  # local t and not account for scr refresh
            dot_centered.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(dot_centered, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'dot_centered.started')
            # update status
            dot_centered.status = STARTED
            dot_centered.setAutoDraw(True)
            ioServer.getDevice('tracker').sendMessage("ET: Start routine 'centered_dot'")
        
        # if dot_centered is active this frame...
        if dot_centered.status == STARTED:
            # update params
            pass
        
        # if dot_centered is stopping this frame...
        if dot_centered.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > dot_centered.tStartRefresh + 5.0-frameTolerance:
                # keep track of stop time/frame for later
                dot_centered.tStop = t  # not accounting for scr refresh
                dot_centered.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'dot_centered.stopped')
                # update status
                dot_centered.status = FINISHED
                dot_centered.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            ioServer.getDevice('tracker').sendMessage("key board escape")
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in centered_dotComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()

    # --- Ending Routine "centered_dot" ---
    for thisComponent in centered_dotComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
    if routineForceEnded:
        routineTimer.reset()
    else:
        routineTimer.addTime(-5.000000)
    thisExp.nextEntry()

# Repeat the dots routine 30 times
for _ in range(30):
    # --- Prepare to start Routine "dots" ---
    continueRoutine = True
    # Begin Routine
    shuffled_positions = list(positions.items())  # Convert dictionary to a list of tuples
    shuffle(shuffled_positions)  # Shuffle the positions
    current_position_index = 0  # Start with the first position
    total_positions = len(shuffled_positions)  # Track the total number of positions
    dot.pos, direction = shuffled_positions[current_position_index]  # Set initial dot position and direction
    time_of_last_change = 0  # Variable to store the time of the last position change
    continueRoutine = True  # Ensure the routine continues
    ioServer.getDevice('tracker').sendMessage("ET: Start routine 'dots'")
    # update component parameters for each repeat
    # keep track of which components have finished
    dotsComponents = [dot]
    for thisComponent in dotsComponents:
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
    
    # --- Run Routine "dots" ---
    routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *dot* updates
        
        # if dot is starting this frame...
        if dot.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            dot.frameNStart = frameN  # exact frame index
            dot.tStart = t  # local t and not account for scr refresh
            dot.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(dot, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'dot.started')
            # update status
            dot.status = STARTED
            dot.setAutoDraw(True)
        
        # if dot is active this frame...
        if dot.status == STARTED:
            # update params
            pass
        
        # Each Frame
        t = routineTimer.getTime()  # Get the current time in the routine

        # Check if enough time has passed since the last position change
        if t - time_of_last_change >= t_dot:
            current_position_index += 1  # Move to the next position
            if current_position_index < total_positions:  # Ensure the index is within bounds
                dot.pos, direction = shuffled_positions[current_position_index]  # Update the dot position and direction
                ioServer.getDevice('tracker').sendMessage(f"ET: dot moved {direction}!")  # Log the direction of movement
                time_of_last_change = t  # Reset the time of the last position change
            else:
                continueRoutine = False  # End the routine when all positions are visited
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            ioServer.getDevice('tracker').sendMessage("key board escape")
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in dotsComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "dots" ---
    for thisComponent in dotsComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # the Routine "dots" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    thisExp.nextEntry()

# Repeat the centered_dot routine 5 times
for _ in range(5):
    # --- Prepare to start Routine "centered_dot" ---
    continueRoutine = True
    # update component parameters for each repeat
    # keep track of which components have finished
    centered_dotComponents = [dot_centered]
    for thisComponent in centered_dotComponents:
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

    # --- Run Routine "centered_dot" ---
    routineForceEnded = not continueRoutine
    while continueRoutine and routineTimer.getTime() < 5.0:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *dot_centered* updates
        
        # if dot_centered is starting this frame...
        if dot_centered.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            dot_centered.frameNStart = frameN  # exact frame index
            dot_centered.tStart = t  # local t and not account for scr refresh
            dot_centered.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(dot_centered, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'dot_centered.started')
            # update status
            dot_centered.status = STARTED
            dot_centered.setAutoDraw(True)
            ioServer.getDevice('tracker').sendMessage("ET: Start routine 'centered_dot'")
        
        # if dot_centered is active this frame...
        if dot_centered.status == STARTED:
            # update params
            pass
        
        # if dot_centered is stopping this frame...
        if dot_centered.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > dot_centered.tStartRefresh + 5.0-frameTolerance:
                # keep track of stop time/frame for later
                dot_centered.tStop = t  # not accounting for scr refresh
                dot_centered.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'dot_centered.stopped')
                # update status
                dot_centered.status = FINISHED
                dot_centered.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            ioServer.getDevice('tracker').sendMessage("key board escape")
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in centered_dotComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()

    # --- Ending Routine "centered_dot" ---
    for thisComponent in centered_dotComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
    if routineForceEnded:
        routineTimer.reset()
    else:
        routineTimer.addTime(-5.000000)
    thisExp.nextEntry()

# --- Prepare to start Routine "end" ---
continueRoutine = True
# update component parameters for each repeat
key_resp_3.keys = []
key_resp_3.rt = []
_key_resp_3_allKeys = []
# Run 'Begin Routine' code from code_3
ioServer.getDevice('tracker').sendMessage("ET: Prepare to start routine 'end'")
# keep track of which components have finished
endComponents = [text, ET_stop, key_resp_3]
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
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *text* updates
    
    # if text is starting this frame...
    if text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        text.frameNStart = frameN  # exact frame index
        text.tStart = t  # local t and not account for scr refresh
        text.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(text, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'text.started')
        # update status
        text.status = STARTED
        text.setAutoDraw(True)
    
    # if text is active this frame...
    if text.status == STARTED:
        # update params
        pass
    
    # *ET_stop* updates
    if ET_stop.status == NOT_STARTED:
        ET_stop.frameNStart = frameN  # exact frame index
        ET_stop.tStart = t  # local t and not account for scr refresh
        ET_stop.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(ET_stop, 'tStartRefresh')  # time at next scr refresh
        ET_stop.status = STARTED
    
    # if ET_stop is stopping this frame...
    if ET_stop.status == STARTED:
        # is it time to stop? (based on local clock)
        if tThisFlip > 1.0-frameTolerance:
            # keep track of stop time/frame for later
            ET_stop.tStop = t  # not accounting for scr refresh
            ET_stop.tStopRefresh = tThisFlipGlobal  # on global time
            ET_stop.frameNStop = frameN  # exact frame index
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'ET_stop.stopped')
            # update status
            ET_stop.status = FINISHED
            ET_stop.stop()
    
    # *key_resp_3* updates
    waitOnFlip = False
    
    # if key_resp_3 is starting this frame...
    if key_resp_3.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        key_resp_3.frameNStart = frameN  # exact frame index
        key_resp_3.tStart = t  # local t and not account for scr refresh
        key_resp_3.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(key_resp_3, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'key_resp_3.started')
        # update status
        key_resp_3.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(key_resp_3.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(key_resp_3.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if key_resp_3.status == STARTED and not waitOnFlip:
        theseKeys = key_resp_3.getKeys(keyList=['t'], waitRelease=False)
        _key_resp_3_allKeys.extend(theseKeys)
        if len(_key_resp_3_allKeys):
            key_resp_3.keys = _key_resp_3_allKeys[-1].name  # just the last key pressed
            key_resp_3.rt = _key_resp_3_allKeys[-1].rt
            # a response ends the routine
            continueRoutine = False
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        ioServer.getDevice('tracker').sendMessage("key board escape")
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
if ET_stop.status != FINISHED:
    ET_stop.status = FINISHED
ioServer.getDevice('tracker').sendMessage("ET: eye-tracker stopped")
# check responses
if key_resp_3.keys in ['', [], None]:  # No response was made
    key_resp_3.keys = None
thisExp.addData('key_resp_3.keys',key_resp_3.keys)
if key_resp_3.keys != None:  # we had a response
    thisExp.addData('key_resp_3.rt', key_resp_3.rt)
thisExp.nextEntry()
# the Routine "end" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# Flip one final time so any remaining win.callOnFlip() 
# and win.timeOnFlip() tasks get executed before quitting
win.flip()

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename + '.csv', delim='auto')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
if eyetracker:
    eyetracker.setConnectionState(False)
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
