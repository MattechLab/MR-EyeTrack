# flake8: noqa
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This experiment was created using PsychoPy3 Experiment Builder 
# (v2022.3.0dev6),
#  on jeu 28 sep 2023 14:03:22
# If you publish work using this script the most relevant publication is:

# Peirce J, Gray JR, Simpson S, MacAskill M, 
# Höchenberger R, Sogo H, Kastman E,Lindeløv JK. (2019) 
# PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
# https://doi.org/10.3758/s13428-018-01193-y
# --- Import packages ---
# https://github.com/TheAxonLab/HCPh-fMRI-tasks/blob/27f5112ef476e35bb8689fb85c0a903c2c2c1cda/task-bht_bold.py#L758
# from psychopy import locale_setup
from psychopy import prefs
from psychopy import plugins
plugins.activatePlugins()
prefs.hardware['audioLib'] = 'ptb'
prefs.hardware['audioLatencyMode'] = '3'
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors, layout, iohub, hardware # 
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED, STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)
import socket
# from hcphsignals import signals
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
psychopyVersion = '2022.3.0dev6'
expName = 'fixation_dots_T1weighted'  # from the Builder filename that created this script
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
    originPath='/home/common/Desktop/4_Protocols/checkerboard_retina_fMRI/fixation_dots_T1weighted_lastrun.py',
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
    monitor='testMonitor', color=[0,0,0], colorSpace='rgb',
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
    'default_native_data_file_name': 'JB1',
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

# --- Initialize components for Routine "trail" --- #yiwei interested
waiting_trigger = visual.TextStim(win=win, name='waiting_trigger',
    text="The program is ready for the scanner trigger. Press 's' to proceed manually.",
    font='Open Sans',
    pos=(0, -0.4), height=0.12, wrapWidth=1.7, ori=0.0, 
    color='white', colorSpace='rgb', opacity=1.0, 
    languageStyle='LTR',
    depth=0.0);
key_resp = keyboard.Keyboard()
fix_desc = visual.TextStim(win=win, name='fix_desc',
    text='In this task you will see a grating patterns on the screen. \nPlease keep your eyes on the fixation circle like the one below.',
    font='Open Sans',
    pos=(0, 0.25), height=0.12, wrapWidth=1.0, ori=0.0, 
    color='white', colorSpace='rgb', opacity=1.0, 
    languageStyle='LTR',
    depth=-2.0);

# --- Initialize components for Routine "start_ET" ---
etRecord = hardware.eyetracker.EyetrackerControl(
    tracker=eyetracker,
    actionType='Start Only'
)#yiwei interested

# --- Initialize components for Routine "dots" ---
dot_1 = visual.ShapeStim(
    win=win, name='dot_1',
    size=(0.04, 0.04), vertices='circle',
    ori=0.0, pos=(0, 0), anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=[1, -0.2, -0.2], fillColor=[1.0000, -1.0000, -1.0000],
    opacity=None, depth=0.0, interpolate=True)
dot_2 = visual.ShapeStim(
    win=win, name='dot_2',
    size=(0.04, 0.04), vertices='circle',
    ori=0.0, pos=(0, 0), anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor='white', fillColor=[-1.0000, 0.2706, -0.9843],
    opacity=None, depth=-1.0, interpolate=True)
dot_3 = visual.ShapeStim(
    win=win, name='dot_3',
    size=(0.04, 0.04), vertices='circle',
    ori=0.0, pos=(0, 0), anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=[1,0,1], fillColor=[1.0000, 0.0745, -0.4667],
    opacity=None, depth=-2.0, interpolate=True)
dot_4 = visual.ShapeStim(
    win=win, name='dot_4',
    size=(0.04, 0.04), vertices='circle',
    ori=0.0, pos=(0, 0), anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=[-1,1,1], fillColor=[-1,1,1],
    opacity=None, depth=-3.0, interpolate=True)
dot_5 = visual.ShapeStim(
    win=win, name='dot_5',
    size=(0.04, 0.04), vertices='circle',
    ori=0.0, pos=(0, 0), anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=[1,1,1], fillColor=[0.4510, 0.4275, 0.4353],
    opacity=None, depth=-4.0, interpolate=True)
dot_6 = visual.ShapeStim(
    win=win, name='dot_6',
    size=(0.04, 0.04), vertices='circle',
    ori=0.0, pos=(0, 0), anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=[1,1,0.5], fillColor=[1.0000, 0.7647, -0.6235],
    opacity=None, depth=-5.0, interpolate=True)
dot_7 = visual.ShapeStim(
    win=win, name='dot_7',
    size=(0.04, 0.04), vertices='circle',
    ori=0.0, pos=(0, 0), anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=[-0.15, -0.1, 1], fillColor=[0.2235, -1.0000, 1.0000],
    opacity=None, depth=-6.0, interpolate=True)
dot_8 = visual.ShapeStim(
    win=win, name='dot_8',
    size=(0.04, 0.04), vertices='circle',
    ori=0.0, pos=(0, 0), anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=[0.5, 0.2, 0.2], fillColor=[0.3412, -0.1294, -0.5059],
    opacity=None, depth=-7.0, interpolate=True)
dot_9 = visual.ShapeStim(
    win=win, name='dot_9',
    size=(0.04, 0.04), vertices='circle',
    ori=0.0, pos=(0, 0), anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=[0.5, 1, 1], fillColor=[-0.4118, 0.5765, 0.8353],
    opacity=None, depth=-8.0, interpolate=True)
dot_10 = visual.ShapeStim(
    win=win, name='dot_10',
    size=(0.04, 0.04), vertices='circle',
    ori=0.0, pos=(0, 0), anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=[0,0,0], fillColor=[0,0,0],
    opacity=None, depth=-9.0, interpolate=True)

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
)#yiwei interested
key_resp_3 = keyboard.Keyboard()

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.Clock()  # to track time remaining of each (possibly non-slip) routine 
# define target for calibration_2
calibration_2Target = visual.TargetStim(win, 
    name='calibration_2Target',
    radius=0.15, fillColor=[-1, -1, -1], borderColor='black', lineWidth=2.0,
    innerRadius=0.07, innerFillColor='green', innerBorderColor='black', innerLineWidth=2.0,
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
calibration_2.run() #yiwei flag
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
        thisExp.timestampOnFlip(win, 'waiting_trigger.started') #yiwei flag
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
        thisExp.timestampOnFlip(win, 'key_resp.started') #yiwei flag
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
            continueRoutine = False #yiwei flag
    
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
    thisExp.addData('key_resp.rt', key_resp.rt) #yiwei flag
thisExp.nextEntry()
# the Routine "trail" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# --- Prepare to start Routine "start_ET" ---#yiwei interest
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

# --- Run Routine "start_ET" ---#yiwei interest
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
        thisExp.timestampOnFlip(win, 'etRecord.started') #yiwei flag
        # update status
        etRecord.status = STARTED
        # Run 'Begin Routine' code from code_channel2
        # send_message((signals.RUN | signals.ET_START_AND_STOP).to_bytes())
        ioServer.getDevice('tracker').sendMessage("Hello tracker")
        # ioServer.getDevice('tracker').sendMessage("start visual stimuli")
    
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
            # ioServer.getDevice('tracker').sendMessage("stop visual stimuli")
    
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
# make sure the eyetracker recording stops
if etRecord.status != FINISHED:
    etRecord.status = FINISHED
# the Routine "start_ET" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
T1w_LIBRE = data.TrialHandler(nReps=720.0, method='random', 
    extraInfo=expInfo, originPath=-1,
    trialList=[None],
    seed=None, name='T1w_LIBRE')
thisExp.addLoop(T1w_LIBRE)  # add the loop to the experiment
thisT1w_LIBRE = T1w_LIBRE.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisT1w_LIBRE.rgb)
if thisT1w_LIBRE != None:
    ioServer.getDevice('tracker').sendMessage("T1w_LIBRE stimuli start")
    for paramName in thisT1w_LIBRE:
        exec('{} = thisT1w_LIBRE[paramName]'.format(paramName))

for thisT1w_LIBRE in T1w_LIBRE:
    currentLoop = T1w_LIBRE
    # abbreviate parameter names if possible (e.g. rgb = thisT1w_LIBRE.rgb)
    if thisT1w_LIBRE != None:
        for paramName in thisT1w_LIBRE:
            exec('{} = thisT1w_LIBRE[paramName]'.format(paramName))
    
    # --- Prepare to start Routine "dots" ---
    continueRoutine = True
    # update component parameters for each repeat
    dot_2.setLineColor([0,1,0])
    # keep track of which components have finished
    dotsComponents = [dot_1, dot_2, dot_3, dot_4, dot_5, dot_6, dot_7, dot_8, dot_9, dot_10]
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
    while continueRoutine and routineTimer.getTime() < 1.0:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *dot_1* updates
        
        # if dot_1 is starting this frame...
        if dot_1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            dot_1.frameNStart = frameN  # exact frame index
            dot_1.tStart = t  # local t and not account for scr refresh
            dot_1.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(dot_1, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'dot_1.started')
            # update status
            dot_1.status = STARTED
            dot_1.setAutoDraw(True)
        
        # if dot_1 is active this frame...
        if dot_1.status == STARTED:
            # update params
            pass
        
        # if dot_1 is stopping this frame...
        if dot_1.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > dot_1.tStartRefresh + 0.1-frameTolerance:
                # keep track of stop time/frame for later
                dot_1.tStop = t  # not accounting for scr refresh
                dot_1.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'dot_1.stopped')
                # update status
                dot_1.status = FINISHED
                dot_1.setAutoDraw(False)
        
        # *dot_2* updates
        
        # if dot_2 is starting this frame...
        if dot_2.status == NOT_STARTED and tThisFlip >= 0.1-frameTolerance:
            # keep track of start time/frame for later
            dot_2.frameNStart = frameN  # exact frame index
            dot_2.tStart = t  # local t and not account for scr refresh
            dot_2.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(dot_2, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'dot_2.started')
            # update status
            dot_2.status = STARTED
            dot_2.setAutoDraw(True)
        
        # if dot_2 is active this frame...
        if dot_2.status == STARTED:
            # update params
            pass
        
        # if dot_2 is stopping this frame...
        if dot_2.status == STARTED:
            # is it time to stop? (based on local clock)
            if tThisFlip > 0.2-frameTolerance:
                # keep track of stop time/frame for later
                dot_2.tStop = t  # not accounting for scr refresh
                dot_2.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'dot_2.stopped')
                # update status
                dot_2.status = FINISHED
                dot_2.setAutoDraw(False)
        
        # *dot_3* updates
        
        # if dot_3 is starting this frame...
        if dot_3.status == NOT_STARTED and tThisFlip >= 0.2-frameTolerance:
            # keep track of start time/frame for later
            dot_3.frameNStart = frameN  # exact frame index
            dot_3.tStart = t  # local t and not account for scr refresh
            dot_3.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(dot_3, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'dot_3.started')
            # update status
            dot_3.status = STARTED
            dot_3.setAutoDraw(True)
        
        # if dot_3 is active this frame...
        if dot_3.status == STARTED:
            # update params
            pass
        
        # if dot_3 is stopping this frame...
        if dot_3.status == STARTED:
            # is it time to stop? (based on local clock)
            if tThisFlip > 0.3-frameTolerance:
                # keep track of stop time/frame for later
                dot_3.tStop = t  # not accounting for scr refresh
                dot_3.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'dot_3.stopped')
                # update status
                dot_3.status = FINISHED
                dot_3.setAutoDraw(False)
        
        # *dot_4* updates
        
        # if dot_4 is starting this frame...
        if dot_4.status == NOT_STARTED and tThisFlip >= 0.3-frameTolerance:
            # keep track of start time/frame for later
            dot_4.frameNStart = frameN  # exact frame index
            dot_4.tStart = t  # local t and not account for scr refresh
            dot_4.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(dot_4, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'dot_4.started')
            # update status
            dot_4.status = STARTED
            dot_4.setAutoDraw(True)
        
        # if dot_4 is active this frame...
        if dot_4.status == STARTED:
            # update params
            pass
        
        # if dot_4 is stopping this frame...
        if dot_4.status == STARTED:
            # is it time to stop? (based on local clock)
            if tThisFlip > 0.4-frameTolerance:
                # keep track of stop time/frame for later
                dot_4.tStop = t  # not accounting for scr refresh
                dot_4.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'dot_4.stopped')
                # update status
                dot_4.status = FINISHED
                dot_4.setAutoDraw(False)
        
        # *dot_5* updates
        
        # if dot_5 is starting this frame...
        if dot_5.status == NOT_STARTED and tThisFlip >= 0.4-frameTolerance:
            # keep track of start time/frame for later
            dot_5.frameNStart = frameN  # exact frame index
            dot_5.tStart = t  # local t and not account for scr refresh
            dot_5.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(dot_5, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'dot_5.started')
            # update status
            dot_5.status = STARTED
            dot_5.setAutoDraw(True)
        
        # if dot_5 is active this frame...
        if dot_5.status == STARTED:
            # update params
            pass
        
        # if dot_5 is stopping this frame...
        if dot_5.status == STARTED:
            # is it time to stop? (based on local clock)
            if tThisFlip > 0.5-frameTolerance:
                # keep track of stop time/frame for later
                dot_5.tStop = t  # not accounting for scr refresh
                dot_5.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'dot_5.stopped')
                # update status
                dot_5.status = FINISHED
                dot_5.setAutoDraw(False)
        
        # *dot_6* updates
        
        # if dot_6 is starting this frame...
        if dot_6.status == NOT_STARTED and tThisFlip >= 0.5-frameTolerance:
            # keep track of start time/frame for later
            dot_6.frameNStart = frameN  # exact frame index
            dot_6.tStart = t  # local t and not account for scr refresh
            dot_6.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(dot_6, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'dot_6.started')
            # update status
            dot_6.status = STARTED
            dot_6.setAutoDraw(True)
        
        # if dot_6 is active this frame...
        if dot_6.status == STARTED:
            # update params
            pass
        
        # if dot_6 is stopping this frame...
        if dot_6.status == STARTED:
            # is it time to stop? (based on local clock)
            if tThisFlip > 0.6-frameTolerance:
                # keep track of stop time/frame for later
                dot_6.tStop = t  # not accounting for scr refresh
                dot_6.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'dot_6.stopped')
                # update status
                dot_6.status = FINISHED
                dot_6.setAutoDraw(False)
        
        # *dot_7* updates
        
        # if dot_7 is starting this frame...
        if dot_7.status == NOT_STARTED and tThisFlip >= 0.6-frameTolerance:
            # keep track of start time/frame for later
            dot_7.frameNStart = frameN  # exact frame index
            dot_7.tStart = t  # local t and not account for scr refresh
            dot_7.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(dot_7, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'dot_7.started')
            # update status
            dot_7.status = STARTED
            dot_7.setAutoDraw(True)
        
        # if dot_7 is active this frame...
        if dot_7.status == STARTED:
            # update params
            pass
        
        # if dot_7 is stopping this frame...
        if dot_7.status == STARTED:
            # is it time to stop? (based on local clock)
            if tThisFlip > 0.7-frameTolerance:
                # keep track of stop time/frame for later
                dot_7.tStop = t  # not accounting for scr refresh
                dot_7.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'dot_7.stopped')
                # update status
                dot_7.status = FINISHED
                dot_7.setAutoDraw(False)
        
        # *dot_8* updates
        
        # if dot_8 is starting this frame...
        if dot_8.status == NOT_STARTED and tThisFlip >= 0.7-frameTolerance:
            # keep track of start time/frame for later
            dot_8.frameNStart = frameN  # exact frame index
            dot_8.tStart = t  # local t and not account for scr refresh
            dot_8.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(dot_8, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'dot_8.started')
            # update status
            dot_8.status = STARTED
            dot_8.setAutoDraw(True)
        
        # if dot_8 is active this frame...
        if dot_8.status == STARTED:
            # update params
            pass
        
        # if dot_8 is stopping this frame...
        if dot_8.status == STARTED:
            # is it time to stop? (based on local clock)
            if tThisFlip > 0.8-frameTolerance:
                # keep track of stop time/frame for later
                dot_8.tStop = t  # not accounting for scr refresh
                dot_8.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'dot_8.stopped')
                # update status
                dot_8.status = FINISHED
                dot_8.setAutoDraw(False)
        
        # *dot_9* updates
        
        # if dot_9 is starting this frame...
        if dot_9.status == NOT_STARTED and tThisFlip >= 0.8-frameTolerance:
            # keep track of start time/frame for later
            dot_9.frameNStart = frameN  # exact frame index
            dot_9.tStart = t  # local t and not account for scr refresh
            dot_9.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(dot_9, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'dot_9.started')
            # update status
            dot_9.status = STARTED
            dot_9.setAutoDraw(True)
        
        # if dot_9 is active this frame...
        if dot_9.status == STARTED:
            # update params
            pass
        
        # if dot_9 is stopping this frame...
        if dot_9.status == STARTED:
            # is it time to stop? (based on local clock)
            if tThisFlip > 0.9-frameTolerance:
                # keep track of stop time/frame for later
                dot_9.tStop = t  # not accounting for scr refresh
                dot_9.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'dot_9.stopped')
                # update status
                dot_9.status = FINISHED
                dot_9.setAutoDraw(False)
        
        # *dot_10* updates
        
        # if dot_10 is starting this frame...
        if dot_10.status == NOT_STARTED and tThisFlip >= 0.9-frameTolerance:
            # keep track of start time/frame for later
            dot_10.frameNStart = frameN  # exact frame index
            dot_10.tStart = t  # local t and not account for scr refresh
            dot_10.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(dot_10, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'dot_10.started')
            # update status
            dot_10.status = STARTED
            dot_10.setAutoDraw(True)
        
        # if dot_10 is active this frame...
        if dot_10.status == STARTED:
            # update params
            pass
        
        # if dot_10 is stopping this frame...
        if dot_10.status == STARTED:
            # is it time to stop? (based on local clock)
            if tThisFlip > 1.0-frameTolerance:
                # keep track of stop time/frame for later
                dot_10.tStop = t  # not accounting for scr refresh
                dot_10.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'dot_10.stopped')
                # update status
                dot_10.status = FINISHED
                dot_10.setAutoDraw(False)
        
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
    # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
    if routineForceEnded:
        ioServer.getDevice('tracker').sendMessage("Routine T1w ended")
        routineTimer.reset()
    else:
        routineTimer.addTime(-1.000000)
    thisExp.nextEntry()
    
# completed 720.0 repeats of 'T1w_LIBRE'


# --- Prepare to start Routine "end" ---
continueRoutine = True
# update component parameters for each repeat
key_resp_3.keys = []
key_resp_3.rt = []
_key_resp_3_allKeys = []
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
    
    # if ET_stop is stopping this frame...
    if ET_stop.status == STARTED:
        # is it time to stop? (based on local clock)
        if tThisFlip > 1.0-frameTolerance:
            # keep track of stop time/frame for later
            ET_stop.tStop = t  # not accounting for scr refresh
            ET_stop.frameNStop = frameN  # exact frame index
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'ET_stop.stopped')
            # update status
            ET_stop.status = FINISHED
            ioServer.getDevice('tracker').sendMessage("stop tracker block")
    
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
# check responses
if key_resp_3.keys in ['', [], None]:  # No response was made
    key_resp_3.keys = None
thisExp.addData('key_resp_3.keys',key_resp_3.keys)
if key_resp_3.keys != None:  # we had a response
    thisExp.addData('key_resp_3.rt', key_resp_3.rt)
thisExp.nextEntry()
# the Routine "end" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# --- End experiment ---
# Flip one final time so any remaining win.callOnFlip() 
# and win.timeOnFlip() tasks get executed before quitting
win.flip()

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename+'.csv', delim='auto')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
if eyetracker:
    eyetracker.setConnectionState(False)
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()

