#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2023.2.3),
    on ven 27 oct 2023 16:45:19
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

# --- Import packages ---
from psychopy import locale_setup
from psychopy import prefs
from psychopy import plugins
plugins.activatePlugins()
prefs.hardware['audioLib'] = 'ptb'
prefs.hardware['audioLatencyMode'] = '0'
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors, layout, iohub, hardware
from psychopy.tools import environmenttools
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER, priority)

import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle, choice as randchoice
import os  # handy system and path functions
import sys  # to get file system encoding

import psychopy.iohub as io
from psychopy.hardware import keyboard

# Run 'Before Experiment' code from code_7
import socket
from hcphsignals import signals
from datetime import datetime
import pwd
def send_message(message, addr="localhost", port=2023):
    client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    client_socket.connect((addr, port))
    client_socket.sendall(message)
    client_socket.close()
current_user = pwd.getpwuid(os.geteuid()).pw_name
folder_path = os.path.dirname('/home/'+current_user +'/workspace/HCPh-data/Psychopy/session-'+data.getDateStr(format="%Y-%m-%d")+'/')
if not os.path.exists(folder_path):
    os.makedirs(folder_path)
ename = 'bht_'+data.getDateStr(format="%Y-%m-%d")
trial = 0
for file in os.listdir(folder_path):
    filename_hour = datetime.strptime(file.split('_')[2],'%Hh%M.%S.%f')
    current_time = datetime.now()

    filename_hour_components = (filename_hour.hour, filename_hour.minute)
    current_time_components = (current_time.hour, current_time.minute)
    hour_difference = abs((current_time_components[0] - filename_hour_components[0]) +
                     (current_time_components[1] - filename_hour_components[1]) / 60.0)

    if ename in file and '.log' in file and hour_difference <1.5:
        trial += 1
# --- Setup global variables (available in all functions) ---
# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
# Store info about the experiment session
psychopyVersion = '2023.2.3'
expName = 'bht'  # from the Builder filename that created this script
expInfo = {
    'trial': trial,
    'session': '9999',
    'date': data.getDateStr(),  # add a simple timestamp
    'expName': expName,
    'psychopyVersion': psychopyVersion,
}


def showExpInfoDlg(expInfo):
    """
    Show participant info dialog.
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    
    Returns
    ==========
    dict
        Information about this experiment.
    """
    # temporarily remove keys which the dialog doesn't need to show
    poppedKeys = {
        'date': expInfo.pop('date', data.getDateStr()),
        'expName': expInfo.pop('expName', expName),
        'psychopyVersion': expInfo.pop('psychopyVersion', psychopyVersion),
    }
    # show participant info dialog
    dlg = gui.DlgFromDict(dictionary=expInfo, sortKeys=False, title=expName)
    if dlg.OK == False:
        core.quit()  # user pressed cancel
    # restore hidden keys
    expInfo.update(poppedKeys)
    # return expInfo
    return expInfo


def setupData(expInfo, dataDir=None):
    """
    Make an ExperimentHandler to handle trials and saving.
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    dataDir : Path, str or None
        Folder to save the data to, leave as None to create a folder in the current directory.    
    Returns
    ==========
    psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    """
    
    # data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
    if dataDir is None:
        dataDir = _thisDir
    filename = u'../HCPh-data/Psychopy/session-'+data.getDateStr(format="%Y-%m-%d")+'/%s_%s_%s_session_%s' %(expName, expInfo['date'],expInfo['trial'],expInfo['session'])
    # make sure filename is relative to dataDir
    if os.path.isabs(filename):
        dataDir = os.path.commonprefix([dataDir, filename])
        filename = os.path.relpath(filename, dataDir)
    
    # an ExperimentHandler isn't essential but helps with data saving
    thisExp = data.ExperimentHandler(
        name=expName, version='',
        extraInfo=expInfo, runtimeInfo=None,
        originPath='/home/common/workspace/HCPh-fMRI-tasks/task-bht_bold.py',
        savePickle=True, saveWideText=True,
        dataFileName=dataDir + os.sep + filename, sortColumns='time'
    )
    thisExp.setPriority('thisRow.t', priority.CRITICAL)
    thisExp.setPriority('expName', priority.LOW)
    # return experiment handler
    return thisExp


def setupLogging(filename):
    """
    Setup a log file and tell it what level to log at.
    
    Parameters
    ==========
    filename : str or pathlib.Path
        Filename to save log file and data files as, doesn't need an extension.
    
    Returns
    ==========
    psychopy.logging.LogFile
        Text stream to receive inputs from the logging system.
    """
    # this outputs to the screen, not a file
    logging.console.setLevel(logging.INFO)
    # save a log file for detail verbose info
    logFile = logging.LogFile(filename+'.log', level=logging.INFO)
    
    return logFile


def setupWindow(expInfo=None, win=None):
    """
    Setup the Window
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    win : psychopy.visual.Window
        Window to setup - leave as None to create a new window.
    
    Returns
    ==========
    psychopy.visual.Window
        Window in which to run this experiment.
    """
    if win is None:
        # if not given a window to setup, make one
        win = visual.Window(
            size=[800, 600], fullscr=True, screen=1,
            winType='pyglet', allowStencil=False,
            monitor='testMonitor', color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb',
            backgroundImage='', backgroundFit='none',
            blendMode='avg', useFBO=True,
            units='norm'
        )
        if expInfo is not None:
            # store frame rate of monitor if we can measure it
            expInfo['frameRate'] = win.getActualFrameRate()
    else:
        # if we have a window, just set the attributes which are safe to set
        win.color = [-1.0000, -1.0000, -1.0000]
        win.colorSpace = 'rgb'
        win.backgroundImage = ''
        win.backgroundFit = 'none'
        win.units = 'norm'
    win.mouseVisible = False
    win.hideMessage()
    return win


def setupInputs(expInfo, thisExp, win):
    """
    Setup whatever inputs are available (mouse, keyboard, eyetracker, etc.)
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    win : psychopy.visual.Window
        Window in which to run this experiment.
    Returns
    ==========
    dict
        Dictionary of input devices by name.
    """
    # --- Setup input devices ---
    inputs = {}
    ioConfig = {}
    
    # Setup eyetracking
    ioConfig['eyetracker.hw.sr_research.eyelink.EyeTracker'] = {
        'name': 'tracker',
        'model_name': 'EYELINK 1000 LONG RANGE',
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
    ioServer = io.launchHubServer(window=win, experiment_code='bht', session_code=ioSession, datastore_name=thisExp.dataFileName, **ioConfig)
    eyetracker = ioServer.getDevice('tracker')
    
    # create a default keyboard (e.g. to check for escape)
    defaultKeyboard = keyboard.Keyboard(backend='iohub')
    # return inputs dict
    return {
        'ioServer': ioServer,
        'defaultKeyboard': defaultKeyboard,
        'eyetracker': eyetracker,
    }

def pauseExperiment(thisExp, inputs=None, win=None, timers=[], playbackComponents=[]):
    """
    Pause this experiment, preventing the flow from advancing to the next routine until resumed.
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    inputs : dict
        Dictionary of input devices by name.
    win : psychopy.visual.Window
        Window for this experiment.
    timers : list, tuple
        List of timers to reset once pausing is finished.
    playbackComponents : list, tuple
        List of any components with a `pause` method which need to be paused.
    """
    # if we are not paused, do nothing
    if thisExp.status != PAUSED:
        return
    
    # pause any playback components
    for comp in playbackComponents:
        comp.pause()
    # prevent components from auto-drawing
    win.stashAutoDraw()
    # run a while loop while we wait to unpause
    while thisExp.status == PAUSED:
        # make sure we have a keyboard
        if inputs is None:
            inputs = {
                'defaultKeyboard': keyboard.Keyboard(backend='ioHub')
            }
        # check for quit (typically the Esc key)
        if inputs['defaultKeyboard'].getKeys(keyList=['escape']):
            endExperiment(thisExp, win=win, inputs=inputs)
        # flip the screen
        win.flip()
    # if stop was requested while paused, quit
    if thisExp.status == FINISHED:
        endExperiment(thisExp, inputs=inputs, win=win)
    # resume any playback components
    for comp in playbackComponents:
        comp.play()
    # restore auto-drawn components
    win.retrieveAutoDraw()
    # reset any timers
    for timer in timers:
        timer.reset()


def run(expInfo, thisExp, win, inputs, globalClock=None, thisSession=None):
    """
    Run the experiment flow.
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    psychopy.visual.Window
        Window in which to run this experiment.
    inputs : dict
        Dictionary of input devices by name.
    globalClock : psychopy.core.clock.Clock or None
        Clock to get global time from - supply None to make a new one.
    thisSession : psychopy.session.Session or None
        Handle of the Session object this experiment is being run from, if any.
    """
    # mark experiment as started
    thisExp.status = STARTED
    # make sure variables created by exec are available globally
    exec = environmenttools.setExecEnvironment(globals())
    # get device handles from dict of input devices
    ioServer = inputs['ioServer']
    defaultKeyboard = inputs['defaultKeyboard']
    eyetracker = inputs['eyetracker']
    # make sure we're running in the directory for this experiment
    os.chdir(_thisDir)
    # get filename from ExperimentHandler for convenience
    filename = thisExp.dataFileName
    frameTolerance = 0.001  # how close to onset before 'same' frame
    endExpNow = False  # flag for 'escape' or other condition => quit the exp
    # get frame duration from frame rate in expInfo
    if 'frameRate' in expInfo and expInfo['frameRate'] is not None:
        frameDur = 1.0 / round(expInfo['frameRate'])
    else:
        frameDur = 1.0 / 60.0  # could not measure, so guess
    
    # Start Code - component code to be run after the window creation
    
    # --- Initialize components for Routine "intro" ---
    text_2 = visual.TextStim(win=win, name='text_2',
        text='The breath-holding task will start shortly.\nPlease breathe at a comfortable pace.',
        font='Open Sans',
        pos=(0, 0), height=0.08, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=0.0);
    key_resp = keyboard.Keyboard()
    
    # --- Initialize components for Routine "instructions" ---
    directions = visual.TextStim(win=win, name='directions',
        text='Breath-holding task\n\nBreathe in and out through your nose during the whole experiment.\n\nPlease DO NOT FOLLOW THE COLOR CUES IN THE FIRST TRIAL (and only the first trial).\n \nA message will let you know when to start holding your breath with the red sign.',
        font='Open Sans',
        pos=(0, 0), height=0.08, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=0.0);
    
    # --- Initialize components for Routine "start_et" ---
    etRecord = hardware.eyetracker.EyetrackerControl(
        tracker=eyetracker,
        actionType='Start Only'
    )
    
    # --- Initialize components for Routine "breath_in" ---
    polygon_4 = visual.Rect(
        win=win, name='polygon_4',
        width=(0.5, 0.02)[0], height=(0.5, 0.02)[1],
        ori=0.0, pos=(0, 0), anchor='center',
        lineWidth=1.0,     colorSpace='rgb',  lineColor=[-1.0000, 0.0039, -1.0000], fillColor=[-1.0000, 0.0039, -1.0000],
        opacity=None, depth=0.0, interpolate=True)
    
    # --- Initialize components for Routine "breath_out" ---
    polygon1 = visual.Rect(
        win=win, name='polygon1',
        width=(0.5, 0.02)[0], height=(0.5, 0.02)[1],
        ori=0.0, pos=(0, 0), anchor='center',
        lineWidth=1.0,     colorSpace='rgb',  lineColor=[1.0000, 0.2941, -1.0000], fillColor=[1.0000,0.2941, -1.0000],
        opacity=None, depth=0.0, interpolate=True)
    
    # --- Initialize components for Routine "breath_in_last" ---
    polygon_6 = visual.Rect(
        win=win, name='polygon_6',
        width=(0.5, 0.02)[0], height=(0.5, 0.02)[1],
        ori=0.0, pos=(0, 0), anchor='center',
        lineWidth=1.0,     colorSpace='rgb',  lineColor=[-0.4000, 0.4539, -0.4000], fillColor=[-0.4000, 0.4539, -0.4000],
        opacity=None, depth=0.0, interpolate=True)
    
    # --- Initialize components for Routine "breath_out_last" ---
    polygon_8 = visual.Rect(
        win=win, name='polygon_8',
        width=(0.5, 0.02)[0], height=(0.5, 0.02)[1],
        ori=0.0, pos=(0, 0), anchor='center',
        lineWidth=1.0,     colorSpace='rgb',  lineColor=[1.0000, 0.1000, -1.0000], fillColor=[1.0000, 0.1000, -1.0000],
        opacity=None, depth=0.0, interpolate=True)
    
    # --- Initialize components for Routine "hold_test" ---
    bh_body = visual.Rect(
        win=win, name='bh_body',
        width=(0.5, 0.02)[0], height=(0.5, 0.02)[1],
        ori=0.0, pos=(0, 0), anchor='center',
        lineWidth=1.0,     colorSpace='rgb',  lineColor=[1.0000, -1.0000, -1.0000], fillColor=[1.0000, -1.0000, -1.0000],
        opacity=None, depth=0.0, interpolate=True)
    bh_end = visual.Rect(
        win=win, name='bh_end',
        width=(0.5, 0.02)[0], height=(0.5, 0.02)[1],
        ori=0.0, pos=(0, 0), anchor='center',
        lineWidth=1.0,     colorSpace='rgb',  lineColor=[1.0000, -0.1200, -0.1200], fillColor=[1.0000, -0.1200, -0.1200],
        opacity=None, depth=-1.0, interpolate=True)
    end_trial_msg = visual.TextStim(win=win, name='end_trial_msg',
        text='End of test trial\n\nPlease engage in the task and HOLD YOUR BREATH during the red colored rectangles.',
        font='Open Sans',
        pos=(0, 0), height=0.05, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-2.0);
    
    # --- Initialize components for Routine "breath_in" ---
    polygon_4 = visual.Rect(
        win=win, name='polygon_4',
        width=(0.5, 0.02)[0], height=(0.5, 0.02)[1],
        ori=0.0, pos=(0, 0), anchor='center',
        lineWidth=1.0,     colorSpace='rgb',  lineColor=[-1.0000, 0.0039, -1.0000], fillColor=[-1.0000, 0.0039, -1.0000],
        opacity=None, depth=0.0, interpolate=True)
    
    # --- Initialize components for Routine "breath_out" ---
    polygon1 = visual.Rect(
        win=win, name='polygon1',
        width=(0.5, 0.02)[0], height=(0.5, 0.02)[1],
        ori=0.0, pos=(0, 0), anchor='center',
        lineWidth=1.0,     colorSpace='rgb',  lineColor=[1.0000, 0.2941, -1.0000], fillColor=[1.0000,0.2941, -1.0000],
        opacity=None, depth=0.0, interpolate=True)
    
    # --- Initialize components for Routine "breath_in_last" ---
    polygon_6 = visual.Rect(
        win=win, name='polygon_6',
        width=(0.5, 0.02)[0], height=(0.5, 0.02)[1],
        ori=0.0, pos=(0, 0), anchor='center',
        lineWidth=1.0,     colorSpace='rgb',  lineColor=[-0.4000, 0.4539, -0.4000], fillColor=[-0.4000, 0.4539, -0.4000],
        opacity=None, depth=0.0, interpolate=True)
    
    # --- Initialize components for Routine "breath_out_last" ---
    polygon_8 = visual.Rect(
        win=win, name='polygon_8',
        width=(0.5, 0.02)[0], height=(0.5, 0.02)[1],
        ori=0.0, pos=(0, 0), anchor='center',
        lineWidth=1.0,     colorSpace='rgb',  lineColor=[1.0000, 0.1000, -1.0000], fillColor=[1.0000, 0.1000, -1.0000],
        opacity=None, depth=0.0, interpolate=True)
    
    # --- Initialize components for Routine "hold" ---
    bh_body_2 = visual.Rect(
        win=win, name='bh_body_2',
        width=(0.5, 0.02)[0], height=(0.5, 0.02)[1],
        ori=0.0, pos=(0, 0), anchor='center',
        lineWidth=1.0,     colorSpace='rgb',  lineColor=[1.0000, -1.0000, -1.0000], fillColor=[1.0000, -1.0000, -1.0000],
        opacity=None, depth=0.0, interpolate=True)
    bh_end_2 = visual.Rect(
        win=win, name='bh_end_2',
        width=(0.5, 0.02)[0], height=(0.5, 0.02)[1],
        ori=0.0, pos=(0, 0), anchor='center',
        lineWidth=1.0,     colorSpace='rgb',  lineColor=[1.0000, -0.1200, -0.1200], fillColor=[1.0000, -0.1200, -0.1200],
        opacity=None, depth=-1.0, interpolate=True)
    bh_end_3 = visual.Rect(
        win=win, name='bh_end_3',
        width=(0.5, 0.02)[0], height=(0.5, 0.02)[1],
        ori=0.0, pos=(0, 0), anchor='center',
        lineWidth=1.0,     colorSpace='rgb',  lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-2.0, interpolate=True)
    
    # --- Initialize components for Routine "done" ---
    finished = visual.TextStim(win=win, name='finished',
        text='Finished!\n\nPlease keep breathing at a comfortable pace',
        font='Open Sans',
        pos=(0, 0), height=0.08, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=0.0);
    
    # --- Initialize components for Routine "stop_et" ---
    etRecord_2 = hardware.eyetracker.EyetrackerControl(
        tracker=eyetracker,
        actionType='Stop Only'
    )
    
    # create some handy timers
    if globalClock is None:
        globalClock = core.Clock()  # to track the time since experiment started
    if ioServer is not None:
        ioServer.syncClock(globalClock)
    logging.setDefaultClock(globalClock)
    routineTimer = core.Clock()  # to track time remaining of each (possibly non-slip) routine
    win.flip()  # flip window to reset last flip timer
    # store the exact time the global clock started
    expInfo['expStart'] = data.getDateStr(format='%Y-%m-%d %Hh%M.%S.%f %z', fractionalSecondDigits=6)
    # define target for calibration_2
    calibration_2Target = visual.TargetStim(win, 
        name='calibration_2Target',
        radius=0.075, fillColor='', borderColor='black', lineWidth=2.0,
        innerRadius=0.035, innerFillColor='green', innerBorderColor='black', innerLineWidth=2.0,
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
    
    # --- Prepare to start Routine "intro" ---
    continueRoutine = True
    # update component parameters for each repeat
    thisExp.addData('intro.started', globalClock.getTime())
    key_resp.keys = []
    key_resp.rt = []
    _key_resp_allKeys = []
    # keep track of which components have finished
    introComponents = [text_2, key_resp]
    for thisComponent in introComponents:
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
    
    # --- Run Routine "intro" ---
    routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *text_2* updates
        
        # if text_2 is starting this frame...
        if text_2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_2.frameNStart = frameN  # exact frame index
            text_2.tStart = t  # local t and not account for scr refresh
            text_2.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_2, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'text_2.started')
            # update status
            text_2.status = STARTED
            text_2.setAutoDraw(True)
        
        # if text_2 is active this frame...
        if text_2.status == STARTED:
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
            theseKeys = key_resp.getKeys(keyList=['t','s'], ignoreKeys=["escape"], waitRelease=False)
            _key_resp_allKeys.extend(theseKeys)
            if len(_key_resp_allKeys):
                key_resp.keys = _key_resp_allKeys[-1].name  # just the last key pressed
                key_resp.rt = _key_resp_allKeys[-1].rt
                key_resp.duration = _key_resp_allKeys[-1].duration
                # a response ends the routine
                continueRoutine = False
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, inputs=inputs, win=win)
            return
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in introComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "intro" ---
    for thisComponent in introComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    thisExp.addData('intro.stopped', globalClock.getTime())
    # check responses
    if key_resp.keys in ['', [], None]:  # No response was made
        key_resp.keys = None
    thisExp.addData('key_resp.keys',key_resp.keys)
    if key_resp.keys != None:  # we had a response
        thisExp.addData('key_resp.rt', key_resp.rt)
        thisExp.addData('key_resp.duration', key_resp.duration)
    thisExp.nextEntry()
    # the Routine "intro" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # --- Prepare to start Routine "instructions" ---
    continueRoutine = True
    # update component parameters for each repeat
    thisExp.addData('instructions.started', globalClock.getTime())
    # keep track of which components have finished
    instructionsComponents = [directions]
    for thisComponent in instructionsComponents:
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
    
    # --- Run Routine "instructions" ---
    routineForceEnded = not continueRoutine
    while continueRoutine and routineTimer.getTime() < 10.0:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *directions* updates
        
        # if directions is starting this frame...
        if directions.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            directions.frameNStart = frameN  # exact frame index
            directions.tStart = t  # local t and not account for scr refresh
            directions.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(directions, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'directions.started')
            # update status
            directions.status = STARTED
            directions.setAutoDraw(True)
        
        # if directions is active this frame...
        if directions.status == STARTED:
            # update params
            pass
        
        # if directions is stopping this frame...
        if directions.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > directions.tStartRefresh + 10.0-frameTolerance:
                # keep track of stop time/frame for later
                directions.tStop = t  # not accounting for scr refresh
                directions.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'directions.stopped')
                # update status
                directions.status = FINISHED
                directions.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, inputs=inputs, win=win)
            return
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in instructionsComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "instructions" ---
    for thisComponent in instructionsComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    thisExp.addData('instructions.stopped', globalClock.getTime())
    # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
    if routineForceEnded:
        routineTimer.reset()
    else:
        routineTimer.addTime(-10.000000)
    
    # --- Prepare to start Routine "start_et" ---
    continueRoutine = True
    # update component parameters for each repeat
    thisExp.addData('start_et.started', globalClock.getTime())
    # Run 'Begin Routine' code from code_channel2
    send_message((signals.RUN | signals.ET_START_AND_STOP).to_bytes())
    ioServer.getDevice('tracker').sendMessage("Hello bht")
    ioServer.getDevice('tracker').sendMessage("start mock block")
    
    # keep track of which components have finished
    start_etComponents = [etRecord]
    for thisComponent in start_etComponents:
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
    
    # --- Run Routine "start_et" ---
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
        if etRecord.status == NOT_STARTED and t >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            etRecord.frameNStart = frameN  # exact frame index
            etRecord.tStart = t  # local t and not account for scr refresh
            etRecord.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(etRecord, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.addData('etRecord.started', t)
            # update status
            etRecord.status = STARTED
        
        # if etRecord is stopping this frame...
        if etRecord.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > etRecord.tStartRefresh + 0-frameTolerance:
                # keep track of stop time/frame for later
                etRecord.tStop = t  # not accounting for scr refresh
                etRecord.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.addData('etRecord.stopped', t)
                # update status
                etRecord.status = FINISHED
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, inputs=inputs, win=win)
            return
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in start_etComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "start_et" ---
    for thisComponent in start_etComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    thisExp.addData('start_et.stopped', globalClock.getTime())
    # make sure the eyetracker recording stops
    if etRecord.status != FINISHED:
        etRecord.status = FINISHED
    # the Routine "start_et" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    normal_respiratory_rate = data.TrialHandler(nReps=4.0, method='sequential', 
        extraInfo=expInfo, originPath=-1,
        trialList=[None],
        seed=None, name='normal_respiratory_rate')
    thisExp.addLoop(normal_respiratory_rate)  # add the loop to the experiment
    thisNormal_respiratory_rate = normal_respiratory_rate.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisNormal_respiratory_rate.rgb)
    if thisNormal_respiratory_rate != None:
        for paramName in thisNormal_respiratory_rate:
            globals()[paramName] = thisNormal_respiratory_rate[paramName]
    
    for thisNormal_respiratory_rate in normal_respiratory_rate:
        currentLoop = normal_respiratory_rate
        thisExp.timestampOnFlip(win, 'thisRow.t')
        # pause experiment here if requested
        if thisExp.status == PAUSED:
            pauseExperiment(
                thisExp=thisExp, 
                inputs=inputs, 
                win=win, 
                timers=[routineTimer], 
                playbackComponents=[]
        )
        # abbreviate parameter names if possible (e.g. rgb = thisNormal_respiratory_rate.rgb)
        if thisNormal_respiratory_rate != None:
            for paramName in thisNormal_respiratory_rate:
                globals()[paramName] = thisNormal_respiratory_rate[paramName]
        
        # --- Prepare to start Routine "breath_in" ---
        continueRoutine = True
        # update component parameters for each repeat
        thisExp.addData('breath_in.started', globalClock.getTime())
        # Run 'Begin Routine' code from channel_4
        send_message(signals.BHT_IN.to_bytes())
        # keep track of which components have finished
        breath_inComponents = [polygon_4]
        for thisComponent in breath_inComponents:
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
        
        # --- Run Routine "breath_in" ---
        routineForceEnded = not continueRoutine
        while continueRoutine and routineTimer.getTime() < 2.7:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *polygon_4* updates
            
            # if polygon_4 is starting this frame...
            if polygon_4.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                polygon_4.frameNStart = frameN  # exact frame index
                polygon_4.tStart = t  # local t and not account for scr refresh
                polygon_4.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(polygon_4, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'polygon_4.started')
                # update status
                polygon_4.status = STARTED
                polygon_4.setAutoDraw(True)
            
            # if polygon_4 is active this frame...
            if polygon_4.status == STARTED:
                # update params
                pass
            
            # if polygon_4 is stopping this frame...
            if polygon_4.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > polygon_4.tStartRefresh + 2.7-frameTolerance:
                    # keep track of stop time/frame for later
                    polygon_4.tStop = t  # not accounting for scr refresh
                    polygon_4.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'polygon_4.stopped')
                    # update status
                    polygon_4.status = FINISHED
                    polygon_4.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, inputs=inputs, win=win)
                return
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in breath_inComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "breath_in" ---
        for thisComponent in breath_inComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        thisExp.addData('breath_in.stopped', globalClock.getTime())
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-2.700000)
        
        # --- Prepare to start Routine "breath_out" ---
        continueRoutine = True
        # update component parameters for each repeat
        thisExp.addData('breath_out.started', globalClock.getTime())
        # Run 'Begin Routine' code from channel_7
        send_message(signals.BHT_OUT.to_bytes())
        # keep track of which components have finished
        breath_outComponents = [polygon1]
        for thisComponent in breath_outComponents:
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
        
        # --- Run Routine "breath_out" ---
        routineForceEnded = not continueRoutine
        while continueRoutine and routineTimer.getTime() < 2.3:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *polygon1* updates
            
            # if polygon1 is starting this frame...
            if polygon1.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                polygon1.frameNStart = frameN  # exact frame index
                polygon1.tStart = t  # local t and not account for scr refresh
                polygon1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(polygon1, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'polygon1.started')
                # update status
                polygon1.status = STARTED
                polygon1.setAutoDraw(True)
            
            # if polygon1 is active this frame...
            if polygon1.status == STARTED:
                # update params
                pass
            
            # if polygon1 is stopping this frame...
            if polygon1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > polygon1.tStartRefresh + 2.3-frameTolerance:
                    # keep track of stop time/frame for later
                    polygon1.tStop = t  # not accounting for scr refresh
                    polygon1.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'polygon1.stopped')
                    # update status
                    polygon1.status = FINISHED
                    polygon1.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, inputs=inputs, win=win)
                return
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in breath_outComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "breath_out" ---
        for thisComponent in breath_outComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        thisExp.addData('breath_out.stopped', globalClock.getTime())
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-2.300000)
        thisExp.nextEntry()
        
        if thisSession is not None:
            # if running in a Session with a Liaison client, send data up to now
            thisSession.sendExperimentData()
    # completed 4.0 repeats of 'normal_respiratory_rate'
    
    # get names of stimulus parameters
    if normal_respiratory_rate.trialList in ([], [None], None):
        params = []
    else:
        params = normal_respiratory_rate.trialList[0].keys()
    # save data for this loop
    normal_respiratory_rate.saveAsText(filename + 'normal_respiratory_rate.csv', delim=',',
        stimOut=params,
        dataOut=['n','all_mean','all_std', 'all_raw'])
    
    # --- Prepare to start Routine "breath_in_last" ---
    continueRoutine = True
    # update component parameters for each repeat
    thisExp.addData('breath_in_last.started', globalClock.getTime())
    # Run 'Begin Routine' code from channel_4_2_
    send_message(signals.BHT_IN.to_bytes())
    # keep track of which components have finished
    breath_in_lastComponents = [polygon_6]
    for thisComponent in breath_in_lastComponents:
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
    
    # --- Run Routine "breath_in_last" ---
    routineForceEnded = not continueRoutine
    while continueRoutine and routineTimer.getTime() < 2.7:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *polygon_6* updates
        
        # if polygon_6 is starting this frame...
        if polygon_6.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            polygon_6.frameNStart = frameN  # exact frame index
            polygon_6.tStart = t  # local t and not account for scr refresh
            polygon_6.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(polygon_6, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'polygon_6.started')
            # update status
            polygon_6.status = STARTED
            polygon_6.setAutoDraw(True)
        
        # if polygon_6 is active this frame...
        if polygon_6.status == STARTED:
            # update params
            pass
        
        # if polygon_6 is stopping this frame...
        if polygon_6.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > polygon_6.tStartRefresh + 2.7-frameTolerance:
                # keep track of stop time/frame for later
                polygon_6.tStop = t  # not accounting for scr refresh
                polygon_6.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'polygon_6.stopped')
                # update status
                polygon_6.status = FINISHED
                polygon_6.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, inputs=inputs, win=win)
            return
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in breath_in_lastComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "breath_in_last" ---
    for thisComponent in breath_in_lastComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    thisExp.addData('breath_in_last.stopped', globalClock.getTime())
    # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
    if routineForceEnded:
        routineTimer.reset()
    else:
        routineTimer.addTime(-2.700000)
    
    # --- Prepare to start Routine "breath_out_last" ---
    continueRoutine = True
    # update component parameters for each repeat
    thisExp.addData('breath_out_last.started', globalClock.getTime())
    # Run 'Begin Routine' code from channel
    send_message(signals.BHT_OUT.to_bytes())
    # keep track of which components have finished
    breath_out_lastComponents = [polygon_8]
    for thisComponent in breath_out_lastComponents:
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
    
    # --- Run Routine "breath_out_last" ---
    routineForceEnded = not continueRoutine
    while continueRoutine and routineTimer.getTime() < 2.3:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *polygon_8* updates
        
        # if polygon_8 is starting this frame...
        if polygon_8.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
            # keep track of start time/frame for later
            polygon_8.frameNStart = frameN  # exact frame index
            polygon_8.tStart = t  # local t and not account for scr refresh
            polygon_8.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(polygon_8, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'polygon_8.started')
            # update status
            polygon_8.status = STARTED
            polygon_8.setAutoDraw(True)
        
        # if polygon_8 is active this frame...
        if polygon_8.status == STARTED:
            # update params
            pass
        
        # if polygon_8 is stopping this frame...
        if polygon_8.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > polygon_8.tStartRefresh + 2.3-frameTolerance:
                # keep track of stop time/frame for later
                polygon_8.tStop = t  # not accounting for scr refresh
                polygon_8.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'polygon_8.stopped')
                # update status
                polygon_8.status = FINISHED
                polygon_8.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, inputs=inputs, win=win)
            return
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in breath_out_lastComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "breath_out_last" ---
    for thisComponent in breath_out_lastComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    thisExp.addData('breath_out_last.stopped', globalClock.getTime())
    # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
    if routineForceEnded:
        routineTimer.reset()
    else:
        routineTimer.addTime(-2.300000)
    
    # --- Prepare to start Routine "hold_test" ---
    continueRoutine = True
    # update component parameters for each repeat
    thisExp.addData('hold_test.started', globalClock.getTime())
    # Run 'Begin Routine' code from channel_8
    send_message(signals.BHT_HOLD_START_AND_STOP.to_bytes())
    
    # keep track of which components have finished
    hold_testComponents = [bh_body, bh_end, end_trial_msg]
    for thisComponent in hold_testComponents:
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
    
    # --- Run Routine "hold_test" ---
    routineForceEnded = not continueRoutine
    while continueRoutine and routineTimer.getTime() < 20.0:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *bh_body* updates
        
        # if bh_body is starting this frame...
        if bh_body.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            bh_body.frameNStart = frameN  # exact frame index
            bh_body.tStart = t  # local t and not account for scr refresh
            bh_body.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(bh_body, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'bh_body.started')
            # update status
            bh_body.status = STARTED
            bh_body.setAutoDraw(True)
        
        # if bh_body is active this frame...
        if bh_body.status == STARTED:
            # update params
            pass
        
        # if bh_body is stopping this frame...
        if bh_body.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > bh_body.tStartRefresh + 13-frameTolerance:
                # keep track of stop time/frame for later
                bh_body.tStop = t  # not accounting for scr refresh
                bh_body.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'bh_body.stopped')
                # update status
                bh_body.status = FINISHED
                bh_body.setAutoDraw(False)
        
        # *bh_end* updates
        
        # if bh_end is starting this frame...
        if bh_end.status == NOT_STARTED and tThisFlip >= 13-frameTolerance:
            # keep track of start time/frame for later
            bh_end.frameNStart = frameN  # exact frame index
            bh_end.tStart = t  # local t and not account for scr refresh
            bh_end.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(bh_end, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'bh_end.started')
            # update status
            bh_end.status = STARTED
            bh_end.setAutoDraw(True)
        
        # if bh_end is active this frame...
        if bh_end.status == STARTED:
            # update params
            pass
        
        # if bh_end is stopping this frame...
        if bh_end.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > bh_end.tStartRefresh + 2-frameTolerance:
                # keep track of stop time/frame for later
                bh_end.tStop = t  # not accounting for scr refresh
                bh_end.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'bh_end.stopped')
                # update status
                bh_end.status = FINISHED
                bh_end.setAutoDraw(False)
        
        # *end_trial_msg* updates
        
        # if end_trial_msg is starting this frame...
        if end_trial_msg.status == NOT_STARTED and tThisFlip >= 15-frameTolerance:
            # keep track of start time/frame for later
            end_trial_msg.frameNStart = frameN  # exact frame index
            end_trial_msg.tStart = t  # local t and not account for scr refresh
            end_trial_msg.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(end_trial_msg, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'end_trial_msg.started')
            # update status
            end_trial_msg.status = STARTED
            end_trial_msg.setAutoDraw(True)
        
        # if end_trial_msg is active this frame...
        if end_trial_msg.status == STARTED:
            # update params
            pass
        
        # if end_trial_msg is stopping this frame...
        if end_trial_msg.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > end_trial_msg.tStartRefresh + 5.0-frameTolerance:
                # keep track of stop time/frame for later
                end_trial_msg.tStop = t  # not accounting for scr refresh
                end_trial_msg.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'end_trial_msg.stopped')
                # update status
                end_trial_msg.status = FINISHED
                end_trial_msg.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, inputs=inputs, win=win)
            return
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in hold_testComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "hold_test" ---
    for thisComponent in hold_testComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    thisExp.addData('hold_test.stopped', globalClock.getTime())
    # Run 'End Routine' code from channel_8
    send_message(signals.BHT_HOLD_START_AND_STOP.to_bytes())
    ioServer.getDevice('tracker').sendMessage("stop mock block")
    ioServer.getDevice('tracker').sendMessage("start bh block")
    
    # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
    if routineForceEnded:
        routineTimer.reset()
    else:
        routineTimer.addTime(-20.000000)
    
    # set up handler to look after randomisation of conditions etc
    trials = data.TrialHandler(nReps=5.0, method='sequential', 
        extraInfo=expInfo, originPath=-1,
        trialList=[None],
        seed=None, name='trials')
    thisExp.addLoop(trials)  # add the loop to the experiment
    thisTrial = trials.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
    if thisTrial != None:
        for paramName in thisTrial:
            globals()[paramName] = thisTrial[paramName]
    
    for thisTrial in trials:
        currentLoop = trials
        thisExp.timestampOnFlip(win, 'thisRow.t')
        # pause experiment here if requested
        if thisExp.status == PAUSED:
            pauseExperiment(
                thisExp=thisExp, 
                inputs=inputs, 
                win=win, 
                timers=[routineTimer], 
                playbackComponents=[]
        )
        # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
        if thisTrial != None:
            for paramName in thisTrial:
                globals()[paramName] = thisTrial[paramName]
        
        # set up handler to look after randomisation of conditions etc
        normal_respiratory_rate_2 = data.TrialHandler(nReps=4.0, method='sequential', 
            extraInfo=expInfo, originPath=-1,
            trialList=[None],
            seed=None, name='normal_respiratory_rate_2')
        thisExp.addLoop(normal_respiratory_rate_2)  # add the loop to the experiment
        thisNormal_respiratory_rate_2 = normal_respiratory_rate_2.trialList[0]  # so we can initialise stimuli with some values
        # abbreviate parameter names if possible (e.g. rgb = thisNormal_respiratory_rate_2.rgb)
        if thisNormal_respiratory_rate_2 != None:
            for paramName in thisNormal_respiratory_rate_2:
                globals()[paramName] = thisNormal_respiratory_rate_2[paramName]
        
        for thisNormal_respiratory_rate_2 in normal_respiratory_rate_2:
            currentLoop = normal_respiratory_rate_2
            thisExp.timestampOnFlip(win, 'thisRow.t')
            # pause experiment here if requested
            if thisExp.status == PAUSED:
                pauseExperiment(
                    thisExp=thisExp, 
                    inputs=inputs, 
                    win=win, 
                    timers=[routineTimer], 
                    playbackComponents=[]
            )
            # abbreviate parameter names if possible (e.g. rgb = thisNormal_respiratory_rate_2.rgb)
            if thisNormal_respiratory_rate_2 != None:
                for paramName in thisNormal_respiratory_rate_2:
                    globals()[paramName] = thisNormal_respiratory_rate_2[paramName]
            
            # --- Prepare to start Routine "breath_in" ---
            continueRoutine = True
            # update component parameters for each repeat
            thisExp.addData('breath_in.started', globalClock.getTime())
            # Run 'Begin Routine' code from channel_4
            send_message(signals.BHT_IN.to_bytes())
            # keep track of which components have finished
            breath_inComponents = [polygon_4]
            for thisComponent in breath_inComponents:
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
            
            # --- Run Routine "breath_in" ---
            routineForceEnded = not continueRoutine
            while continueRoutine and routineTimer.getTime() < 2.7:
                # get current time
                t = routineTimer.getTime()
                tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                
                # *polygon_4* updates
                
                # if polygon_4 is starting this frame...
                if polygon_4.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    polygon_4.frameNStart = frameN  # exact frame index
                    polygon_4.tStart = t  # local t and not account for scr refresh
                    polygon_4.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(polygon_4, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'polygon_4.started')
                    # update status
                    polygon_4.status = STARTED
                    polygon_4.setAutoDraw(True)
                
                # if polygon_4 is active this frame...
                if polygon_4.status == STARTED:
                    # update params
                    pass
                
                # if polygon_4 is stopping this frame...
                if polygon_4.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > polygon_4.tStartRefresh + 2.7-frameTolerance:
                        # keep track of stop time/frame for later
                        polygon_4.tStop = t  # not accounting for scr refresh
                        polygon_4.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'polygon_4.stopped')
                        # update status
                        polygon_4.status = FINISHED
                        polygon_4.setAutoDraw(False)
                
                # check for quit (typically the Esc key)
                if defaultKeyboard.getKeys(keyList=["escape"]):
                    thisExp.status = FINISHED
                if thisExp.status == FINISHED or endExpNow:
                    endExperiment(thisExp, inputs=inputs, win=win)
                    return
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    routineForceEnded = True
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in breath_inComponents:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # --- Ending Routine "breath_in" ---
            for thisComponent in breath_inComponents:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            thisExp.addData('breath_in.stopped', globalClock.getTime())
            # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
            if routineForceEnded:
                routineTimer.reset()
            else:
                routineTimer.addTime(-2.700000)
            
            # --- Prepare to start Routine "breath_out" ---
            continueRoutine = True
            # update component parameters for each repeat
            thisExp.addData('breath_out.started', globalClock.getTime())
            # Run 'Begin Routine' code from channel_7
            send_message(signals.BHT_OUT.to_bytes())
            # keep track of which components have finished
            breath_outComponents = [polygon1]
            for thisComponent in breath_outComponents:
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
            
            # --- Run Routine "breath_out" ---
            routineForceEnded = not continueRoutine
            while continueRoutine and routineTimer.getTime() < 2.3:
                # get current time
                t = routineTimer.getTime()
                tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                
                # *polygon1* updates
                
                # if polygon1 is starting this frame...
                if polygon1.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    polygon1.frameNStart = frameN  # exact frame index
                    polygon1.tStart = t  # local t and not account for scr refresh
                    polygon1.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(polygon1, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'polygon1.started')
                    # update status
                    polygon1.status = STARTED
                    polygon1.setAutoDraw(True)
                
                # if polygon1 is active this frame...
                if polygon1.status == STARTED:
                    # update params
                    pass
                
                # if polygon1 is stopping this frame...
                if polygon1.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > polygon1.tStartRefresh + 2.3-frameTolerance:
                        # keep track of stop time/frame for later
                        polygon1.tStop = t  # not accounting for scr refresh
                        polygon1.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'polygon1.stopped')
                        # update status
                        polygon1.status = FINISHED
                        polygon1.setAutoDraw(False)
                
                # check for quit (typically the Esc key)
                if defaultKeyboard.getKeys(keyList=["escape"]):
                    thisExp.status = FINISHED
                if thisExp.status == FINISHED or endExpNow:
                    endExperiment(thisExp, inputs=inputs, win=win)
                    return
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    routineForceEnded = True
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in breath_outComponents:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # --- Ending Routine "breath_out" ---
            for thisComponent in breath_outComponents:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            thisExp.addData('breath_out.stopped', globalClock.getTime())
            # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
            if routineForceEnded:
                routineTimer.reset()
            else:
                routineTimer.addTime(-2.300000)
            thisExp.nextEntry()
            
            if thisSession is not None:
                # if running in a Session with a Liaison client, send data up to now
                thisSession.sendExperimentData()
        # completed 4.0 repeats of 'normal_respiratory_rate_2'
        
        # get names of stimulus parameters
        if normal_respiratory_rate_2.trialList in ([], [None], None):
            params = []
        else:
            params = normal_respiratory_rate_2.trialList[0].keys()
        # save data for this loop
        normal_respiratory_rate_2.saveAsText(filename + 'normal_respiratory_rate_2.csv', delim=',',
            stimOut=params,
            dataOut=['n','all_mean','all_std', 'all_raw'])
        
        # --- Prepare to start Routine "breath_in_last" ---
        continueRoutine = True
        # update component parameters for each repeat
        thisExp.addData('breath_in_last.started', globalClock.getTime())
        # Run 'Begin Routine' code from channel_4_2_
        send_message(signals.BHT_IN.to_bytes())
        # keep track of which components have finished
        breath_in_lastComponents = [polygon_6]
        for thisComponent in breath_in_lastComponents:
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
        
        # --- Run Routine "breath_in_last" ---
        routineForceEnded = not continueRoutine
        while continueRoutine and routineTimer.getTime() < 2.7:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *polygon_6* updates
            
            # if polygon_6 is starting this frame...
            if polygon_6.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                polygon_6.frameNStart = frameN  # exact frame index
                polygon_6.tStart = t  # local t and not account for scr refresh
                polygon_6.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(polygon_6, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'polygon_6.started')
                # update status
                polygon_6.status = STARTED
                polygon_6.setAutoDraw(True)
            
            # if polygon_6 is active this frame...
            if polygon_6.status == STARTED:
                # update params
                pass
            
            # if polygon_6 is stopping this frame...
            if polygon_6.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > polygon_6.tStartRefresh + 2.7-frameTolerance:
                    # keep track of stop time/frame for later
                    polygon_6.tStop = t  # not accounting for scr refresh
                    polygon_6.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'polygon_6.stopped')
                    # update status
                    polygon_6.status = FINISHED
                    polygon_6.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, inputs=inputs, win=win)
                return
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in breath_in_lastComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "breath_in_last" ---
        for thisComponent in breath_in_lastComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        thisExp.addData('breath_in_last.stopped', globalClock.getTime())
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-2.700000)
        
        # --- Prepare to start Routine "breath_out_last" ---
        continueRoutine = True
        # update component parameters for each repeat
        thisExp.addData('breath_out_last.started', globalClock.getTime())
        # Run 'Begin Routine' code from channel
        send_message(signals.BHT_OUT.to_bytes())
        # keep track of which components have finished
        breath_out_lastComponents = [polygon_8]
        for thisComponent in breath_out_lastComponents:
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
        
        # --- Run Routine "breath_out_last" ---
        routineForceEnded = not continueRoutine
        while continueRoutine and routineTimer.getTime() < 2.3:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *polygon_8* updates
            
            # if polygon_8 is starting this frame...
            if polygon_8.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                polygon_8.frameNStart = frameN  # exact frame index
                polygon_8.tStart = t  # local t and not account for scr refresh
                polygon_8.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(polygon_8, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'polygon_8.started')
                # update status
                polygon_8.status = STARTED
                polygon_8.setAutoDraw(True)
            
            # if polygon_8 is active this frame...
            if polygon_8.status == STARTED:
                # update params
                pass
            
            # if polygon_8 is stopping this frame...
            if polygon_8.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > polygon_8.tStartRefresh + 2.3-frameTolerance:
                    # keep track of stop time/frame for later
                    polygon_8.tStop = t  # not accounting for scr refresh
                    polygon_8.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'polygon_8.stopped')
                    # update status
                    polygon_8.status = FINISHED
                    polygon_8.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, inputs=inputs, win=win)
                return
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in breath_out_lastComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "breath_out_last" ---
        for thisComponent in breath_out_lastComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        thisExp.addData('breath_out_last.stopped', globalClock.getTime())
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-2.300000)
        
        # --- Prepare to start Routine "hold" ---
        continueRoutine = True
        # update component parameters for each repeat
        thisExp.addData('hold.started', globalClock.getTime())
        # Run 'Begin Routine' code from channel_8_2_
        send_message(signals.BHT_HOLD_START_AND_STOP.to_bytes())
        
        # keep track of which components have finished
        holdComponents = [bh_body_2, bh_end_2, bh_end_3]
        for thisComponent in holdComponents:
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
        
        # --- Run Routine "hold" ---
        routineForceEnded = not continueRoutine
        while continueRoutine and routineTimer.getTime() < 25.0:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *bh_body_2* updates
            
            # if bh_body_2 is starting this frame...
            if bh_body_2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                bh_body_2.frameNStart = frameN  # exact frame index
                bh_body_2.tStart = t  # local t and not account for scr refresh
                bh_body_2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(bh_body_2, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'bh_body_2.started')
                # update status
                bh_body_2.status = STARTED
                bh_body_2.setAutoDraw(True)
            
            # if bh_body_2 is active this frame...
            if bh_body_2.status == STARTED:
                # update params
                pass
            
            # if bh_body_2 is stopping this frame...
            if bh_body_2.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > bh_body_2.tStartRefresh + 13-frameTolerance:
                    # keep track of stop time/frame for later
                    bh_body_2.tStop = t  # not accounting for scr refresh
                    bh_body_2.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'bh_body_2.stopped')
                    # update status
                    bh_body_2.status = FINISHED
                    bh_body_2.setAutoDraw(False)
            
            # *bh_end_2* updates
            
            # if bh_end_2 is starting this frame...
            if bh_end_2.status == NOT_STARTED and tThisFlip >= 13-frameTolerance:
                # keep track of start time/frame for later
                bh_end_2.frameNStart = frameN  # exact frame index
                bh_end_2.tStart = t  # local t and not account for scr refresh
                bh_end_2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(bh_end_2, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'bh_end_2.started')
                # update status
                bh_end_2.status = STARTED
                bh_end_2.setAutoDraw(True)
            
            # if bh_end_2 is active this frame...
            if bh_end_2.status == STARTED:
                # update params
                pass
            
            # if bh_end_2 is stopping this frame...
            if bh_end_2.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > bh_end_2.tStartRefresh + 2-frameTolerance:
                    # keep track of stop time/frame for later
                    bh_end_2.tStop = t  # not accounting for scr refresh
                    bh_end_2.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'bh_end_2.stopped')
                    # update status
                    bh_end_2.status = FINISHED
                    bh_end_2.setAutoDraw(False)
            
            # *bh_end_3* updates
            
            # if bh_end_3 is starting this frame...
            if bh_end_3.status == NOT_STARTED and tThisFlip >= 15-frameTolerance:
                # keep track of start time/frame for later
                bh_end_3.frameNStart = frameN  # exact frame index
                bh_end_3.tStart = t  # local t and not account for scr refresh
                bh_end_3.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(bh_end_3, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'bh_end_3.started')
                # update status
                bh_end_3.status = STARTED
                bh_end_3.setAutoDraw(True)
            
            # if bh_end_3 is active this frame...
            if bh_end_3.status == STARTED:
                # update params
                pass
            
            # if bh_end_3 is stopping this frame...
            if bh_end_3.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > bh_end_3.tStartRefresh + 10-frameTolerance:
                    # keep track of stop time/frame for later
                    bh_end_3.tStop = t  # not accounting for scr refresh
                    bh_end_3.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'bh_end_3.stopped')
                    # update status
                    bh_end_3.status = FINISHED
                    bh_end_3.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, inputs=inputs, win=win)
                return
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in holdComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "hold" ---
        for thisComponent in holdComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        thisExp.addData('hold.stopped', globalClock.getTime())
        # Run 'End Routine' code from channel_8_2_
        send_message(signals.BHT_HOLD_START_AND_STOP.to_bytes())
        
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-25.000000)
        thisExp.nextEntry()
        
        if thisSession is not None:
            # if running in a Session with a Liaison client, send data up to now
            thisSession.sendExperimentData()
    # completed 5.0 repeats of 'trials'
    
    # get names of stimulus parameters
    if trials.trialList in ([], [None], None):
        params = []
    else:
        params = trials.trialList[0].keys()
    # save data for this loop
    trials.saveAsText(filename + 'trials.csv', delim=',',
        stimOut=params,
        dataOut=['n','all_mean','all_std', 'all_raw'])
    
    # --- Prepare to start Routine "done" ---
    continueRoutine = True
    # update component parameters for each repeat
    thisExp.addData('done.started', globalClock.getTime())
    # Run 'Begin Routine' code from code_6
    ioServer.getDevice('tracker').sendMessage("stop bh block")
    # keep track of which components have finished
    doneComponents = [finished]
    for thisComponent in doneComponents:
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
    
    # --- Run Routine "done" ---
    routineForceEnded = not continueRoutine
    while continueRoutine and routineTimer.getTime() < 25.0:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *finished* updates
        
        # if finished is starting this frame...
        if finished.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            finished.frameNStart = frameN  # exact frame index
            finished.tStart = t  # local t and not account for scr refresh
            finished.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(finished, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'finished.started')
            # update status
            finished.status = STARTED
            finished.setAutoDraw(True)
        
        # if finished is active this frame...
        if finished.status == STARTED:
            # update params
            pass
        
        # if finished is stopping this frame...
        if finished.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > finished.tStartRefresh + 25.0-frameTolerance:
                # keep track of stop time/frame for later
                finished.tStop = t  # not accounting for scr refresh
                finished.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'finished.stopped')
                # update status
                finished.status = FINISHED
                finished.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, inputs=inputs, win=win)
            return
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in doneComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "done" ---
    for thisComponent in doneComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    thisExp.addData('done.stopped', globalClock.getTime())
    # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
    if routineForceEnded:
        routineTimer.reset()
    else:
        routineTimer.addTime(-25.000000)
    
    # --- Prepare to start Routine "stop_et" ---
    continueRoutine = True
    # update component parameters for each repeat
    thisExp.addData('stop_et.started', globalClock.getTime())
    # Run 'Begin Routine' code from channel_2_2
    send_message(signals.ET_START_AND_STOP.to_bytes())
    
    ioServer.getDevice('tracker').sendMessage("bye bht")
    # keep track of which components have finished
    stop_etComponents = [etRecord_2]
    for thisComponent in stop_etComponents:
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
    
    # --- Run Routine "stop_et" ---
    routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        # *etRecord_2* updates
        
        # if etRecord_2 is stopping this frame...
        if etRecord_2.status == STARTED:
            if bool(True):
                # keep track of stop time/frame for later
                etRecord_2.tStop = t  # not accounting for scr refresh
                etRecord_2.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.addData('etRecord_2.stopped', t)
                # update status
                etRecord_2.status = FINISHED
        # Run 'Each Frame' code from channel_2_2
        if etRecord_2.status == NOT_STARTED:
            etRecord_2.status = STARTED
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, inputs=inputs, win=win)
            return
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in stop_etComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "stop_et" ---
    for thisComponent in stop_etComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    thisExp.addData('stop_et.stopped', globalClock.getTime())
    # make sure the eyetracker recording stops
    if etRecord_2.status != FINISHED:
        etRecord_2.status = FINISHED
    # the Routine "stop_et" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # mark experiment as finished
    endExperiment(thisExp, win=win, inputs=inputs)


def saveData(thisExp):
    """
    Save data from this experiment
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    """
    filename = thisExp.dataFileName
    # these shouldn't be strictly necessary (should auto-save)
    thisExp.saveAsWideText(filename + '.csv', delim='auto')
    thisExp.saveAsPickle(filename)


def endExperiment(thisExp, inputs=None, win=None):
    """
    End this experiment, performing final shut down operations.
    
    This function does NOT close the window or end the Python process - use `quit` for this.
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    inputs : dict
        Dictionary of input devices by name.
    win : psychopy.visual.Window
        Window for this experiment.
    """
    if win is not None:
        # remove autodraw from all current components
        win.clearAutoDraw()
        # Flip one final time so any remaining win.callOnFlip() 
        # and win.timeOnFlip() tasks get executed
        win.flip()
    # mark experiment handler as finished
    thisExp.status = FINISHED
    # shut down eyetracker, if there is one
    if inputs is not None:
        if 'eyetracker' in inputs and inputs['eyetracker'] is not None:
            inputs['eyetracker'].setConnectionState(False)
    logging.flush()


def quit(thisExp, win=None, inputs=None, thisSession=None):
    """
    Fully quit, closing the window and ending the Python process.
    
    Parameters
    ==========
    win : psychopy.visual.Window
        Window to close.
    inputs : dict
        Dictionary of input devices by name.
    thisSession : psychopy.session.Session or None
        Handle of the Session object this experiment is being run from, if any.
    """
    thisExp.abort()  # or data files will save again on exit
    # make sure everything is closed down
    if win is not None:
        # Flip one final time so any remaining win.callOnFlip() 
        # and win.timeOnFlip() tasks get executed before quitting
        win.flip()
        win.close()
    if inputs is not None:
        if 'eyetracker' in inputs and inputs['eyetracker'] is not None:
            inputs['eyetracker'].setConnectionState(False)
    logging.flush()
    if thisSession is not None:
        thisSession.stop()
    # terminate Python process
    core.quit()


# if running this experiment as a script...
if __name__ == '__main__':
    # call all functions in order
    expInfo = showExpInfoDlg(expInfo=expInfo)
    thisExp = setupData(expInfo=expInfo)
    logFile = setupLogging(filename=thisExp.dataFileName)
    win = setupWindow(expInfo=expInfo)
    inputs = setupInputs(expInfo=expInfo, thisExp=thisExp, win=win)
    run(
        expInfo=expInfo, 
        thisExp=thisExp, 
        win=win, 
        inputs=inputs
    )
    saveData(thisExp=thisExp)
    quit(thisExp=thisExp, win=win, inputs=inputs)
