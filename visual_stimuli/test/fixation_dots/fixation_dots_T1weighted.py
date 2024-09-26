#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2024.2.1),
    on September 19, 2024, at 15:08
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
prefs.hardware['audioLatencyMode'] = '3'
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors, layout, hardware, iohub
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

# --- Setup global variables (available in all functions) ---
# create a device manager to handle hardware (keyboards, mice, mirophones, speakers, etc.)
deviceManager = hardware.DeviceManager()
# ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
# store info about the experiment session
psychopyVersion = '2024.2.1'
expName = 'fixation_dots_T1weighted'  # from the Builder filename that created this script
# information about this experiment
expInfo = {
    'participant': f"{randint(0, 999999):06.0f}",
    'session': '001',
    'date|hid': data.getDateStr(),
    'expName|hid': expName,
    'psychopyVersion|hid': psychopyVersion,
}

# --- Define some variables which will change depending on pilot mode ---
'''
To run in pilot mode, either use the run/pilot toggle in Builder, Coder and Runner, 
or run the experiment with `--pilot` as an argument. To change what pilot 
#mode does, check out the 'Pilot mode' tab in preferences.
'''
# work out from system args whether we are running in pilot mode
PILOTING = core.setPilotModeFromArgs()
# start off with values from experiment settings
_fullScr = True
_winSize = [800, 600]
# if in pilot mode, apply overrides according to preferences
if PILOTING:
    # force windowed mode
    if prefs.piloting['forceWindowed']:
        _fullScr = False
        # set window size
        _winSize = prefs.piloting['forcedWindowSize']

def showExpInfoDlg(expInfo):
    """
    Show participant info dialog.
    Parameters
    ==========
    expInfo : dict
        Information about this experiment.
    
    Returns
    ==========
    dict
        Information about this experiment.
    """
    # show participant info dialog
    dlg = gui.DlgFromDict(
        dictionary=expInfo, sortKeys=False, title=expName, alwaysOnTop=True
    )
    if dlg.OK == False:
        core.quit()  # user pressed cancel
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
    # remove dialog-specific syntax from expInfo
    for key, val in expInfo.copy().items():
        newKey, _ = data.utils.parsePipeSyntax(key)
        expInfo[newKey] = expInfo.pop(key)
    
    # data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
    if dataDir is None:
        dataDir = _thisDir
    filename = u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])
    # make sure filename is relative to dataDir
    if os.path.isabs(filename):
        dataDir = os.path.commonprefix([dataDir, filename])
        filename = os.path.relpath(filename, dataDir)
    
    # an ExperimentHandler isn't essential but helps with data saving
    thisExp = data.ExperimentHandler(
        name=expName, version='',
        extraInfo=expInfo, runtimeInfo=None,
        originPath='C:\\yiwei\\1_Pilot_MREye_Data\\1_anatomical-Protocol\\Code\\fixation_dots_T1weighted.py',
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
    # set how much information should be printed to the console / app
    if PILOTING:
        logging.console.setLevel(
            prefs.piloting['pilotConsoleLoggingLevel']
        )
    else:
        logging.console.setLevel('warning')
    # save a log file for detail verbose info
    logFile = logging.LogFile(filename+'.log')
    if PILOTING:
        logFile.setLevel(
            prefs.piloting['pilotLoggingLevel']
        )
    else:
        logFile.setLevel(
            logging.getLevel('exp')
        )
    
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
    if PILOTING:
        logging.debug('Fullscreen settings ignored as running in pilot mode.')
    
    if win is None:
        # if not given a window to setup, make one
        win = visual.Window(
            size=_winSize, fullscr=_fullScr, screen=0,
            winType='pyglet', allowStencil=False,
            monitor='testMonitor', color=[0,0,0], colorSpace='rgb',
            backgroundImage='', backgroundFit='none',
            blendMode='avg', useFBO=True,
            units='norm', 
            checkTiming=False  # we're going to do this ourselves in a moment
        )
    else:
        # if we have a window, just set the attributes which are safe to set
        win.color = [0,0,0]
        win.colorSpace = 'rgb'
        win.backgroundImage = ''
        win.backgroundFit = 'none'
        win.units = 'norm'
    if expInfo is not None:
        # get/measure frame rate if not already in expInfo
        if win._monitorFrameRate is None:
            win._monitorFrameRate = win.getActualFrameRate(infoMsg='Attempting to measure frame rate of screen, please wait...')
        expInfo['frameRate'] = win._monitorFrameRate
    win.mouseVisible = False
    win.hideMessage()
    # show a visual indicator if we're in piloting mode
    if PILOTING and prefs.piloting['showPilotingIndicator']:
        win.showPilotingIndicator()
    
    return win


def setupDevices(expInfo, thisExp, win):
    """
    Setup whatever devices are available (mouse, keyboard, speaker, eyetracker, etc.) and add them to 
    the device manager (deviceManager)
    
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
    bool
        True if completed successfully.
    """
    # --- Setup input devices ---
    ioConfig = {}
    
    # Setup eyetracking
    ioConfig['eyetracker.eyelink.EyeTracker'] = {
        'name': 'tracker',
        'model_name': 'EYELINK 1000 DESKTOP',
        'simulation_mode': False,
        'network_settings': '100.1.1.1',
        'default_native_data_file_name': 'EXPFILE',
        'runtime_settings': {
            'sampling_rate': 1000.0,
            'track_eyes': 'RIGHT_EYE',
            'sample_filtering': {
                'FILTER_FILE': 'FILTER_LEVEL_2',
                'FILTER_ONLINE': 'FILTER_LEVEL_OFF',
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
    
    # Setup iohub experiment
    ioConfig['Experiment'] = dict(filename=thisExp.dataFileName)
    
    # Start ioHub server
    ioServer = io.launchHubServer(window=win, **ioConfig)
    
    # store ioServer object in the device manager
    deviceManager.ioServer = ioServer
    deviceManager.devices['eyetracker'] = ioServer.getDevice('tracker')
    
    # create a default keyboard (e.g. to check for escape)
    if deviceManager.getDevice('defaultKeyboard') is None:
        deviceManager.addDevice(
            deviceClass='keyboard', deviceName='defaultKeyboard', backend='iohub'
        )
    if deviceManager.getDevice('key_resp') is None:
        # initialise key_resp
        key_resp = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='key_resp',
        )
    if deviceManager.getDevice('key_resp_3') is None:
        # initialise key_resp_3
        key_resp_3 = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='key_resp_3',
        )
    # return True if completed successfully
    return True

def pauseExperiment(thisExp, win=None, timers=[], playbackComponents=[]):
    """
    Pause this experiment, preventing the flow from advancing to the next routine until resumed.
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
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
    
    # start a timer to figure out how long we're paused for
    pauseTimer = core.Clock()
    # pause any playback components
    for comp in playbackComponents:
        comp.pause()
    # make sure we have a keyboard
    defaultKeyboard = deviceManager.getDevice('defaultKeyboard')
    if defaultKeyboard is None:
        defaultKeyboard = deviceManager.addKeyboard(
            deviceClass='keyboard',
            deviceName='defaultKeyboard',
            backend='ioHub',
        )
    # run a while loop while we wait to unpause
    while thisExp.status == PAUSED:
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=['escape']):
            endExperiment(thisExp, win=win)
        # sleep 1ms so other threads can execute
        clock.time.sleep(0.001)
    # if stop was requested while paused, quit
    if thisExp.status == FINISHED:
        endExperiment(thisExp, win=win)
    # resume any playback components
    for comp in playbackComponents:
        comp.play()
    # reset any timers
    for timer in timers:
        timer.addTime(-pauseTimer.getTime())


def run(expInfo, thisExp, win, globalClock=None, thisSession=None):
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
    ioServer = deviceManager.ioServer
    # get/create a default keyboard (e.g. to check for escape)
    defaultKeyboard = deviceManager.getDevice('defaultKeyboard')
    if defaultKeyboard is None:
        deviceManager.addDevice(
            deviceClass='keyboard', deviceName='defaultKeyboard', backend='ioHub'
        )
    eyetracker = deviceManager.getDevice('eyetracker')
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
    
    # --- Initialize components for Routine "trail" ---
    waiting_trigger = visual.TextStim(win=win, name='waiting_trigger',
        text="The program is ready for the scanner trigger. Press 's' to proceed manually.",
        font='Open Sans',
        pos=(0, -0.4), draggable=False, height=0.12, wrapWidth=1.7, ori=0.0, 
        color='white', colorSpace='rgb', opacity=1.0, 
        languageStyle='LTR',
        depth=0.0);
    key_resp = keyboard.Keyboard(deviceName='key_resp')
    fix_desc = visual.TextStim(win=win, name='fix_desc',
        text='In this task you will see a grating patterns on the screen. \nPlease keep your eyes on the fixation circle like the one below.',
        font='Open Sans',
        pos=(0, 0.25), draggable=False, height=0.12, wrapWidth=1.0, ori=0.0, 
        color='white', colorSpace='rgb', opacity=1.0, 
        languageStyle='LTR',
        depth=-2.0);
    
    # --- Initialize components for Routine "start_ET" ---
    etRecord = hardware.eyetracker.EyetrackerControl(
        tracker=eyetracker,
        actionType='Start Only'
    )
    
    # --- Initialize components for Routine "dots" ---
    dot_1 = visual.ShapeStim(
        win=win, name='dot_1',
        size=(0.04, 0.04), vertices='circle',
        ori=0.0, pos=(0, 0), draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1, -0.2, -0.2], fillColor=[1.0000, -1.0000, -1.0000],
        opacity=None, depth=0.0, interpolate=True)
    dot_2 = visual.ShapeStim(
        win=win, name='dot_2',
        size=(0.04, 0.04), vertices='circle',
        ori=0.0, pos=(0, 0), draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor='white', fillColor=[-1.0000, 0.2706, -0.9843],
        opacity=None, depth=-1.0, interpolate=True)
    dot_3 = visual.ShapeStim(
        win=win, name='dot_3',
        size=(0.04, 0.04), vertices='circle',
        ori=0.0, pos=(0, 0), draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1,0,1], fillColor=[1.0000, 0.0745, -0.4667],
        opacity=None, depth=-2.0, interpolate=True)
    dot_4 = visual.ShapeStim(
        win=win, name='dot_4',
        size=(0.04, 0.04), vertices='circle',
        ori=0.0, pos=(0, 0), draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[-1,1,1], fillColor=[-1,1,1],
        opacity=None, depth=-3.0, interpolate=True)
    dot_5 = visual.ShapeStim(
        win=win, name='dot_5',
        size=(0.04, 0.04), vertices='circle',
        ori=0.0, pos=(0, 0), draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1,1,1], fillColor=[0.4510, 0.4275, 0.4353],
        opacity=None, depth=-4.0, interpolate=True)
    dot_6 = visual.ShapeStim(
        win=win, name='dot_6',
        size=(0.04, 0.04), vertices='circle',
        ori=0.0, pos=(0, 0), draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1,1,0.5], fillColor=[1.0000, 0.7647, -0.6235],
        opacity=None, depth=-5.0, interpolate=True)
    dot_7 = visual.ShapeStim(
        win=win, name='dot_7',
        size=(0.04, 0.04), vertices='circle',
        ori=0.0, pos=(0, 0), draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[-0.15, -0.1, 1], fillColor=[0.2235, -1.0000, 1.0000],
        opacity=None, depth=-6.0, interpolate=True)
    dot_8 = visual.ShapeStim(
        win=win, name='dot_8',
        size=(0.04, 0.04), vertices='circle',
        ori=0.0, pos=(0, 0), draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[0.5, 0.2, 0.2], fillColor=[0.3412, -0.1294, -0.5059],
        opacity=None, depth=-7.0, interpolate=True)
    dot_9 = visual.ShapeStim(
        win=win, name='dot_9',
        size=(0.04, 0.04), vertices='circle',
        ori=0.0, pos=(0, 0), draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[0.5, 1, 1], fillColor=[-0.4118, 0.5765, 0.8353],
        opacity=None, depth=-8.0, interpolate=True)
    dot_10 = visual.ShapeStim(
        win=win, name='dot_10',
        size=(0.04, 0.04), vertices='circle',
        ori=0.0, pos=(0, 0), draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[0,0,0], fillColor=[0,0,0],
        opacity=None, depth=-9.0, interpolate=True)
    
    # --- Initialize components for Routine "end" ---
    text = visual.TextStim(win=win, name='text',
        text="End press 't'",
        font='Open Sans',
        pos=(0, 0), draggable=False, height=0.12, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=0.0);
    ET_stop = hardware.eyetracker.EyetrackerControl(
        tracker=eyetracker,
        actionType='Stop Only'
    )
    key_resp_3 = keyboard.Keyboard(deviceName='key_resp_3')
    
    # create some handy timers
    
    # global clock to track the time since experiment started
    if globalClock is None:
        # create a clock if not given one
        globalClock = core.Clock()
    if isinstance(globalClock, str):
        # if given a string, make a clock accoridng to it
        if globalClock == 'float':
            # get timestamps as a simple value
            globalClock = core.Clock(format='float')
        elif globalClock == 'iso':
            # get timestamps in ISO format
            globalClock = core.Clock(format='%Y-%m-%d_%H:%M:%S.%f%z')
        else:
            # get timestamps in a custom format
            globalClock = core.Clock(format=globalClock)
    if ioServer is not None:
        ioServer.syncClock(globalClock)
    logging.setDefaultClock(globalClock)
    # routine timer to track time remaining of each (possibly non-slip) routine
    routineTimer = core.Clock()
    win.flip()  # flip window to reset last flip timer
    # store the exact time the global clock started
    expInfo['expStart'] = data.getDateStr(
        format='%Y-%m-%d %Hh%M.%S.%f %z', fractionalSecondDigits=6
    )
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
    calibration_2.run()
    # clear any keypresses from during calibration_2 so they don't interfere with the experiment
    defaultKeyboard.clearEvents()
    thisExp.nextEntry()
    # the Routine "calibration_2" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # --- Prepare to start Routine "trail" ---
    # create an object to store info about Routine trail
    trail = data.Routine(
        name='trail',
        components=[waiting_trigger, key_resp, fix_desc],
    )
    trail.status = NOT_STARTED
    continueRoutine = True
    # update component parameters for each repeat
    # create starting attributes for key_resp
    key_resp.keys = []
    key_resp.rt = []
    _key_resp_allKeys = []
    # store start times for trail
    trail.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
    trail.tStart = globalClock.getTime(format='float')
    trail.status = STARTED
    thisExp.addData('trail.started', trail.tStart)
    trail.maxDuration = None
    # keep track of which components have finished
    trailComponents = trail.components
    for thisComponent in trail.components:
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
    trail.forceEnded = routineForceEnded = not continueRoutine
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
            theseKeys = key_resp.getKeys(keyList=['s'], ignoreKeys=["escape"], waitRelease=False)
            _key_resp_allKeys.extend(theseKeys)
            if len(_key_resp_allKeys):
                key_resp.keys = _key_resp_allKeys[-1].name  # just the last key pressed
                key_resp.rt = _key_resp_allKeys[-1].rt
                key_resp.duration = _key_resp_allKeys[-1].duration
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
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, win=win)
            return
        # pause experiment here if requested
        if thisExp.status == PAUSED:
            pauseExperiment(
                thisExp=thisExp, 
                win=win, 
                timers=[routineTimer], 
                playbackComponents=[]
            )
            # skip the frame we paused on
            continue
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            trail.forceEnded = routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in trail.components:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "trail" ---
    for thisComponent in trail.components:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # store stop times for trail
    trail.tStop = globalClock.getTime(format='float')
    trail.tStopRefresh = tThisFlipGlobal
    thisExp.addData('trail.stopped', trail.tStop)
    # check responses
    if key_resp.keys in ['', [], None]:  # No response was made
        key_resp.keys = None
    thisExp.addData('key_resp.keys',key_resp.keys)
    if key_resp.keys != None:  # we had a response
        thisExp.addData('key_resp.rt', key_resp.rt)
        thisExp.addData('key_resp.duration', key_resp.duration)
    thisExp.nextEntry()
    # the Routine "trail" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # --- Prepare to start Routine "start_ET" ---
    # create an object to store info about Routine start_ET
    start_ET = data.Routine(
        name='start_ET',
        components=[etRecord],
    )
    start_ET.status = NOT_STARTED
    continueRoutine = True
    # update component parameters for each repeat
    # store start times for start_ET
    start_ET.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
    start_ET.tStart = globalClock.getTime(format='float')
    start_ET.status = STARTED
    thisExp.addData('start_ET.started', start_ET.tStart)
    start_ET.maxDuration = None
    # keep track of which components have finished
    start_ETComponents = start_ET.components
    for thisComponent in start_ET.components:
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
    start_ET.forceEnded = routineForceEnded = not continueRoutine
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
            etRecord.start()
        if etRecord.status == STARTED:
            etRecord.tStop = t  # not accounting for scr refresh
            etRecord.tStopRefresh = tThisFlipGlobal  # on global time
            etRecord.frameNStop = frameN  # exact frame index
            etRecord.status = FINISHED
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, win=win)
            return
        # pause experiment here if requested
        if thisExp.status == PAUSED:
            pauseExperiment(
                thisExp=thisExp, 
                win=win, 
                timers=[routineTimer], 
                playbackComponents=[]
            )
            # skip the frame we paused on
            continue
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            start_ET.forceEnded = routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in start_ET.components:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "start_ET" ---
    for thisComponent in start_ET.components:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # store stop times for start_ET
    start_ET.tStop = globalClock.getTime(format='float')
    start_ET.tStopRefresh = tThisFlipGlobal
    thisExp.addData('start_ET.stopped', start_ET.tStop)
    thisExp.nextEntry()
    # the Routine "start_ET" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    T1w_LIBRE = data.TrialHandler2(
        name='T1w_LIBRE',
        nReps=720.0, 
        method='random', 
        extraInfo=expInfo, 
        originPath=-1, 
        trialList=[None], 
        seed=None, 
    )
    thisExp.addLoop(T1w_LIBRE)  # add the loop to the experiment
    thisT1w_LIBRE = T1w_LIBRE.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisT1w_LIBRE.rgb)
    if thisT1w_LIBRE != None:
        for paramName in thisT1w_LIBRE:
            globals()[paramName] = thisT1w_LIBRE[paramName]
    if thisSession is not None:
        # if running in a Session with a Liaison client, send data up to now
        thisSession.sendExperimentData()
    
    for thisT1w_LIBRE in T1w_LIBRE:
        currentLoop = T1w_LIBRE
        thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
        if thisSession is not None:
            # if running in a Session with a Liaison client, send data up to now
            thisSession.sendExperimentData()
        # abbreviate parameter names if possible (e.g. rgb = thisT1w_LIBRE.rgb)
        if thisT1w_LIBRE != None:
            for paramName in thisT1w_LIBRE:
                globals()[paramName] = thisT1w_LIBRE[paramName]
        
        # --- Prepare to start Routine "dots" ---
        # create an object to store info about Routine dots
        dots = data.Routine(
            name='dots',
            components=[dot_1, dot_2, dot_3, dot_4, dot_5, dot_6, dot_7, dot_8, dot_9, dot_10],
        )
        dots.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        dot_2.setLineColor([0,1,0])
        # store start times for dots
        dots.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        dots.tStart = globalClock.getTime(format='float')
        dots.status = STARTED
        thisExp.addData('dots.started', dots.tStart)
        dots.maxDuration = None
        # keep track of which components have finished
        dotsComponents = dots.components
        for thisComponent in dots.components:
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
        # if trial has changed, end Routine now
        if isinstance(T1w_LIBRE, data.TrialHandler2) and thisT1w_LIBRE.thisN != T1w_LIBRE.thisTrial.thisN:
            continueRoutine = False
        dots.forceEnded = routineForceEnded = not continueRoutine
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
                    dot_1.tStopRefresh = tThisFlipGlobal  # on global time
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
                    dot_2.tStopRefresh = tThisFlipGlobal  # on global time
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
                    dot_3.tStopRefresh = tThisFlipGlobal  # on global time
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
                    dot_4.tStopRefresh = tThisFlipGlobal  # on global time
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
                    dot_5.tStopRefresh = tThisFlipGlobal  # on global time
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
                    dot_6.tStopRefresh = tThisFlipGlobal  # on global time
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
                    dot_7.tStopRefresh = tThisFlipGlobal  # on global time
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
                    dot_8.tStopRefresh = tThisFlipGlobal  # on global time
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
                    dot_9.tStopRefresh = tThisFlipGlobal  # on global time
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
                    dot_10.tStopRefresh = tThisFlipGlobal  # on global time
                    dot_10.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'dot_10.stopped')
                    # update status
                    dot_10.status = FINISHED
                    dot_10.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, win=win)
                return
            # pause experiment here if requested
            if thisExp.status == PAUSED:
                pauseExperiment(
                    thisExp=thisExp, 
                    win=win, 
                    timers=[routineTimer], 
                    playbackComponents=[]
                )
                # skip the frame we paused on
                continue
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                dots.forceEnded = routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in dots.components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "dots" ---
        for thisComponent in dots.components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # store stop times for dots
        dots.tStop = globalClock.getTime(format='float')
        dots.tStopRefresh = tThisFlipGlobal
        thisExp.addData('dots.stopped', dots.tStop)
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if dots.maxDurationReached:
            routineTimer.addTime(-dots.maxDuration)
        elif dots.forceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-1.000000)
        thisExp.nextEntry()
        
    # completed 720.0 repeats of 'T1w_LIBRE'
    
    if thisSession is not None:
        # if running in a Session with a Liaison client, send data up to now
        thisSession.sendExperimentData()
    
    # --- Prepare to start Routine "end" ---
    # create an object to store info about Routine end
    end = data.Routine(
        name='end',
        components=[text, ET_stop, key_resp_3],
    )
    end.status = NOT_STARTED
    continueRoutine = True
    # update component parameters for each repeat
    # create starting attributes for key_resp_3
    key_resp_3.keys = []
    key_resp_3.rt = []
    _key_resp_3_allKeys = []
    # store start times for end
    end.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
    end.tStart = globalClock.getTime(format='float')
    end.status = STARTED
    thisExp.addData('end.started', end.tStart)
    end.maxDuration = None
    # keep track of which components have finished
    endComponents = end.components
    for thisComponent in end.components:
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
    end.forceEnded = routineForceEnded = not continueRoutine
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
            theseKeys = key_resp_3.getKeys(keyList=['t'], ignoreKeys=["escape"], waitRelease=False)
            _key_resp_3_allKeys.extend(theseKeys)
            if len(_key_resp_3_allKeys):
                key_resp_3.keys = _key_resp_3_allKeys[-1].name  # just the last key pressed
                key_resp_3.rt = _key_resp_3_allKeys[-1].rt
                key_resp_3.duration = _key_resp_3_allKeys[-1].duration
                # a response ends the routine
                continueRoutine = False
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, win=win)
            return
        # pause experiment here if requested
        if thisExp.status == PAUSED:
            pauseExperiment(
                thisExp=thisExp, 
                win=win, 
                timers=[routineTimer], 
                playbackComponents=[]
            )
            # skip the frame we paused on
            continue
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            end.forceEnded = routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in end.components:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "end" ---
    for thisComponent in end.components:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # store stop times for end
    end.tStop = globalClock.getTime(format='float')
    end.tStopRefresh = tThisFlipGlobal
    thisExp.addData('end.stopped', end.tStop)
    # check responses
    if key_resp_3.keys in ['', [], None]:  # No response was made
        key_resp_3.keys = None
    thisExp.addData('key_resp_3.keys',key_resp_3.keys)
    if key_resp_3.keys != None:  # we had a response
        thisExp.addData('key_resp_3.rt', key_resp_3.rt)
        thisExp.addData('key_resp_3.duration', key_resp_3.duration)
    thisExp.nextEntry()
    # the Routine "end" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # mark experiment as finished
    endExperiment(thisExp, win=win)


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


def endExperiment(thisExp, win=None):
    """
    End this experiment, performing final shut down operations.
    
    This function does NOT close the window or end the Python process - use `quit` for this.
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    win : psychopy.visual.Window
        Window for this experiment.
    """
    if win is not None:
        # remove autodraw from all current components
        win.clearAutoDraw()
        # Flip one final time so any remaining win.callOnFlip() 
        # and win.timeOnFlip() tasks get executed
        win.flip()
    # return console logger level to WARNING
    logging.console.setLevel(logging.WARNING)
    # mark experiment handler as finished
    thisExp.status = FINISHED
    logging.flush()


def quit(thisExp, win=None, thisSession=None):
    """
    Fully quit, closing the window and ending the Python process.
    
    Parameters
    ==========
    win : psychopy.visual.Window
        Window to close.
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
    setupDevices(expInfo=expInfo, thisExp=thisExp, win=win)
    run(
        expInfo=expInfo, 
        thisExp=thisExp, 
        win=win,
        globalClock='float'
    )
    saveData(thisExp=thisExp)
    quit(thisExp=thisExp, win=win)
