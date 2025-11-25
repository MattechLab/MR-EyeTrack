#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Custom version: lastrun structure with composite dot stimulus and movement variables.
"""

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
import numpy as np
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle, choice as randchoice
import os
import sys
import math
import socket
import psychopy.iohub as io
from psychopy.hardware import keyboard

# device manager and paths
from psychopy import hardware as _hardware
deviceManager = _hardware.DeviceManager()
_thisDir = os.path.dirname(os.path.abspath(__file__))

psychopyVersion = '2024.2.1'
expName = 'fixed_dot-16_grid_T1w_custom'
expInfo = {
    'participant': f"{randint(0, 999):03.0f}",
    'session': '001',
    'date|hid': data.getDateStr(),
    'expName|hid': expName,
    'psychopyVersion|hid': psychopyVersion,
}

PILOTING = core.setPilotModeFromArgs()
_fullScr = True
_winSize = [800, 600]
if PILOTING and prefs.piloting['forceWindowed']:
    _fullScr = False
    _winSize = prefs.piloting['forcedWindowSize']


def showExpInfoDlg(expInfo):
    dlg = gui.DlgFromDict(dictionary=expInfo, sortKeys=False, title=expName, alwaysOnTop=True)
    if dlg.OK == False:
        core.quit()
    return expInfo


def setupData(expInfo, dataDir=None):
    for key, val in expInfo.copy().items():
        newKey, _ = data.utils.parsePipeSyntax(key)
        expInfo[newKey] = expInfo.pop(key)
    if dataDir is None:
        dataDir = _thisDir
    filename = u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])
    if os.path.isabs(filename):
        dataDir = os.path.commonprefix([dataDir, filename])
        filename = os.path.relpath(filename, dataDir)
    thisExp = data.ExperimentHandler(
        name=expName, version='',
        extraInfo=expInfo, runtimeInfo=None,
        originPath=__file__,
        savePickle=True, saveWideText=True,
        dataFileName=dataDir + os.sep + filename, sortColumns='time'
    )
    return thisExp


def setupLogging(filename):
    if PILOTING:
        logging.console.setLevel(prefs.piloting['pilotConsoleLoggingLevel'])
    else:
        logging.console.setLevel('warning')
    logFile = logging.LogFile(filename+'.log')
    if PILOTING:
        logFile.setLevel(prefs.piloting['pilotLoggingLevel'])
    else:
        logFile.setLevel(logging.getLevel('exp'))
    return logFile


def setupWindow(expInfo=None, win=None):
    if win is None:
        win = visual.Window(
            size=_winSize, fullscr=_fullScr, screen=0,
            winType='pyglet', allowStencil=False,
            monitor='testMonitor', color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb',
            backgroundImage='', backgroundFit='none',
            blendMode='avg', useFBO=True,
            units='norm', checkTiming=False
        )
    else:
        win.color = [-1.0, -1.0, -1.0]
        win.units = 'norm'
    if expInfo is not None:
        if win._monitorFrameRate is None:
            win._monitorFrameRate = win.getActualFrameRate()
        expInfo['frameRate'] = win._monitorFrameRate
    win.mouseVisible = False
    return win


def setupDevices(expInfo, thisExp, win):
    ioConfig = {}
    ioConfig['eyetracker.eyelink.EyeTracker'] = {
        'name': 'tracker',
        'model_name': 'EYELINK 1000 DESKTOP',
        'simulation_mode': True,
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
    ioConfig['Keyboard'] = dict(use_keymap='psychopy')
    ioConfig['Experiment'] = dict(filename=thisExp.dataFileName)
    ioServer = io.launchHubServer(window=win, **ioConfig)
    deviceManager.ioServer = ioServer
    deviceManager.devices['eyetracker'] = ioServer.getDevice('tracker')
    if deviceManager.getDevice('defaultKeyboard') is None:
        deviceManager.addDevice(deviceClass='keyboard', deviceName='defaultKeyboard', backend='iohub')
    # also create key_resp devices used in lastrun
    deviceManager.addDevice(deviceClass='keyboard', deviceName='key_resp')
    deviceManager.addDevice(deviceClass='keyboard', deviceName='key_resp_3')
    return True


def run(expInfo, thisExp, win, globalClock=None, thisSession=None):
    thisExp.status = STARTED
    exec = environmenttools.setExecEnvironment(globals())
    ioServer = deviceManager.ioServer
    defaultKeyboard = deviceManager.getDevice('defaultKeyboard')
    if defaultKeyboard is None:
        deviceManager.addDevice(deviceClass='keyboard', deviceName='defaultKeyboard', backend='ioHub')
    eyetracker = deviceManager.getDevice('eyetracker')
    os.chdir(_thisDir)
    filename = thisExp.dataFileName
    frameTolerance = 0.001
    endExpNow = False
    if 'frameRate' in expInfo and expInfo['frameRate'] is not None:
        frameDur = 1.0 / round(expInfo['frameRate'])
    else:
        frameDur = 1.0 / 60.0

    # --- Initialize components for Routine "trail" ---
    waiting_trigger = visual.TextStim(win=win, name='waiting_trigger', text="The program is ready for the scanner trigger. Press 's' to proceed manually.", pos=(0, -0.4), height=0.12, color='white')
    key_resp = keyboard.Keyboard(deviceName='key_resp')
    fix_desc = visual.TextStim(win=win, name='fix_desc', text='In this task you will see a dot moving randomly within different positions on the screen. You have to follow the dot :)', pos=(0, 0.25), height=0.12, color='white')

    # --- Initialize components for Routine "start_ET" ---
    etRecord = hardware.eyetracker.EyetrackerControl(tracker=eyetracker, actionType='Start Only')

    def send_message(message, addr="localhost", port=2023):
        client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        client_socket.connect((addr, port))
        client_socket.sendall(message)
        client_socket.close()

    # Geometry for composite dot (copied from original)
    outer_radius = 0.1 / 2
    inner_radius = 0.02 / 2
    cross_length = outer_radius * 2.1
    cross_thickness = inner_radius * 2.1

    aspect = float(win.size[0]) / float(win.size[1]) if hasattr(win, 'size') else 1.0
    x_scale = 1.0 / aspect

    def _circle_vertices(radius, n=64, x_scale=1.0):
        verts = []
        for i in range(n):
            theta = (2.0 * math.pi * i) / n
            x = math.cos(theta) * radius * x_scale
            y = math.sin(theta) * radius
            verts.append((x, y))
        return verts

    outer_verts = _circle_vertices(outer_radius, n=64, x_scale=x_scale)
    inner_verts = _circle_vertices(inner_radius, n=64, x_scale=x_scale)

    # Composite dot: a list of components (circles + cross rects)
    dot = []
    dot.extend([
        visual.ShapeStim(win=win, units='norm', vertices=outer_verts, fillColor=(-1, -1, -1), lineColor=(-1, -1, -1), interpolate=True),
        visual.Rect(win=win, units='norm', width=cross_length * x_scale, height=cross_thickness, fillColor=(0, 0, 0), lineColor=(0, 0, 0)),
        visual.Rect(win=win, units='norm', width=cross_thickness * x_scale, height=cross_length, fillColor=(0, 0, 0), lineColor=(0, 0, 0)),
        visual.ShapeStim(win=win, units='norm', vertices=inner_verts, fillColor=(-1, -1, -1), lineColor=(-1, -1, -1), interpolate=True),
    ])

    # Begin Experiment variables for movement
    grid_size = 4 # 4x4 grid
    dot_size = 0.05 # Size of the dot
    t_dot = 5*6.2/8.01 # Seconds of showing the dot per position. TR is decreased from 8.1 to 6.2ms

    # Get the screen dimensions
    # In norm units, screen goes from -1 to +1 vertically, and aspect-ratio-scaled horizontally.
    # So on an 800×600 display, full x-range is from -800/600 = -1.333 to +1.333
    # We normalize the pixel offsets to norm coordinates.
    screen_width, screen_height = win.size
    x_offset_norm = 1.33*2/3 # of 1.33
    y_offset_norm = 1/2 # of 1.00

    positions = {
        (0, y_offset_norm): "up",
        (0, -y_offset_norm): "down",
        (-x_offset_norm, 0): "left",
        (x_offset_norm, 0): "right"
    }

    ioServer.getDevice('tracker').sendMessage("ET: Start experiment 'dots'")

    # --- Initialize components for Routine "end" ---
    text = visual.TextStim(win=win, name='text', text="End press 't'", pos=(0, 0), height=0.12, color='white')
    ET_stop = hardware.eyetracker.EyetrackerControl(tracker=eyetracker, actionType='Stop Only')
    key_resp_3 = keyboard.Keyboard(deviceName='key_resp_3')

    # clocks
    if globalClock is None:
        globalClock = core.Clock()
    if isinstance(globalClock, str):
        if globalClock == 'float':
            globalClock = core.Clock(format='float')
        else:
            globalClock = core.Clock(format=globalClock)
    if ioServer is not None:
        ioServer.syncClock(globalClock)
    logging.setDefaultClock(globalClock)
    routineTimer = core.Clock()
    win.flip()
    expInfo['expStart'] = data.getDateStr(format='%Y-%m-%d %Hh%M.%S.%f %z', fractionalSecondDigits=6)

    # calibration (keep behavior same as lastrun)
    calibration_2Target = visual.TargetStim(win, name='calibration_2Target', radius=0.15, fillColor=[0.5,0.5,0.5], borderColor=[0.5,0.5,0.5], lineWidth=2.0, innerRadius=0.07, innerFillColor=[0.5,0.5,0.5])
    calibration_2 = hardware.eyetracker.EyetrackerCalibration(win, eyetracker, calibration_2Target, progressMode='time', targetDur=1.5, expandScale=1.5, targetLayout='FIVE_POINTS', randomisePos=True, textColor='white', movementAnimation=True, targetDelay=1.0)
    calibration_2.run()
    defaultKeyboard.clearEvents()
    routineTimer.reset()

    # --- trail routine (keep simplified) ---
    continueRoutine = True
    key_resp.keys = []
    key_resp.rt = []
    _key_resp_allKeys = []
    trailComponents = [waiting_trigger, key_resp, fix_desc]
    for thisComponent in trailComponents:
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
        if waiting_trigger.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            waiting_trigger.frameNStart = frameN
            waiting_trigger.tStart = t
            waiting_trigger.tStartRefresh = tThisFlipGlobal
            win.timeOnFlip(waiting_trigger, 'tStartRefresh')
            thisExp.timestampOnFlip(win, 'waiting_trigger.started')
            waiting_trigger.status = STARTED
            waiting_trigger.setAutoDraw(True)
        if key_resp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            key_resp.frameNStart = frameN
            key_resp.tStart = t
            key_resp.tStartRefresh = tThisFlipGlobal
            win.timeOnFlip(key_resp, 'tStartRefresh')
            thisExp.timestampOnFlip(win, 'key_resp.started')
            key_resp.status = STARTED
            win.callOnFlip(key_resp.clock.reset)
            win.callOnFlip(key_resp.clearEvents, eventType='keyboard')
        if key_resp.status == STARTED:
            theseKeys = key_resp.getKeys(keyList=['s'], ignoreKeys=["escape"], waitRelease=False)
            if len(theseKeys):
                key_resp.keys = theseKeys[-1].name
                key_resp.rt = theseKeys[-1].rt
                continueRoutine = False
        if fix_desc.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            fix_desc.frameNStart = frameN
            fix_desc.tStart = t
            fix_desc.tStartRefresh = tThisFlipGlobal
            win.timeOnFlip(fix_desc, 'tStartRefresh')
            thisExp.timestampOnFlip(win, 'fix_desc.started')
            fix_desc.status = STARTED
            fix_desc.setAutoDraw(True)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            return
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
    for thisComponent in trailComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    if key_resp.keys in ['', [], None]:
        key_resp.keys = None
    thisExp.addData('key_resp.keys',key_resp.keys)
    if key_resp.keys != None:
        thisExp.addData('key_resp.rt', key_resp.rt)
    thisExp.nextEntry()
    routineTimer.reset()

    # --- start_ET routine: preserve etRecord.start() ---
    continueRoutine = True
    start_ETComponents = [etRecord]
    for thisComponent in start_ETComponents:
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
        if etRecord.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            etRecord.frameNStart = frameN
            etRecord.tStart = t
            etRecord.tStartRefresh = tThisFlipGlobal
            win.timeOnFlip(etRecord, 'tStartRefresh')
            thisExp.timestampOnFlip(win, 'etRecord.started')
            etRecord.status = STARTED
            # Start eye-tracker recording
            eyetracker.sendMessage("ET: Start recording")
            etRecord.start()
        if etRecord.status == STARTED:
            etRecord.tStop = t
            etRecord.tStopRefresh = tThisFlipGlobal
            etRecord.frameNStop = frameN
            etRecord.status = FINISHED
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            return
        if not continueRoutine:
            routineForceEnded = True
            break
        continueRoutine = False
        for thisComponent in start_ETComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break
        if continueRoutine:
            win.flip()
    for thisComponent in start_ETComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    thisExp.nextEntry()
    routineTimer.reset()

    # --- Quick example of one "dots" routine using our composite dot list ---
    # We'll run a single iteration here as an example: the original has loops that reuse this pattern.
    shuffled_positions = list(positions.items())
    shuffle(shuffled_positions)
    current_position_index = 0
    total_positions = len(shuffled_positions)
    position, direction = shuffled_positions[current_position_index]
    for comp in dot:
        comp.setPos(position)
    time_of_last_change = routineTimer.getTime()
    continueRoutine = True
    dotComponents = [*dot]
    for thisComponent in dotComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    while continueRoutine:
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1
        # start all dot components together if not started
        if dotComponents and getattr(dotComponents[0], 'status', NOT_STARTED) == NOT_STARTED:
            for comp in dotComponents:
                comp.setAutoDraw(True)
        # move positions every t_dot seconds
        if t - time_of_last_change >= t_dot:
            current_position_index += 1
            if current_position_index < total_positions:
                position, direction = shuffled_positions[current_position_index]
                for comp in dot:
                    comp.setPos(position)
                ioServer.getDevice('tracker').sendMessage(f"ET: dot moved {direction}!")
                time_of_last_change = t
            else:
                continueRoutine = False
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            return
        if continueRoutine:
            win.flip()
    for comp in dot:
        comp.setAutoDraw(False)
    routineTimer.reset()

    # --- Ending sequence similar to lastrun ---
    ioServer.getDevice('tracker').sendMessage("ET: eye-tracker stopped")
    text.setAutoDraw(True)
    win.flip()
    # stop eyetracker if needed
    ET_stop.status = FINISHED
    eyetracker.sendMessage("ET: Stop eye tracker")
    thisExp.saveAsWideText(filename + '.csv', delim='auto')
    thisExp.saveAsPickle(filename)
    if eyetracker:
        eyetracker.setConnectionState(False)
        eyetracker.setRecordingState(False)
    thisExp.abort()
    win.close()
    core.quit()


if __name__ == '__main__':
    expInfo = showExpInfoDlg(expInfo=expInfo)
    thisExp = setupData(expInfo=expInfo)
    logFile = setupLogging(filename=thisExp.dataFileName)
    win = setupWindow(expInfo=expInfo)
    setupDevices(expInfo=expInfo, thisExp=thisExp, win=win)
    run(expInfo=expInfo, thisExp=thisExp, win=win, globalClock='float')
    # final save/quit handled inside run
