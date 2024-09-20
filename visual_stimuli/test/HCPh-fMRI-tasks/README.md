# Functional MRI tasks of the Human Connectome PHantom (HCPh) study

This repository contains three functional MRI tasks implemented as *Psychopy3* experiments:

* Breath-holding task (BHT): [`task-bht_bold.psyexp`](https://github.com/TheAxonLab/HCPh-fMRI-tasks/blob/master/task-bht_bold.psyexp)
* Quality-control task (QCT): [`task-qct_bold.psyexp`](https://github.com/TheAxonLab/HCPh-fMRI-tasks/blob/master/task-qct_bold.psyexp)
* Resting-state fMRI (rest): [`task-rest_bold.psyexp`](https://github.com/TheAxonLab/HCPh-fMRI-tasks/blob/master/task-rest_bold.psyexp)

## Quality control task
The task contains a visual trial block, an eye-movement trial block, and a finger-tapping trial block.
The task was adapted from «Version A» of the task proposed by Harvey et al.:

     Harvey J, Demetriou L, McGonigle J, Wall MB. 2018. A short, robust brain activation
     control task optimised for pharmacological fMRI studies. PeerJ 6:e5540
     doi:https://doi.org/10.7717/peerj.5540

The original implementation of their task is found at [mattwall1103/fMRI-Control-Task](https://github.com/mattwall1103/fMRI-Control-Task).

The most prominent changes to the original task are:

- Elimination of the auditive block (since we do not plan the parcipant to wear headphones)
- Substitution of the motor task (button pressing) with a simple finger-tapping paradigm
- Enabling eye-tracking with our SR instruments device.

### Detailed description of QCT

The QCT consists of four paradigms whose presentation order and realization are randomly selected.
The background for all tasks is black (HEX #000000).
The units reported here correspond to *Psychopy*'s "normalized units".

1. Visual gratting pattern: Visual trials consist of a centrally-presented sine-wave grating subtending approximately 10° of visual angle and with a spatial frequency of 1.2 cycles/degree.
    The grating drifts laterally at a rate of 6 cycles per second, and the direction of the drift is reversed every 0.5s (Harvey et al. 2018).
    A small green, circular fixation point is also displayed in the center of the screen.
    The total duration of the stimulus is 3s.

2. Fingertapping: Motor trials correspond to a simple finger-tapping task of the left or right hand.
    The participant is instructed to tap their thumb on each of their other four fingers of the hand indicated on the screen (presenting the words 'LEFT' or 'RIGHT'), sequentially with all fingers and reversing the direction at the extremes (from pointer to pinkie).
    The words 'LEFT' or 'RIGHT', are presented in random order, in white color, and 0.5 units to the left from the center of the screen, respectively to the right, for 5s.

3. Gaze movement: Cognitive trials involve a series of fixation points moving across the screen.
    Each fixation point comprises a small green dot (HEX #00ff00) surrounded by a larger concentric gray circle (radius ratio is ~2.2).
    Participants are instructed to focus their gaze on the center of the fixation point and follow it with their eyes while avoiding head movements.
    The fixation point moves to six different locations corresponding to the compass directions North–East, East, South–East, North–West, West, and South–West.
    These points are mapped on a circle with a radius of approximately 8.75° of visual angle.
    Each location is maintained for 0.5 s, and all six are presented (in random order) in each three-second trial (Harvey et al. 2018).

4. Fixation: Blank trials consist of a fixation point presented at the center of the screen for 3s.
    The fixation point is built like in the gaze movement task.

# License

These tasks are released under the terms of the Apache 2.0, in order to abide by the [NiPreps licensing principles](https://www.nipreps.org/community/licensing/). See ``NOTICE`` file for further details.
