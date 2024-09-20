# Readme

## Experiment

The experiment is implemented using the builder of Psychopy-2023.1.2
The experiment starts with the S key (it is fundamental to synchronize the stimuli with the scans).

The audio stimulation consists of short sounds of objects of everyday life (e.g., a vacuum cleaner's sound, a child's crying, a guitar song, a bell ringing, etc.) or animal noises. With these sounds, 7 audio tracks (using Audicity software) were created sequentially and repeated during the experiment.
The experimental protocol consists of 15 seconds of audio stimulation using the above tracks (corresponding to the on period) followed by 25 seconds of pauses in which we provided a pink noise (corresponding to the off period). An on/off period represents a trial. In total, we provided 49 trials.

For the entire experiment period, the participant takes their eye open on a fixation point (a red cross) at the center of a gray environment.

## SOPs

1. Assets: contain the material needed for the stimulation. Importantly: when you downloaded the stim_audio folder, you need to change the new path in the experimental semantic_audio_fMRI psychopy envirorment.
2. Data_collection: add 5 files

    - participant_preparation
    - pre_session
    - preliminary
    - scanning
    - tear_down

3. Subject_Management: add 2 files

    - Recruitment
    - Scheduling

Aggiungere come files:

1. changes
2. data_storage
3. index
4. preprocessing
5. release
