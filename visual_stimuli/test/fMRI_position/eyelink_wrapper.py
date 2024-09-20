import os
import pylink


class EyelinkWrapper:
    def __init__(self, screen_size, output_file_name, color_bits=24, offline=False):
        self.screen_size = screen_size
        self.output_file_name = output_file_name
        self.color_bits = color_bits
        self.offline = offline
        self.eyelink = self._init_eyelink(screen_size, color_bits, offline)

    @staticmethod
    def _init_eyelink(screen_size, color_bits, offline):
        eyelink = pylink.EyeLink()
        if offline:
            eyelink.setOfflineMode()
        eyelink.sendCommand("screen_pixel_coords = 0 0 %d %d" % screen_size)
        eyelink.sendMessage("DISPLAY_COORDS  0 0 %d %d" % screen_size)
        eyelink.sendCommand("select_parser_configuration 0")
        eyelink.sendCommand("scene_camera_gazemap = NO")
        eyelink.sendCommand("pupil_size_diameter = %s" % ("YES"))
        pylink.openGraphics(screen_size, color_bits)
        pylink.setCalibrationColors((255, 255, 255), (0, 0, 0))
        pylink.setTargetSize(int(screen_size[0] / 70), int(screen_size[1] / 300))
        pylink.setCalibrationSounds("", "", "")
        pylink.setDriftCorrectSounds("", "off", "off")
        eyelink.doTrackerSetup(width=screen_size[0], height=screen_size[1])
        pylink.closeGraphics()
        return eyelink

    def start_recording(self):
        self.eyelink.openDataFile(self.output_file_name)
        self.eyelink.startRecording(1, 1, 1, 1)
        self.eyelink.sendMessage("Start")

    def finish_recording(self, output_dir):
        output_file_path = os.path.realpath(
            os.path.join(output_dir, self.output_file_name)
        )
        self.eyelink.stopRecording()
        self.eyelink.closeDataFile()
        os.makedirs(output_dir, exist_ok=True)
        self.eyelink.receiveDataFile(self.output_file_name, output_file_path)
        self.eyelink.close()
