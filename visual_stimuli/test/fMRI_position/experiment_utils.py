from copy import copy
import os
from psychopy import core, gui


def prompt_for_params(*fields, title=""):
    dlg = gui.Dlg(title=title)
    for field in fields:
        dlg.addField(*field[1:])
    dlg.show()
    if not dlg.OK:
        raise Exception("Experiment canceled")
    return dict(zip((f[0] for f in fields), dlg.data))


class ExperimentLogCSV:
    def __init__(self, name, params):
        self.name = name
        self.params = params
        self.output_path = os.path.realpath(
            os.path.join(
                os.path.dirname(__file__),
                "results",
                "%s_%s.csv" % (self.name, self.params["subject_id"]),
            )
        )
        self.timer = core.Clock()
        self.buffer = []

    def start(self):
        if os.path.isfile(self.output_path):
            raise Exception("Output file already exists at: " + self.output_path)
        if not os.path.exists(os.path.dirname(self.output_path)):
            os.makedirs(os.path.dirname(self.output_path))
        self.file = open(self.output_path, "w")
        self.write_file_header()
        self.timer.reset()

    def write_file_header(self):
        params = copy(self.params)
        self.buffer.append(("Experiment", self.name))
        self.buffer.append(("Subject Id", params.pop("subject_id")))
        for param in params.items():
            self.buffer.append(param)
        self.buffer.append((" ",))
        self.flush()

    def log(self, event, *details):
        self.buffer.append(
            (self.timer.getTime(), event, *(detail for detail in details))
        )

    def flush(self):
        for line in self.buffer:
            replace_args = ('"', r"\"")
            self.file.write(
                ",".join(f'"{str(cell).replace(*replace_args)}"' for cell in line)
                + "\r\n"
            )
        self.file.flush()
        self.buffer = []

    def finish(self):
        self.flush()
        self.file.close()
