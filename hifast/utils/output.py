

__all__ = ['CustomStdout', 'set_output']


import os
import sys
from .colors import Colors as c


class CustomStdout:
    def __init__(self, original_stdout, additional_info=''):
        self.original_stdout = original_stdout
        self.additional_info = additional_info
        self.additional_info = f'{c.bold}{self.additional_info}{c.end}'
        self.lineends = True

    def _write_with_additional_info(self, message):
        if message.strip() in ['\n', ''] or message.startswith('File exists') or message.startswith('Saved to'):
            self.original_stdout.write(message)
        else:
            self.original_stdout.write(f'{self.additional_info}{message}')

    def write(self, message):
        # tqdm start with '\r' but not with '\n'
        if message.startswith('\r'):
            message = message[1:]  # Properly handle messages that intend to overwrite the line
            self.original_stdout.write(f'\r{self.additional_info}{message}')
            # or not add additional info?
            self.lineends = True
        else:
            if not message.endswith('\n'):
                if self.lineends:
                    self._write_with_additional_info(message)
                else:
                    self.original_stdout.write(message)
                self.lineends = False
            else:
                if self.lineends:
                    self._write_with_additional_info(message)
                else:
                    self.original_stdout.write(message)
                self.lineends = True


    def flush(self):
        self.original_stdout.flush()

    def __getattr__(self, attr):
        """ Delegate access to other attributes and methods to the original stdout """
        return getattr(self.original_stdout, attr)

def set_output(additional_info):
    # if HIFAST_MAKE_OUTPUT_CLEAR is set in environment variable and 1
    if 'HIFAST_MAKE_OUTPUT_CLEAR' in os.environ and os.environ['HIFAST_MAKE_OUTPUT_CLEAR'] == '1':
        from tqdm import tqdm
        from functools import partialmethod
        # after add the additional to the progress bar in tqdm, the linewith will exceed the terminal width, so we need to set nclos to 0, a small number
        # or import shutil;shutil.get_terminal_size().columns - len(additional_info) - extra_to_make_sure_it_will_not_exceed_the_terminal_width
        tqdm.__init__ = partialmethod(tqdm.__init__, ncols=0, mininterval=5)

        sys.stdout = CustomStdout(sys.stdout, additional_info)
        sys.stderr = CustomStdout(sys.stderr, additional_info)
