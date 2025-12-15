

__all__ = ['ObsLogParser', 'get_sep_para_from_log', 'search_log_fpath']


import re
import io


class ObsLogParser:
    def __init__(self, fpath_or_content):
        """
        fpath_or_string: str
           If `fpath_or_string` ends with `.log` or `.txt`, it will be treated as a file path. Otherwise, the content will be treated as a log.
        """
        if fpath_or_content.endswith('.log') or fpath_or_content.endswith('.txt'):
            self.fpath = fpath_or_content
            # Read the contents of the file
            with open(self.fpath, 'r') as file:
                text = file.read()
        else:
            self.fpath = None
            text = fpath_or_content


        # Initialize the dictionary
        self.D = {}
        self.text = text

    def __call__(self):
        return self.parse_log()

    def parse_log(self):

        self.extract_type1_info()
        self.extract_table_info()
        self.extract_obs_time()
        self.extract_obs_mode()

        return self.D

    def extract_type1_info(self):
        text = self.text
        D = self.D
        # Extract information using regular expressions
        pattern = r'(\w+)\s*[:：=]\s*([\-\w\d.]+)\s*'
        matches = re.findall(pattern, text)

        # Assign the extracted information to the dictionary
        for match in matches:
            key = match[0]
            value = match[1]
            D[key] = value

    def extract_table_info(self):
        text = self.text
        D = self.D
        # Extract information using regular expressions
        table_pattern = r'后端\s+采样时间\s+波束\s+极化\s+谱线窄带中心频率\s*[\r\n]+([^\r\n]+)'
        table_match = re.search(table_pattern, text)

        if table_match:
            table_text = table_match.group(0)

            # Extract rows from the table and remove leading/trailing spaces
            rows = [row.strip() for row in table_text.split('\n')]

            # Extract header row and remove leading/trailing spaces
            header = rows[0].split()
            values = rows[1].split()
            for i in range(len(header)):
                D[header[i]] = values[i]

    def extract_obs_time(self):

        pattern = r"(\d{4}/\d{2}/\d{2} \d{2}:\d{2}:\d{2}) - (\d{4}/\d{2}/\d{2} \d{2}:\d{2}:\d{2})"

        matches = re.findall(pattern, self.text)
        try:
            self.D['obs_time_range'] = matches[0]
        except err:
            print(err)
            self.D['obs_time_range'] = None

    def extract_obs_mode(self):
        pattern = r"(\w+)\s+(模式观测)\s+(\w+)"
        matches = re.findall(pattern, self.text)
        try:
            self.D['obs_mode'] = matches[0]
        except:
            self.D['obs_mode'] = None


def get_sep_para_from_log(fpath_or_content):
    """
    fpath_or_string: str
        If `fpath_or_string` ends with `.log` or `.txt`, it will be treated as a file path. Otherwise, the content will be treated as a log.
    """
    obs_log = ObsLogParser(fpath_or_content)()
    # unit as 4ns
    t_smaple = float(obs_log['采样时间'][:-1]) * 251658240
    t_smaple = int(t_smaple)

    para = {}
    n_delay = int(obs_log['Delay']) / t_smaple
    n_on = int(obs_log['On']) / t_smaple
    n_off = int(obs_log['Off']) / t_smaple

    if not n_delay.is_integer():
        raise(ValueError('The delay time should be an integer multiple of the sample time.'))
    if not n_on.is_integer():
        raise(ValueError('The Noise-On time should be an integer multiple of the sample time.'))
    if not n_off.is_integer():
        raise(ValueError('The Noise-Off should be an integer multiple of the sample time.'))

    para['n_delay'] = int(n_delay)
    para['n_on'] = int(n_on)
    para['n_off'] = int(n_off)
    para['noise_mode'] = obs_log['噪声强度']

    command = ''
    for key in para.keys():
        command += f' --{key} {para[key]}'

    return para, command


import os
from glob import glob
from .io import get_project

def search_log_fpath(fpath_fits, log_dir):
    project = get_project(fpath_fits)
    project_b = project.rsplit('_',1)[0]

    date = os.path.basename(os.path.dirname(fpath_fits))
    date = f"{date[:4]}-{date[4:6]}-{date[6:8]}"

    pattern = f"{log_dir}/{project_b}..{date}..*.log.txt"
    logfiles = glob(pattern)
    if len(logfiles) == 0:
        return None
    elif len(logfiles) == 1:
        return logfiles[0]
    else:
        return None
