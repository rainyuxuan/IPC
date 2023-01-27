import os
import re
import glob

REGEX_TO_CLEANUP = [
    '^.*\.obj$',
    '^info\d+.txt$',
    '^status\d+$',
]

DIR_TO_EXCLUDE = []

def filepath(abs_path: str=''):
    """Get the path to a file relative to the current working directory.

    :param abs_path: path relative to the root directory.
    :return:
    """
    return '../output/' + abs_path


if __name__ == '__main__':
    # dir_list = glob.glob(filepath(''))
    dir_list = [f for f in os.listdir(filepath()) if f != '__images']

    for directory in dir_list:
        dir_files = os.listdir(filepath(directory))
        for file in dir_files:
            for regex in REGEX_TO_CLEANUP:
                if re.match(regex, file):
                    path = filepath(f'{directory}/{file}')
                    os.remove(path)
