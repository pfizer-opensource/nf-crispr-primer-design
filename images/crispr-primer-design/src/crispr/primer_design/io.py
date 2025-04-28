"""Input and output file handling."""
import logging
from pathlib import Path
import sys

from .exceptions import PrimerDesignIOError


_logger = logging.getLogger(__name__)


class InputReader:
    """An iterator for reading input data that stores non-file input."""
    
    def __init__(self, input_file, prompt=None):
        self.input_file = input_file
        self.prompt = prompt
        self._stored_input = None
    
    def __iter__(self):
        if self.input_file == '-':
            if self._stored_input:
                for line in self._stored_input:
                    yield line
            else:
                if self.prompt:
                    print(self.prompt, file=sys.stderr)
                self._stored_input = []
                for line in sys.stdin:
                    self._stored_input.append(line)
                    yield line
        else:
            try:
                with open(self.input_file, 'r') as input:
                    for line in input:
                        yield line
            except OSError as e:
                raise PrimerDesignIOError(
                    f'Could not open input file "{self.input_file}": {e}'
                ) from None


def prep_output_dir(output_file):
    output_dir = Path(output_file).parent
    try:
        output_dir.mkdir(exist_ok=True)
    except OSError as e:
        raise PrimerDesignIOError(
            f'Directory "{output_dir}" for output file is not accessible: {e}'
        ) from None
    return output_dir

def write_tsv_table(output_file, table):
    try:
        with open(output_file, 'w') as output:
            for row in table:
                print(*row, sep='\t', file=output)
    except OSError as e:
        raise PrimerDesignIOError(
            f'Could not open output file "{output_file}": {e}'
        ) from None
