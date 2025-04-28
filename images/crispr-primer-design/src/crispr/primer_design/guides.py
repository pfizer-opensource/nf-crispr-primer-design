"""Guide processing for primer design."""
import logging
import re
import sys

from .exceptions import PrimerDesignGuideError
from .subcommands import guide_bed_cmd, guide_flank_cmd
from .utils import reverse_complement


_logger = logging.getLogger(__name__)

def prepare_guides(input, reference_files, config):
    """Prepare guide input for primer design."""
    guides = process_guide_input(input,
                                 config['input_format'])
    _logger.info('Finding guide sequences in reference genome...')
    get_guide_coordinates(guides, reference_files, config)
    _logger.info('Calculating guide cut sites and genomic sequences...')
    for id, guide in guides.items():
        annotate_guide_site(guide, id)
    _logger.info('Extracting guide flanking genomic sequences...')
    get_guide_flanking_seq(guides, reference_files, config)
    return guides

def process_guide_input(input, input_format):
    """Process input data into a dict of guides."""
    input_rows = iter(input)
    header = next(input_rows)
    guides = {}
    for line in input_rows:
        fields = line.rstrip('\r\n').split('\t')
        try:
            id = re.sub(r'\s+', '', fields[input_format['id_col'] - 1])
            seq = re.sub(r'\s+', '', fields[input_format['seq_col'] - 1])
        except IndexError:
            _logger.warning('Input row could not be parsed: "%s". Skipping.', id)
            continue
        if not seq:
            _logger.warning('No target sequence found for guide "%s". Skipping.', id)
            continue
        guides[id] = {
            'seq': seq,
            'do_align': True,
        }
        try:
            guide_chr = fields[input_format['chr_col'] - 1].strip()
            guide_pos = fields[input_format['pos_col'] - 1].strip()
            guide_strand = fields[input_format['strand_col'] - 1].strip()
        except (KeyError, TypeError):
            pass
        else:
            if (guide_chr and guide_pos and guide_strand and
                  re.match(r'^[-+1]+$', guide_strand)):
                guides[id]['chr'] = guide_chr
                guides[id]['strand'] = '+' if guide_strand in ('+', '1') else '-'
                guide_pos_parts = guide_pos.split('-', 1)
                if len(guide_pos_parts) == 2:
                    guides[id]['start'] = int(guide_pos_parts[0])
                    guides[id]['end'] = int(guide_pos_parts[1])
                elif guide_strand in ('+', '1'):
                    guides[id]['start'] = int(guide_pos_parts[0]) - len(guide_seq) - 2
                    guides[id]['end'] = int(guide_pos_parts[0]) + 2
                else:
                    guides[id]['start'] = int(guide_pos_parts[0]) - 4
                    guides[id]['end'] = int(guide_pos_parts[0]) + len(guide_seq) - 4
                guides[id]['do_align'] = False
    _logger.debug('Processed guides: "%s"', guides)
    return guides

def get_guide_coordinates(guides, reference_files, config):
    """Get the genomic coordinates of each guide."""
    guide_bed = guide_bed_cmd(reference_files, config)
    guide_input = (f'>{id}\n{guide["seq"]}\n'
                   for id, guide in guides.items()
                   if guide['do_align'])
    for bed in guide_bed.communicate(''.join(guide_input)):
        fields = bed.rstrip('\r\n').split('\t')
        guides[fields[3]].update({
            'chr': fields[0],
            'start': int(fields[1]),
            'end': int(fields[2]),
            'strand': fields[5],
            'do_align': False,
        })
    if guide_bed.returncode != 0:
        raise PrimerDesignGuideError('Guide alignment failed.')

def annotate_guide_site(guide, id=None):
    """Determine the cut site and genomic sequence of a guide."""
    try:
        if guide['strand'] == '+':
            guide['cut_site'] = guide['end'] - 3
            guide['genomic_seq'] = guide['seq'].upper()
        else:
            guide['cut_site'] = guide['start'] + 3
            guide['genomic_seq'] = reverse_complement(guide['seq'].upper())
    except KeyError:
        _logger.warning(
            'Could not determine cut site or genomic sequence for guide: %s',
            id
        )

def get_guide_flanking_seq(guides, reference_files, config):
    """Get the flanking genomic sequence around each guide."""
    if config['parameters']['flanking'] < 0:
        # Flanking basepairs must not be negative
        config['parameters']['flanking'] = 0
    guide_flank = guide_flank_cmd(reference_files, config)
    bed_input = (f'{guide["chr"]}\t{guide["start"]}\t{guide["end"]}\t{id}\n'
                   for id, guide in guides.items()
                   if guide.get('chr'))
    # Store flanking sequence (first left flank, then right flank for each guide)
    for flanking in guide_flank.communicate(''.join(bed_input)):
        fields = flanking.rstrip('\r\n').split('\t')
        id = fields[0]
        if not guides[id].get('L_flank_seq'):
            guides[id].update({
                'L_flank_seq': fields[1].upper(),
                'flank_start': guides[id]['start'] - len(fields[1])
            })
        else:
            guides[id]['R_flank_seq'] = fields[1].upper()
    if guide_flank.returncode != 0:
        raise PrimerDesignGuideError('Extracting flanking sequence failed')
