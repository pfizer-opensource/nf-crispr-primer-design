"""Primer design and processing."""
import logging

import primer3


_logger = logging.getLogger(__name__)

designer = f'Primer3 v{primer3.__version__}'

def design_primers(guides, params):
    """Design primers for each guide sequence."""
    _logger.info('Designing primers...')
    primers = {}
    
    size_range = [list(map(int, sr.split('-'))) for sr in params['size_range']]
    if params['num_pairs'] < 1:
        params['num_pairs'] = 1
        _logger.warning('Number of primer pairs must be at least 1. '
                        'Designing 1 primer pair per target site.')
    primer3.setP3Globals({
        'PRIMER_OPT_SIZE': params['primer_opt'],
        'PRIMER_MIN_SIZE': params['primer_min'],
        'PRIMER_MAX_SIZE': params['primer_max'],
        'PRIMER_PRODUCT_OPT_SIZE': params['size_opt'],
        'PRIMER_PRODUCT_SIZE_RANGE': size_range,
        'PRIMER_NUM_RETURN': params['num_pairs'],
        'PRIMER_MAX_POLY_X': params['max_poly_x'],
    })
    
    for id, guide in guides.items():
        primers[id] = {'guide': guide}
        if not (guide.get('L_flank_seq') and guide.get('R_flank_seq')):
            primers[id]['pairs'] = []
            continue
        # Sequence to search for primers
        seq = ''.join((guide['L_flank_seq'],
                       guide['genomic_seq'],
                       guide['R_flank_seq']))
        # Target region is centered on the guide sequence. The window cannot
        # be a negative number.
        target_window = params['target'] if params['target'] > 0 else 0
        target_flank = int((target_window - len(guide['seq'])) / 2)
        target  = (len(guide['L_flank_seq']) - target_flank, target_window)
        primer_info = primer3.designPrimers({
            'SEQUENCE_ID': id,
            'SEQUENCE_TEMPLATE': seq,
            'SEQUENCE_TARGET': target,
        })
        primer_pairs = process_primers(primer_info, params['num_pairs'])
        calc_primer_coordinates(primer_pairs, guide['flank_start'], guide['cut_site'])
        primers[id]['pairs'] = primer_pairs
    flag_primer_duplicates(primers)
    return primers

def process_primers(primers, max_num_pairs):
    """Process primer3 output into a list of primer pairs with needed info."""
    pairs = []
    for i in range(0, max_num_pairs):
        try:
            pairs.append({
                'pair_number': i + 1,
                'total_pairs': primers['PRIMER_PAIR_NUM_RETURNED'],
                'F_primer_seq': primers[f'PRIMER_LEFT_{i}_SEQUENCE'],
                'F_primer_offset': primers[f'PRIMER_LEFT_{i}'][0],
                'F_primer_length': primers[f'PRIMER_LEFT_{i}'][1],
                'F_primer_tm': primers[f'PRIMER_LEFT_{i}_TM'],
                'F_primer_gc': primers[f'PRIMER_LEFT_{i}_GC_PERCENT'],
                'R_primer_seq': primers[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                'R_primer_offset': primers[f'PRIMER_RIGHT_{i}'][0],
                'R_primer_length': primers[f'PRIMER_RIGHT_{i}'][1],
                'R_primer_tm': primers[f'PRIMER_RIGHT_{i}_TM'],
                'R_primer_gc': primers[f'PRIMER_RIGHT_{i}_GC_PERCENT'],
                'amplicon_size': primers[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                'penalty': primers[f'PRIMER_PAIR_{i}_PENALTY'],
            })
        except KeyError:
            pass
    return pairs

def calc_primer_coordinates(primer_pairs, genome_offset, cut_site):
    """Calculate primer genomic coordinates based on a genome offset value."""
    for pair in primer_pairs:
        # Start/end are genomic coordinates in BED format (ie, 0-indexed
        # and end does not include last base)
        pair['F_primer_start'] = genome_offset + pair['F_primer_offset']
        pair['F_primer_end'] = pair['F_primer_start'] + pair['F_primer_length']
        pair['R_primer_end'] = genome_offset + pair['R_primer_offset'] + 1
        pair['R_primer_start'] = pair['R_primer_end'] - pair['R_primer_length']
        # Calculate minimum distance from end of primers to guide cut site
        forward_distance = cut_site - pair['F_primer_start']
        reverse_distance = pair['R_primer_end'] - cut_site
        pair['dist_to_cut'] = min(forward_distance, reverse_distance)

def flag_primer_duplicates(primers):
    duplicates = {}
    for guide_primers in primers.values():
        for pair in guide_primers['pairs']:
            duplicates[pair['F_primer_seq']] = duplicates.get(pair['F_primer_seq']) is not None
            duplicates[pair['R_primer_seq']] = duplicates.get(pair['R_primer_seq']) is not None
    for guide_primers in primers.values():
        for pair in guide_primers['pairs']:
            pair['duplicates'] = '{}{}'.format(
                'F' if duplicates[pair['F_primer_seq']] else '',
                'R' if duplicates[pair['R_primer_seq']] else ''
            )

def make_primers_table(primers, max_num_pairs, input=None, add_header=True):
    """Return a generator that provides each row of the primer table as a
    list of fields."""
    if input is None:
        input_rows = None
        header_names = ['Guide ID', 'Guide seq']
    else:
        input_rows = iter(input)
        header_names = next(input_rows).rstrip('\r\n').split('\t')
    if max_num_pairs > 1:
        header_names.append('Primer pair')
    header_names.extend(('F primer seq', 'F primer Tm', 'F primer %GC',
                         'R primer seq', 'R primer Tm', 'R primer %GC',
                         'Amplicon size', 'Penalty score',
                         'Min dist to cut', 'Flag duplicates'))
    if add_header:
        yield header_names
    for id, guide_primers in primers.items():
        if input_rows is None:
            row_start = [id, guide_primers['guide']['seq']]
        else:
            row_start = next(input_rows).rstrip('\r\n').split('\t')
        if guide_primers['pairs']:
            for pair in guide_primers['pairs']:
                row = row_start[:]
                if max_num_pairs > 1:
                    row.append(pair['pair_number'])
                row.extend((pair['F_primer_seq'], f'{pair["F_primer_tm"]:.3f}',
                            f'{pair["F_primer_gc"]:.3f}', pair['R_primer_seq'],
                            f'{pair["R_primer_tm"]:.3f}', f'{pair["R_primer_gc"]:.3f}',
                            pair['amplicon_size'], f'{pair["penalty"]:.6f}',
                            pair['dist_to_cut'], pair['duplicates']))
                yield row
        else:
            row_start.extend([''] * (len(header_names) - len(row_start)))
            yield row_start

def get_amplicon_name(guide_id, primer_pair):
    """Return the amplicon name based on the guide ID and primer pair number.
    
    Amplicon name is the same as guide ID if only one pair was designed, or
    it will be "[guide ID]-pair[pair number]" if multiple pairs were designed.
    """
    if primer_pair['total_pairs'] > 1:
        return f'{guide_id}-pair{primer_pair["pair_number"]}'
    else:
        return guide_id

def make_primers_bed_table(primers, header=True):
    """Return a generator that provides each row of a BED-formatted file
    for each primer pair."""
    if header:
        # Header for UCSC genome browser
        yield ['track name="CRISPR primers" '
               'description="CRISPR guides and sequencing primers" '
               'visibility=pack useScore=1']
    # BED output for each guide and amplicon
    for id, guide_primers in primers.items():
        for pair in guide_primers['pairs']:
            # Calculate offsets for guide and R primer from start of amplicon
            guide_offset = guide_primers['guide']['start'] - pair['F_primer_start']
            R_primer_offset = pair['R_primer_start'] - pair['F_primer_start']
            amplicon_name = get_amplicon_name(id, pair)
            # First 3 fields are the genomic coordinates of entire amplicon
            row = [
                guide_primers['guide']['chr'],
                str(pair['F_primer_start']),
                str(pair['R_primer_end'])
            ]
            # Next field is the amplicon name
            row.append(amplicon_name)
            # Convert the penalty score (lower is better) into shading
            # from 0-1000 (1000 is best)
            row.append(str(round(1000 - round(pair['penalty'], 6) * 50, 6)))
            # The next field is for strand and is not used
            row.append('.')
            # Set thick start and end to the guide location, primers will appear thin
            row.extend((
                str(guide_primers['guide']['start']),
                str(guide_primers['guide']['end'])
            ))
            # The next field is for itemRGB and is not used
            row.append('0')
            # Display the two primers and the guide as 3 separate blocks
            row.append('3')  # blockCount
            row.append(','.join(map(str, (pair['F_primer_length'],
                                          len(guide_primers['guide']['seq']),
                                          pair['R_primer_length']))))  # blockSizes
            row.append(f'0,{guide_offset},{R_primer_offset}')  # blockStarts
            yield row
