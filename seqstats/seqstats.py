"""
Given a reference sequence and a pairs of zero-based coordinates, 
compute and report various statistics about the sequences.
"""
import sys
from math import log 
from argparse import ArgumentParser
import logging

from Bio import SeqIO


DEFAULT_LOG_FILE = "stderr"
DEFAULT_CONTEXT_SIZE = 500
DEFAULT_WINDOW_SIZE = 50


UNNAMED = "Unnamed"

#################### Argument Parser, Logging and Main #######################

def parse_args():
    "Parse command line argumetns for the program"
    parser = ArgumentParser(description="Generate statistics for a sequence "
                                        "and its context in a reference "
                                        "sequence")
    parser.add_argument("--log", metavar="FILE",
                        type=str,
                        default=DEFAULT_LOG_FILE,
                        help="A file used for program logging. Default: %s"
                             %DEFAULT_LOG_FILE)
    parser.add_argument("--ref_file", metavar="FILE",
                        type=str,
                        help="A fasta file containing the reference sequence "
                             "Alternatively, a single reference sequence can "
                             "specified by the --ref_seq option")
    parser.add_argument("--ref_seq", metavar="SEQ",
                        type=str,
                        help="The reference sequence used. Alternatively, "
                             "A list of reference sequence can be specified "
                             "by passing a fasta file via the --ref_file "
                             "option.")
    parser.add_argument("--coords", metavar="COORD",
                        type=int,
                        nargs=2,
                        help="Zero based coordinate of the sequence to " 
                             "be analysed. Either this specify coords or "
                             "specify a file containing all the "
                             "coordinate of the sequence to be analysed "
                             "via the --bed option.")
    parser.add_argument("--bed", metavar="FILE",
                        type=str,
                        help="A bed file containing coordinates "
                             "of sequence to "
                             "be analysed. e.g. \n"
                             "seq\t 00500\t100600\nseq\t115500\t11560\n.\n.\n"
                             "Alternatively, a single pair of "
                             "coordinate can be specified via the --coords "
                             "option.")
    parser.add_argument("--context_size", metavar="N",
                        type=int,
                        default=DEFAULT_CONTEXT_SIZE,
                        help="The number of sequence bases before and after "
                             "the specified the sequence in the reference "
                             "sequence to be considered as 'sequence context' "
                             "in the analysis. \n"
                             "Default: %s"%DEFAULT_CONTEXT_SIZE)
    #!! "Whenever used" need to be defined afterwards.
    parser.add_argument("--window_size", metavar="N",
                        type=int,
                        default=DEFAULT_WINDOW_SIZE,
                        help="The size of the sliding window (whenever used) "
                             "in the analysis. \n"
                             "Default: %s"%DEFAULT_WINDOW_SIZE)

    return parser.parse_args()


def start_log(log, level=logging.DEBUG):
    """
    Initiate program logging. If no log file is specified then log output
    goes to stderr.
    """
    if log == "stderr":                                                
        log = sys.stderr                                                       
    logging.basicConfig(stream = log,                                          
                        level = level,                                         
                        filemode = 'w',                                        
                        format = "%(asctime)s %(message)s",                    
                        datefmt = "[%m/%d/%Y %H:%M:%S] ")                      
    logging.info("Bug Sensor Program started")                                 
    logging.info("Command line: {0}\n".format(' '.join(sys.argv)))

def main():
    """
    main function. Entry point for the program.
    """
    user_inputs = parse_args()
    start_log(user_inputs.log)

    # Obtain a list of coordinates from user_inputs
    coords_dict = {}
    if user_inputs.coords:
        start, end = user_inputs.coords
        name = UNNAMED 
        coords_dict[UNNAMED] = (start, end)
    if user_inputs.bed:
        with open(user_inputs.bed) as bed_file:
            reader = csv.reader(bed_file, delimiter='\t')
            for line in reader:
                chrom, start, end = line
                start = int(start)
                end = int(start)
                name = "%s_%i_%i"%(chrom, start, end)
                coord_dict[name] = (start, end)

    # obtain all reference sequences from user_inputs
    ref_dict = {}
    if user_inputs.ref_seq:
        ref_dict[UNNAMED] = user_inputs.ref_seq
    if user_inputs.ref_file:
        name_seq_dict = ref_seq_from_file(user_inputs.ref_file)
        if UNNAMED in name_seq_dict:
            logging.warning("WARNING: Reference sequence naming clash")
        ref_dict.update(name_seq_dict)

    if set(ref_dict.keys()) != set(coords_dict.keys()):
        logging.warning("WARNING: inconsistent naming accross "
                        "reference sequences and coordinate names.")

    for name in ref_dict:
        ref_seq = ref_dict[name]
        start, end = coords_dict[name]
        seq_stats = Seqstats(ref_seq, 
                             start, 
                             end, 
                             name,
                             user_inputs.context_size, 
                             user_inputs.window_size)
        seq_stats.sequence_entropy_stats()
        seq_stats.context_entropy_stats()
        seq_stats.print_stats()

    return 0 




########################## Helper functions #################################
    

def entropy(seq, case_sensitive=False):
    """
    Calculate the (information theoretic) entropy of the given string.
    """
    if not case_sensitive:
        seq = seq.upper()
    length = float(len(seq))
    probs = [seq.count(char) / length for char in set(seq)]
    return sum(p*log(p) for p in probs)

def ref_seq_from_file(ref_file):
    """
    Return a dictionary of name:reference_sequence by parsing the
    input fasta file "ref_file".
    """
    ref_parser = SeqIO.parse(ref_file, 'fasta')
    name_seq_dict = {}
    unnamed_counter = 0
    for read in ref_parser:
        name = read.name
        if name:
            name_seq_dict[name] = read.seq
        else:
            unnamed_counter += 1
            name_seq_dict["%s%i"%(UNNAMED, unnamed_counter)] = read.seq
            logging.warning("WARNING: Unnamed reference sequence.")
    return name_seq_dict

def window_generator(sequence, window_size):
    seq_len = len(sequence)
    if window_size > seq_len:
        logging.critical("CRITICAL: Sequence size smaller than window size")
    for i in range(seq_len - window_size):
        yield sequence[i : i + window_size]

def mean(num_list):
    """
    Return the average of a list of number. Return 0.0 if empty list.
    """
    if not num_list:
        return 0.0
    else:
        return sum(num_list) / float(len(num_list))


def get_context(ref_seq, start, end, context_size):
    """
    Given the reference sequence and the coordinates of the 
    embedded sequence of concern, obtain the sequence before 
    and after the sequence (i.e. the flanking "context").
    The size of the context depends on the specifed "context size".
    """
    if start <= context_size:
        logging.warning("WARNING: prior context will be shorter "
                        "than context_size specified")
    ref_len = len(ref_seq)
    if end + context_size > ref_len:
        logging.warning("WARNING: posterior context will be shorter than "
                        "context_size specified")
    prior_start = max(start - context_size, 0)
    prior_seq = ref_seq[prior_start : start]
    posterior_end = min(end + context_size, ref_len)
    posterior_seq = ref_seq[end: posterior_end]
    return (prior_seq, posterior_seq)

########################### Seqstats Object ##################################

class Seqstats(object):
    """
    #!! Documentation
    """
    def __init__(self, ref_seq, start, end, name, context_size, window_size):
        self.ref_seq = ref_seq
        self.coords = (start, end)
        self.name = name
        self.window_size = window_size
        self.context_size = context_size
        
        self.sequence = self.ref_seq[start : end]
        self.context = get_context(self.ref_seq, 
                                   self.coords[0],
                                   self.coords[1],
                                   self.context_size)
        self.stats_dict = {}


    def sequence_entropy_stats(self):
        stats_dict = {} 
        stats_dict["overall_entropy"] = entropy(self.sequence)
        windows = window_generator(self.sequence, self.window_size)
        window_entropies = [entropy(seq) for seq in windows]
        avg_entropy = mean(window_entropies)
        stats_dict["avg_sequence_windowed_entropy"] = avg_entropy
        self.stats_dict.update(stats_dict)

    def context_entropy_stats(self):
        stats_dict = {}
        prior, posterior = self.context
        stats_dict["overall_prior_context_entropy"] = entropy(prior)
        stats_dict["overall_posterior_context_entropy"] = entropy(posterior)
        prior_windows = window_generator(prior, self.window_size)
        avg_prior_entropy = mean([entropy(prior_seq) 
                                  for prior_seq in prior_windows])
        posterior_windows = window_generator(posterior, self.window_size)
        avg_posterior_entropy = mean([entropy(posterior_seq) 
                                      for posterior_seq in posterior_windows])
        stats_dict["avg_prior_context_windowed_entropy"] = avg_prior_entropy
        stats_dict["avg_posterior_context"
                  "_windowed_entropy"] = avg_posterior_entropy
        self.stats_dict.update(stats_dict)
    
    def print_stats(self):
        """
        Print a report for the statistics gathered in Seqstats.stats_dict
        """
        report_string = "Report for %s \n"%self.name
        for field_name, value in self.stats_dict.items():
            report_string += "%s\t= %s\n"%(field_name, str(value))
        print report_string


if __name__ == "__main__":
    main()
