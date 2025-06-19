import logging
from Bio.Blast import NCBIXML # For parsing

logger = logging.getLogger(__name__) # cds_processor.rbh_handler

class RBHHit:
    """Stores information about a single BLAST hit relevant for RBH."""
    def __init__(self, query_id, subject_id, pident, length, evalue, bitscore,
                 query_len=0, subject_len=0, query_start=0, query_end=0,
                 subject_start=0, subject_end=0):
        self.query_id = str(query_id)
        self.subject_id = str(subject_id)
        self.pident = float(pident)
        self.length = int(length)
        self.evalue = float(evalue)
        self.bitscore = float(bitscore)
        self.query_len = int(query_len) # Original length of the query sequence
        self.subject_len = int(subject_len) # Original length of the subject sequence
        self.query_start = int(query_start)
        self.query_end = int(query_end)
        self.subject_start = int(subject_start)
        self.subject_end = int(subject_end)

    def __repr__(self):
        return (f"RBHHit(q='{self.query_id}', s='{self.subject_id}', "
                f"e={self.evalue:.2e}, score={self.bitscore})")

    def __lt__(self, other):
        """Sort by bitscore (descending), then by e-value (ascending)."""
        if self.bitscore != other.bitscore:
            return self.bitscore > other.bitscore # Higher bitscore is better
        return self.evalue < other.evalue # Lower e-value is better

class RBHHitGroup:
    """Stores all hits for a given query and helps find the best one(s)."""
    def __init__(self, query_id):
        self.query_id = query_id
        self.hits = []
        self._best_hits_cache = None

    def add_hit(self, hit: RBHHit):
        if hit.query_id != self.query_id:
            raise ValueError("Hit's query_id does not match group's query_id")
        self.hits.append(hit)
        self._best_hits_cache = None # Invalidate cache

    @property
    def best_hits(self):
        """Returns a list of the best hit(s) based on score and e-value.
        Handles ties by returning all hits that are equally best."""
        if not self.hits:
            return []
        if self._best_hits_cache is not None:
            return self._best_hits_cache

        self.hits.sort() # Uses RBHHit.__lt__

        best_score = self.hits[0].bitscore
        best_evalue = self.hits[0].evalue

        self._best_hits_cache = [
            h for h in self.hits
            if h.bitscore == best_score and h.evalue == best_evalue
        ]
        return self._best_hits_cache

    def __repr__(self):
        return f"RBHHitGroup(query='{self.query_id}', num_hits={len(self.hits)})"

# Further functions for BLAST execution, parsing, and RBH logic will be added below.
import os
import uuid
from Bio import SeqIO # For get_sequence_lengths_from_fasta

# (Keep existing logger and class definitions for RBHHit, RBHHitGroup)

def get_sequence_lengths_from_fasta(fasta_file_path):
    """Reads a FASTA file and returns a dictionary of sequence IDs to their lengths."""
    lengths = {}
    try:
        for record in SeqIO.parse(fasta_file_path, "fasta"):
            lengths[str(record.id)] = len(record.seq)
    except FileNotFoundError:
        logger.error(f"FASTA file not found for length extraction: {fasta_file_path}")
    except Exception as e:
        logger.error(f"Error reading FASTA file {fasta_file_path} for lengths: {e}", exc_info=True)
    return lengths

def perform_blast_for_rbh(query_fasta_path, target_fasta_path,
                          blastn_exe, makeblastdb_exe,
                          blast_handler_module, # Pass the actual module
                          temp_dir_base, # Base directory for temp files for this RBH run
                          utils_module): # For hash_file_content if create_or_get_blast_db is used
    """
    Performs BLASTN search of query against target.
    Uses functions from the provided blast_handler_module.
    Returns path to BLAST XML output, or None on failure.
    """
    logger.info(f"RBH BLAST: Querying {query_fasta_path} against {target_fasta_path}")

    if not all([query_fasta_path, target_fasta_path, blastn_exe, makeblastdb_exe, blast_handler_module, temp_dir_base, utils_module]):
        logger.error("Missing one or more critical parameters for perform_blast_for_rbh.")
        return None

    try:
        # Create a specific subdirectory for this BLAST run's temporary files to avoid collisions
        run_temp_dir = os.path.join(temp_dir_base, f"blast_run_{uuid.uuid4().hex}")
        os.makedirs(run_temp_dir, exist_ok=True)

        db_temp_dir = os.path.join(run_temp_dir, "db")
        blast_out_temp_dir = os.path.join(run_temp_dir, "output")
        os.makedirs(db_temp_dir, exist_ok=True)
        os.makedirs(blast_out_temp_dir, exist_ok=True)

        # Create BLAST database from target_fasta_path
        # Using create_or_get_blast_db for potential caching, though for RBH, DBs might be transient
        # The blast_handler.create_or_get_blast_db expects a utils_module for hashing.
        db_prefix = blast_handler_module.create_or_get_blast_db(
            target_fasta_path=target_fasta_path,
            temp_dir=db_temp_dir, # Store this specific DB here
            makeblastdb_executable_path=makeblastdb_exe,
            utils_module=utils_module # Pass the utils module here
        )
        if not db_prefix:
            logger.error(f"Failed to create/get BLAST database for {target_fasta_path}")
            return None

        # Run BLASTN
        blast_xml_output_path = blast_handler_module.run_blastn(
            query_sequence_path=query_fasta_path,
            db_prefix_path=db_prefix,
            blastn_executable_path=blastn_exe,
            output_dir=blast_out_temp_dir # Store BLAST output here
            # config can be passed if specific BLAST params like e-value, word_size are needed here
            # For now, using defaults in blast_handler.run_blastn
        )

        if not blast_xml_output_path:
            logger.error(f"BLASTN run failed for query {query_fasta_path} against DB {db_prefix}")
            return None

        logger.info(f"BLASTN for RBH completed. Output: {blast_xml_output_path}")
        # Note: The run_temp_dir and its contents (DB, output) should be cleaned up by the caller eventually.
        return blast_xml_output_path

    except Exception as e:
        logger.error(f"Exception in perform_blast_for_rbh: {e}", exc_info=True)
        return None

def parse_blast_xml(blast_xml_path, query_seq_lengths, subject_seq_lengths):
    """
    Parses BLAST XML output and returns a dictionary of query IDs to RBHHitGroup objects.
    query_seq_lengths: dict mapping query_id to its length.
    subject_seq_lengths: dict mapping subject_id to its length.
    """
    hits_by_query = {}
    if not blast_xml_path or not os.path.exists(blast_xml_path):
        logger.error(f"BLAST XML file not found: {blast_xml_path}")
        return hits_by_query

    try:
        with open(blast_xml_path) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for record in blast_records:
                query_id = str(record.query_id) # Ensure string type, consistent with FASTA parsing

                if query_id not in hits_by_query:
                    hits_by_query[query_id] = RBHHitGroup(query_id)

                q_len = query_seq_lengths.get(query_id, record.query_length) # Fallback to record's query_length if not in map

                for alignment in record.alignments:
                    subject_id = str(alignment.hit_id) # Ensure string type
                    s_len = subject_seq_lengths.get(subject_id, alignment.length) # Fallback to alignment.length

                    for hsp in alignment.hsps:
                        hit = RBHHit(
                            query_id=query_id,
                            subject_id=subject_id,
                            pident=(hsp.identities * 100.0 / hsp.align_length) if hsp.align_length > 0 else 0.0,
                            length=hsp.align_length,
                            evalue=hsp.expect,
                            bitscore=hsp.score,
                            query_len=q_len,
                            subject_len=s_len,
                            query_start=hsp.query_start,
                            query_end=hsp.query_end,
                            subject_start=hsp.sbjct_start,
                            subject_end=hsp.sbjct_end
                        )
                        hits_by_query[query_id].add_hit(hit)
        logger.info(f"Parsed {len(hits_by_query)} queries from {blast_xml_path}")
    except Exception as e:
        logger.error(f"Error parsing BLAST XML {blast_xml_path}: {e}", exc_info=True)

    return hits_by_query

# Core RBH logic function will be added next
# (Keep existing logger, classes, and functions: get_sequence_lengths_from_fasta, perform_blast_for_rbh, parse_blast_xml)

def find_rbh_pairs(hits_A_vs_B, hits_B_vs_A, strict_reciprocity=True):
    """
    Finds Reciprocal Best Hit (RBH) pairs from two sets of BLAST results.

    Args:
        hits_A_vs_B (dict): {query_A_id: RBHHitGroup} from BLAST of A (queries) vs B (subjects).
        hits_B_vs_A (dict): {query_B_id: RBHHitGroup} from BLAST of B (queries) vs A (subjects).
        strict_reciprocity (bool): If True, requires that for a pair (a,b) to be RBH,
                                   'a' must be among b's best hits AND 'b' must be among a's best hits,
                                   considering all tied best hits. If False, a more relaxed
                                   definition might be used (e.g., if any of a's best hits b_i has a
                                   as a best hit, even if not all b_i do).
                                   For now, this implementation implies strict handling of ties:
                                   A is best hit for B if A is in B's best_hits list.

    Returns:
        list: A list of tuples, where each tuple is (id_A, id_B, rbh_hit_A_to_B, rbh_hit_B_to_A).
              rbh_hit_A_to_B is one of the RBHHit objects from A to B that qualifies.
              rbh_hit_B_to_A is one of the RBHHit objects from B to A that qualifies.
    """
    rbh_pairs = []
    processed_A_ids = set() # To avoid duplicate pairs if A has multiple best hits to B that are RBH

    if not hits_A_vs_B or not hits_B_vs_A:
        logger.warning("One or both BLAST hit dictionaries are empty. Cannot find RBH pairs.")
        return rbh_pairs

    for query_A_id, group_A_vs_B in hits_A_vs_B.items():
        if query_A_id in processed_A_ids:
            continue

        best_hits_A_to_B = group_A_vs_B.best_hits
        if not best_hits_A_to_B:
            logger.debug(f"Query {query_A_id} has no hits in B. Skipping.")
            continue

        # Iterate through all best hits of A in B (handles ties)
        for hit_A_to_B in best_hits_A_to_B:
            subject_B_id = hit_A_to_B.subject_id # This is the ID from genome B

            # Now, check if this subject_B_id has query_A_id as a best hit in return
            group_B_vs_A = hits_B_vs_A.get(subject_B_id)
            if not group_B_vs_A:
                logger.debug(f"Subject {subject_B_id} (hit from {query_A_id}) not found as a query in B_vs_A hits. Skipping.")
                continue

            best_hits_B_to_A = group_B_vs_A.best_hits
            if not best_hits_B_to_A:
                logger.debug(f"Subject {subject_B_id} (hit from {query_A_id}) has no hits back in A. Skipping.")
                continue

            # Check if our original query_A_id is among the best hits for subject_B_id
            is_reciprocal_best = False
            qualifying_hit_B_to_A = None
            for hit_B_to_A in best_hits_B_to_A:
                if hit_B_to_A.subject_id == query_A_id:
                    is_reciprocal_best = True
                    qualifying_hit_B_to_A = hit_B_to_A
                    break # Found reciprocity for this specific hit_A_to_B

            if is_reciprocal_best:
                logger.info(f"RBH found: {query_A_id} <-> {subject_B_id}")
                rbh_pairs.append((query_A_id, subject_B_id, hit_A_to_B, qualifying_hit_B_to_A))
                # Mark query_A_id as processed to avoid adding it multiple times
                # if it has other best hits in B that also form RBHs.
                # This ensures one entry per A-ID in the output list, based on the first B-ID that forms an RBH.
                # If multiple B-IDs form RBH with A1, all will be added if hit_A_to_B iteration continues.
                # The current logic: if A1->B1 is best, and A1->B2 is also best (tie),
                # and B1->A1 is best, and B2->A1 is best, then (A1,B1) and (A1,B2) are both RBHs.
                # This seems correct. `processed_A_ids` might not be strictly necessary with this loop structure.
                # Let's remove processed_A_ids for now to allow one-to-many RBHs if ties occur on both sides.
                # Example: A1 best hits are B1, B2. B1 best hit is A1. B2 best hit is A1. -> (A1,B1), (A1,B2)

    # Remove duplicate pairs that might arise if, e.g. (A1, B1) is found, and later when processing B1 as a query (if inputs were symmetric)
    # However, the input dictionaries are hits_A_vs_B and hits_B_vs_A, so we only iterate A's queries.
    # Duplicates like (A1,B1) and (B1,A1) are not possible unless inputs are swapped.
    # What if A1->B1, A1->B2 are best hits, and B1->A1, B2->A1 are best hits.
    # Loop for A1:
    #   hit_A_to_B = (A1,B1,...) -> check B1. B1 best hits include A1. Add (A1,B1).
    #   hit_A_to_B = (A1,B2,...) -> check B2. B2 best hits include A1. Add (A1,B2).
    # This is the desired behavior for one-to-many or many-to-many RBHs due to ties.

    logger.info(f"Found {len(rbh_pairs)} RBH pairs.")
    return rbh_pairs
